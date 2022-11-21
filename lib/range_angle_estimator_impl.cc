/* -*- c++ -*- */
/*
 * Copyright 2022 Ceyhun D. Ozkaptan @ The Ohio State University <ozkaptan.1@osu.edu>
 *
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include "range_angle_estimator_impl.h"

namespace gr {
  namespace mimo_ofdm_jrc {

    range_angle_estimator::sptr
    range_angle_estimator::make(int vlen, 
                        std::vector<float> range_bins,
                        std::vector<float> angle_bins,
                        float noise_discard_range_m,
                        float noise_discard_angle_deg,
                        float snr_threshold, 
                        float power_threshold,
                        const std::string& stats_path, 
                        bool stats_record, 
                        const std::string& len_key,
                        bool debug)
    {
      return gnuradio::get_initial_sptr
        (new range_angle_estimator_impl(vlen, 
                                        range_bins,
                                        angle_bins,
                                        noise_discard_range_m,
                                        noise_discard_angle_deg,
                                        snr_threshold, 
                                        power_threshold, 
                                        stats_path,
                                        stats_record, 
                                        len_key,
                                        debug));
    }


    /*
     * The private constructor
     */
    range_angle_estimator_impl::range_angle_estimator_impl(int vlen, 
                                                            std::vector<float> range_bins,
                                                            std::vector<float> angle_bins,
                                                            float noise_discard_range_m,
                                                            float noise_discard_angle_deg,
                                                            float snr_threshold, 
                                                            float power_threshold, 
                                                            const std::string& stats_path,
                                                            bool stats_record, 
                                                            const std::string& len_key,
                                                            bool debug)
      : gr::tagged_stream_block("range_angle_estimator",
                gr::io_signature::make(1, 1, sizeof(gr_complex)*vlen),
                gr::io_signature::make(0, 0, 0), 
                len_key),
                d_vlen(vlen),
                d_range_bins(range_bins),
                d_angle_bins(angle_bins),
                d_noise_discard_range_m(noise_discard_range_m),
                d_noise_discard_angle_deg(noise_discard_angle_deg),
                d_snr_threshold(snr_threshold), 
                d_power_threshold(power_threshold), 
                d_stats_path(stats_path),
                d_stats_record(stats_record), 
                d_debug(debug),
                file_stream(d_stats_path, std::ofstream::app)
    {
        message_port_register_out(pmt::mp("params"));
        // TODO sanity checks


        if (!file_stream.is_open())
        {
            std::cerr << "[RANGE-ANGLE ESTIMATOR] Could not open log file at " << d_stats_path << std::endl;
        }
        file_stream.close();

        d_new_stat_started = false;
    }

    /*
     * Our virtual destructor.
     */
    range_angle_estimator_impl::~range_angle_estimator_impl()
    {
        if (file_stream.is_open())
        {
            file_stream.flush();
            file_stream.close();
        }
    }

    int
    range_angle_estimator_impl::calculate_output_stream_length(const gr_vector_int &ninput_items)
    {
        int noutput_items = 0;
        return noutput_items;
    }

    int
    range_angle_estimator_impl::work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
        const gr_complex* in = (const gr_complex*)input_items[0];

        // Get mean power
        int n_inputs = ninput_items[0]; // size of i_range axis
        float mean_power = 0;
        float peak_power = -1;
        float curr_power;
        int peak_range_idx = -1;
        int peak_angle_idx = -1;

        for (int i_range = 0; i_range < n_inputs; i_range++)  // go through range axis
        {
            for (int i_angle = 0; i_angle < d_vlen; i_angle++) // go through angle axis
            {   
                curr_power = std::pow(std::abs( in[i_angle + d_vlen * i_range] ), 2);
                mean_power += curr_power;

                if(curr_power > peak_power)
                {
                    peak_power = curr_power;
                    peak_range_idx = i_range;
                    peak_angle_idx = i_angle; 
                }
            }
        }
        float angle_val = d_angle_bins[peak_angle_idx];
        float range_val = d_range_bins[peak_range_idx];

        float angle_null = angle_val + 90;

        if (angle_null >= 90)
        {
            angle_null = angle_null - 180;
        }

        int angle_null_idx;
        auto iter_geq = std::lower_bound(
            d_angle_bins.begin(), 
            d_angle_bins.end(), 
            angle_null
        );

        double a = *(iter_geq - 1);
        double b = *(iter_geq);

        if (iter_geq == d_angle_bins.begin()) {
            angle_null_idx = 0;
        }
        else if (fabs(angle_null - a) < fabs(angle_null - b)) {
            angle_null_idx = iter_geq - d_angle_bins.begin() - 1;
        }
        else{
            angle_null_idx = iter_geq - d_angle_bins.begin();
        }

        dout << "[RANGE-ANGLE ESTIMATOR] angle_null: " << d_angle_bins[angle_null_idx]  << std::endl;

        if ( angle_null_idx == d_angle_bins.size()-1 )
        {
            angle_null_idx = d_angle_bins.size()-2;
        }

        int discard_range_idx = d_noise_discard_range_m / (d_range_bins[1]-d_range_bins[0]);
        int discard_angle_idx = d_noise_discard_angle_deg / (d_angle_bins[(angle_null_idx+1)%d_angle_bins.size()]-d_angle_bins[angle_null_idx]);

        if (discard_angle_idx <= 0)
        {
            discard_angle_idx = 1;
        }

        int start_range_idx = (peak_range_idx+d_range_bins.size()/2-discard_range_idx);
        int end_range_idx = (peak_range_idx+d_range_bins.size()/2+discard_range_idx);

        int start_angle_idx = angle_null_idx - discard_angle_idx;
        int end_angle_idx = angle_null_idx + discard_angle_idx;

        dout << "[RANGE-ANGLE ESTIMATOR] start_range_idx: " << start_range_idx << ", end_range_idx:" << end_range_idx  << std::endl;
        dout << "[RANGE-ANGLE ESTIMATOR] start_angle_idx: " << start_angle_idx << ", end_angle_idx:" << end_angle_idx  << std::endl;

        dout << "[RANGE-ANGLE ESTIMATOR] start_range: " << d_range_bins[((start_range_idx % n_inputs) + n_inputs)%n_inputs ] << ", end_range:" << d_range_bins[((end_range_idx % n_inputs) + n_inputs)%n_inputs] << std::endl;
        dout << "[RANGE-ANGLE ESTIMATOR] start_angle: " << d_angle_bins[((start_angle_idx % d_vlen) + d_vlen)%d_vlen] << ", end_angle:" << d_angle_bins[((end_angle_idx % d_vlen) + d_vlen)%d_vlen]  << std::endl;

        float noise_power = 0;
        int n_noise_samples = 0;
        for (int i_range = start_range_idx; i_range < end_range_idx; i_range++)  // go through range axis
        {
            int r_idx = ((i_range % n_inputs) + n_inputs)%n_inputs;
            
            for (int i_angle = start_angle_idx; i_angle < end_angle_idx ; i_angle++) // go through angle axis
            {   
                int a_idx = ((i_angle % d_vlen) + d_vlen)%d_vlen;
                noise_power += std::pow(std::abs( in[a_idx + d_vlen * r_idx] ), 2);;
                n_noise_samples++;
            }
        }

        dout << "[RANGE-ANGLE ESTIMATOR] n_noise_samples: " << n_noise_samples << std::endl;


        noise_power = noise_power/n_noise_samples;
        float snr_est = 10 * std::log10( peak_power / noise_power );

        dout << "[RANGE-ANGLE ESTIMATOR] range: " << d_range_bins[peak_range_idx] << ", angle: " << d_angle_bins[peak_angle_idx] << std::endl;
        dout << "[RANGE-ANGLE ESTIMATOR] peak_power: " << peak_power << ", noise_power: " << noise_power << std::endl;
        dout << "[RANGE-ANGLE ESTIMATOR] snr_est: " << snr_est << std::endl;


        if (snr_est >= d_snr_threshold && peak_power >= d_power_threshold)
        {
            pmt::pmt_t d_range_key = pmt::string_to_symbol("range");          // identifier 
            pmt::pmt_t d_range_value = pmt::init_f32vector(1,&range_val); // pmt
            pmt::pmt_t d_range_pack = pmt::list2(d_range_key, d_range_value); // make list 

            pmt::pmt_t d_angle_key = pmt::string_to_symbol("angle");          // identifier 
            pmt::pmt_t d_angle_value = pmt::init_f32vector(1,&angle_val); // pmt
            pmt::pmt_t d_angle_pack = pmt::list2(d_angle_key, d_angle_value); // make list 
            
            pmt::pmt_t d_power_key = pmt::string_to_symbol("power");          // identifier 
            pmt::pmt_t d_power_value = pmt::init_f32vector(1, &peak_power); // pmt
            pmt::pmt_t d_power_pack = pmt::list2(d_power_key, d_power_value); // make list 
            
            pmt::pmt_t d_snr_key = pmt::string_to_symbol("snr");          // identifier 
            pmt::pmt_t d_snr_value = pmt::init_f32vector(1, &snr_est); // pmt
            pmt::pmt_t d_snr_pack = pmt::list2(d_snr_key, d_snr_value); // make list 

            pmt::pmt_t stats_msg = pmt::list4(d_range_pack, d_angle_pack, d_power_pack, d_snr_pack);
            message_port_pub(pmt::mp("params"), stats_msg); // publish message

            if(d_stats_record)
            {
                file_stream.open(d_stats_path, std::ofstream::app);

                dout << "[RANGE-ANGLE ESTIMATOR] d_stats_path:" << d_stats_path << ", " << file_stream.is_open() << std::endl;
                // file_stream.open(d_stats_path, std::ofstream::app);

                if (file_stream.is_open())
                {
                    if(!d_new_stat_started)
                    {
                        file_stream << "\n NEW RECORD - " << current_date_time() << "\n";
                        d_new_stat_started = true;
                    }
                    file_stream << current_date_time2() << ", \t" << peak_power << ", \t" << snr_est << ", \t" << range_val << ", \t" << angle_val << "\n" ;
                    file_stream.flush();
                    file_stream.close();
                }
                else
                {
                    dout << "[RANGE-ANGLE ESTIMATOR] d_stats_path:" << d_stats_path << ", " << file_stream.is_open() << std::endl;

                    throw std::runtime_error("[STREAM DECODER] Could not open file!!");
                }
            }
        }
 

        return 0;
    }

    void range_angle_estimator_impl::set_snr_threshold(float snr_threshold) { d_snr_threshold = snr_threshold; }
    void range_angle_estimator_impl::set_power_threshold(float power_threshold) { d_power_threshold = power_threshold; }
    
    void range_angle_estimator_impl::set_stats_record(bool stats_record) 
    { 
        d_stats_record = stats_record; 
        d_new_stat_started = false;
            
        if(stats_record)
        {
            dout << "[RANGE-ANGLE ESTIMATOR] Recording Stats --> Enabled" << std::endl;
        }
        else
        {
            dout << "[RANGE-ANGLE ESTIMATOR] Recording Stats --> Disabled" << std::endl;
        }
    }

  } /* namespace mimo_ofdm_jrc */
} /* namespace gr */

