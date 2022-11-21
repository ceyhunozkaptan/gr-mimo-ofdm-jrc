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
#include "target_simulator_impl.h"
#include <volk/volk.h>
#include <gnuradio/math.h>


namespace gr {
namespace mimo_ofdm_jrc {
    const double FOUR_PI_CUBED_SQRT = 44.54662397465366;

    target_simulator::sptr
    target_simulator::make(std::vector<float> range,
                                 std::vector<float> velocity,
                                 std::vector<float> rcs,
                                 std::vector<float> azimuth,
                                 std::vector<float> position_rx,
                                 int samp_rate,
                                 float center_freq,
                                 float self_coupling_db,
                                 bool rndm_phaseshift,
                                 bool self_coupling,
                                 const std::string& len_key,
                                 bool debug)
    {
      return gnuradio::get_initial_sptr
        (new target_simulator_impl(range,
                                    velocity,
                                    rcs,
                                    azimuth,
                                    position_rx,
                                    samp_rate,
                                    center_freq,
                                    self_coupling_db,
                                    rndm_phaseshift,
                                    self_coupling,
                                    len_key,
                                    debug));
    }


    /*
     * The private constructor
     */
    target_simulator_impl::target_simulator_impl(std::vector<float> range,
                                                std::vector<float> velocity,
                                                std::vector<float> rcs,
                                                std::vector<float> azimuth,
                                                std::vector<float> position_rx,
                                                int samp_rate,
                                                float center_freq,
                                                float self_coupling_db,
                                                bool rndm_phaseshift,
                                                bool self_coupling,
                                                const std::string& len_key,
                                                bool debug)
      : gr::tagged_stream_block("target_simulator",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(position_rx.size(), position_rx.size(), sizeof(gr_complex)), len_key)
    {
        d_debug = debug;

        d_buff_size = 2;

        d_fft = new fft::fft_complex(d_buff_size, true);
        d_ifft = new fft::fft_complex(d_buff_size, false);

        d_buff_time = (gr_complex*)volk_malloc(sizeof(gr_complex) * d_buff_size, volk_get_alignment());
        d_buff_freq = (gr_complex*)volk_malloc(sizeof(gr_complex) * d_buff_size, volk_get_alignment());
        memset(d_buff_time, 0, d_buff_size * sizeof(gr_complex));
        memset(d_buff_freq, 0, d_buff_size * sizeof(gr_complex));

        setup_targets(range,
                velocity,
                rcs,
                azimuth,
                position_rx,
                samp_rate,
                center_freq,
                self_coupling_db,
                rndm_phaseshift,
                self_coupling);
    }

    /*
     * Our virtual destructor.
     */
    target_simulator_impl::~target_simulator_impl()
    {
        volk_free(d_buff_time);
        volk_free(d_buff_freq);

        delete d_fft;
        delete d_ifft;
    }

    int
    target_simulator_impl::calculate_output_stream_length(const gr_vector_int &ninput_items)
    {
      int noutput_items = ninput_items[0];
      return noutput_items ;
    }

    void target_simulator_impl::setup_targets(std::vector<float> range,
                                                    std::vector<float> velocity,
                                                    std::vector<float> rcs,
                                                    std::vector<float> azimuth,
                                                    std::vector<float> position_rx,
                                                    int samp_rate,
                                                    float center_freq,
                                                    float self_coupling_db,
                                                    bool rndm_phaseshift,
                                                    bool self_coupling)
    {
        gr::thread::scoped_lock lock(d_setlock);

        d_range = range;
        d_velocity = velocity;
        d_rcs = rcs;
        d_azimuth = azimuth;
        d_position_rx = position_rx;
        d_center_freq = center_freq; // center frequency of simulated hardware for doppler estimation
        d_samp_rate = samp_rate;
        d_rndm_phaseshift = rndm_phaseshift;
        d_self_coupling = self_coupling;
        d_self_coupling_db = self_coupling_db;
        
        d_new_channel = true;

        // Setup rx_time tag
        d_key = pmt::string_to_symbol("rx_time");
        d_srcid = pmt::string_to_symbol("stat_targ_sim");

        // Get num targets
        d_num_targets = range.size(); // FIXME: throw exceptions for len(range)!=len(velocity)!=...

        // Get doppler frequencies
        d_doppler.resize(d_num_targets);
        d_filt_doppler.resize(d_num_targets);
        for (int k = 0; k < d_num_targets; k++){
            d_doppler[k] = 2 * d_velocity[k] * d_center_freq / c_light;

        }

        // Get timeshifts
        d_timeshift_azimuth.resize(d_position_rx.size());
        d_filt_time_azimuth.resize(d_position_rx.size());
        for (int l = 0; l < d_position_rx.size(); l++) {
            
            d_filt_time_azimuth[l].resize(d_num_targets);
            d_timeshift_azimuth[l].resize(d_num_targets);
            
            for (int k = 0; k < d_num_targets; k++) {
                d_timeshift_azimuth[l][k] = (2.0 * d_range[k] - d_position_rx[l] * std::sin(d_azimuth[k] * GR_M_PI / 180.0)) / c_light;
            }
        }

        
        // Get signal amplitude of reflection with free space path loss and rcs (radar
        // equation)
        d_scale_ampl.resize(d_num_targets);
        for (int k = 0; k < d_num_targets; k++) 
        {
            // Factor out all terms out of the sqrt except the RCS:
            d_scale_ampl[k] = c_light * std::sqrt(d_rcs[k]) / FOUR_PI_CUBED_SQRT / (d_range[k] * d_range[k]) / d_center_freq;
        }
        
        if (d_rndm_phaseshift) 
        {
            // Resize phase shift filter
            d_filt_phase.resize(d_num_targets);

            // Setup random numbers
            std::srand(std::time(NULL)); // initial with time
        }
    }

    int
    target_simulator_impl::work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
        gr::thread::scoped_lock lock(d_setlock);
        dout << "================================================================================" << std::endl;        

        const gr_complex* in = (const gr_complex*)input_items[0];
        gr_complex* out;
        
        dout << "[TARGET SIM] ninput_items[0] = " << ninput_items[0] << ", noutput_items: " << noutput_items << std::endl;

        // Set output items to tagged stream length
        // noutput_items = ninput_items[0];

        int n_input = ninput_items[0];

        dout << "[TARGET SIM] Current buffer size: " << d_buff_size << ", n_input: " << n_input << std::endl;

        if (d_buff_size != n_input) 
        {
            volk_free(d_buff_time);
            volk_free(d_buff_freq);

            // Setup buffer length
            d_buff_time = (gr_complex*)volk_malloc(sizeof(gr_complex) * n_input, volk_get_alignment());
            d_buff_freq = (gr_complex*)volk_malloc(sizeof(gr_complex) * n_input, volk_get_alignment());
            memset(d_buff_time, 0, n_input * sizeof(gr_complex));
            memset(d_buff_freq, 0, n_input * sizeof(gr_complex));

            if (d_buff_time == 0 || d_buff_freq == 0) 
            {
                volk_free(d_buff_freq);
                volk_free(d_buff_time);
                throw std::runtime_error("[TARGET SIM] volk_malloc");
            }

            dout << "[TARGET SIM] Setup buffers " << std::endl;

            delete d_fft;
            d_fft = new fft::fft_complex(n_input, true);
            delete d_ifft;
            d_ifft = new fft::fft_complex(n_input, false);

            dout << "[TARGET SIM] Setup FFTs " << std::endl;
        }
        
        // Check if new filter, buffer or fft plan is necessary
        if (d_buff_size != n_input || d_new_channel) 
        {
            dout << "[TARGET SIM] New channel response computed " << std::endl;

            // Setup frequency vector for shift in frequency domain
            if (d_buff_size != n_input)
            {
                d_freq.resize(n_input);
            }

            for (int i = 0; i < n_input; i++) 
            {
                if (i < n_input / 2)
                    d_freq[i] = i * (float)d_samp_rate / (float)n_input; // zero to samp_rate/2
                else
                    d_freq[i] = i * (float)d_samp_rate / (float)n_input - (float)d_samp_rate; // -samp_rate/2 to zero
            }

            // Setup freq and time shift filter, resize phase shift filter
            for (int k = 0; k < d_num_targets; k++) 
            {
                if (d_buff_size != n_input)
                {
                    d_filt_doppler[k].resize(n_input);
                    if (d_rndm_phaseshift){
                        d_filt_phase[k].resize(n_input);
                    }
                }

                d_phase_doppler = 0;
                for (int i = 0; i < n_input; i++) 
                {
                    // Doppler shift filter and rescaling amplitude with rcs
                    d_filt_doppler[k][i] = std::exp(d_phase_doppler) * d_scale_ampl[k];
                    d_phase_doppler = gr_complex(0, std::fmod(std::imag(d_phase_doppler) + 2 * GR_M_PI * d_doppler[k] / (float)d_samp_rate, 2 * GR_M_PI)); // integrate phase (with plus!)
                }

                // Do time shift filter with azimuth and position, there are two
                // time shift filters to avoid problems with significant digits of float
                for (int l = 0; l < d_position_rx.size(); l++) 
                {
                    if (d_buff_size != n_input)
                    {
                        d_filt_time_azimuth[l][k].resize(n_input);
                    }
                    d_phase_time = 0;
                    for (int i = 0; i < n_input; i++) 
                    {
                        // Time shift filter, uses azimuth and RX position
                        d_phase_time = gr_complex(0, std::fmod(2 * GR_M_PI *(d_timeshift_azimuth[l][k]) * (d_freq[i] + d_center_freq), 2 * GR_M_PI)); // integrate phase (with minus!)
                        d_filt_time_azimuth[l][k][i] = std::exp(-d_phase_time) / (float)n_input; 
                        // dout << "[TARGET SIM] d_phase_time = " << d_phase_time << std::endl;
                    }
                }
            }
            // Resize hold of n_input
            d_new_channel = false;
            dout << "[TARGET SIM] Resize buffer size to " << n_input << " from " << d_buff_size << std::endl;
            d_buff_size = n_input;
        }

        // Setup random phase shift
        if (d_rndm_phaseshift) 
        {
            gr_complex phase_random_hold;
            for (int k = 0; k < d_num_targets; k++) 
            {
                phase_random_hold = gr_complex(0, 2 * GR_M_PI * float((std::rand() % 1000 + 1) / 1000.0));
                d_phase_random = std::exp(phase_random_hold);
                std::fill_n(&d_filt_phase[k][0], n_input, d_phase_random);
            }
        }
        
        dout << "[TARGET SIM] Computing for all RX antennas " << std::endl;

        // Go through RXs
        for (int l = 0; l < d_position_rx.size(); l++) 
        {
            // Setup pointer on output buffer
            out = (gr_complex*) output_items[l];

            // Set rx_time tag
            d_time_sec = nitems_written(l) / d_samp_rate;
            d_time_frac_sec = nitems_written(l) / (float)d_samp_rate - d_time_sec;
            d_val = pmt::make_tuple(pmt::from_uint64(d_time_sec), pmt::from_double(d_time_frac_sec)); // FIXME: correct implementation?
            add_item_tag(l, nitems_written(l), d_key, d_val, d_srcid);

            // Set output to zero
            std::memset(out, 0, n_input * sizeof(gr_complex));

            // Go through targets and apply filters
            for (int k = 0; k < d_num_targets; k++) 
            {
                
                // Add doppler shift
                volk_32fc_x2_multiply_32fc(&d_buff_time[0], in, &d_filt_doppler[k][0], n_input); // add doppler shift with rescaled amplitude

                // Go to freq domain
                memcpy(d_fft->get_inbuf(), &d_buff_time[0], sizeof(gr_complex)*n_input);
                d_fft->execute();

                // Add timeshift with multiply exp-func in freq domain (rx position with azimuth and range)
                volk_32fc_x2_multiply_32fc(&d_buff_freq[0], d_fft->get_outbuf(), &d_filt_time_azimuth[l][k][0], n_input); 

                // Go back to time domain
                memcpy(d_ifft->get_inbuf(), d_buff_freq, sizeof(gr_complex)*n_input);
                d_ifft->execute();

                if (d_rndm_phaseshift) 
                {
                    // Add random phase shift
                    volk_32fc_x2_multiply_32fc(out, d_ifft->get_outbuf(), &d_filt_phase[k][0], n_input);
                }
                else
                {
                    memcpy(out, d_ifft->get_outbuf(), sizeof(gr_complex) * n_input);
                }

            }

            // Add self coupling
            if (d_self_coupling) 
            {
                for (int i = 0; i < n_input; i++)
                {
                    out[i] += (gr_complex)pow(10, d_self_coupling_db / 20.0) * in[i]; // d_self_coupling_db gives scaling of power
                }
            }
        }

        dout << "[TARGET SIM] DONE! " << std::endl;

        // Tell runtime system how many output items we produced.
        return n_input;
    }

  } /* namespace mimo_ofdm_jrc */
} /* namespace gr */

