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
#include "frame_detector_impl.h"

namespace gr {
  namespace mimo_ofdm_jrc {

    frame_detector::sptr
    frame_detector::make(int fft_len, int cp_len, double threshold, unsigned int min_n_peaks, unsigned int ignore_gap, bool debug)
    {
      return gnuradio::get_initial_sptr
        (new frame_detector_impl(fft_len, cp_len, threshold, min_n_peaks, ignore_gap, debug));
    }


    /*
     * The private constructor
     */
    frame_detector_impl::frame_detector_impl(int fft_len, int cp_len, double threshold, unsigned int min_n_peaks, unsigned int ignore_gap,  bool debug)
      : gr::block("frame_detector",
              gr::io_signature::make3(3, 3, sizeof(gr_complex), sizeof(gr_complex), sizeof(float)),
              gr::io_signature::make(1, 1, sizeof(gr_complex))),
        curr_state(SEARCH),
        n_peaks(0),
        coarse_cfo_est(0),
        copied_samples(0),
        d_fft_len(fft_len),
        d_cp_len(cp_len),
        d_min_n_peaks(min_n_peaks),
        d_threshold(threshold),
        d_debug(debug),
        d_ignore_gap(ignore_gap),  // PREAMBLE_SYMBOLS*(fft_len+cp_len)
        MAX_PEAK_VALUE(2.0),
        MAX_PEAK_DISTANCE( 2*(fft_len+cp_len) ), // 2 STF symbols
        MAX_SAMPLES( 540*(fft_len+cp_len) )
    {
    	set_tag_propagation_policy(block::TPP_DONT);
    }

    /*
     * Our virtual destructor.
     */
    frame_detector_impl::~frame_detector_impl()
    {
    }
    
    int frame_detector_impl::general_work(int noutput_items,
                                        gr_vector_int &ninput_items,
                                        gr_vector_const_void_star &input_items,
                                        gr_vector_void_star &output_items)
    {
        gr::thread::scoped_lock lock(d_mutex);
        gr::block::set_thread_priority(50);

        const gr_complex *in = (const gr_complex *)input_items[0];
        const gr_complex *in_abs = (const gr_complex *)input_items[1];
        const float *in_cor = (const float *)input_items[2];
        gr_complex *out = (gr_complex *)output_items[0];

        int noutput = noutput_items;
        int ninput = std::min(std::min(ninput_items[0], ninput_items[1]), ninput_items[2]);

        switch (curr_state)
        {

        case SEARCH:
        {
            int n_in;

            for (n_in = 0; n_in < ninput; n_in++)
            {
                if (in_cor[n_in] > d_threshold && in_cor[n_in] < MAX_PEAK_VALUE)
                {
                    if (n_peaks < d_min_n_peaks)
                    {
                        n_peaks++;

                        if (n_peaks == 1){
                            first_peak_ind = nitems_read(0) + n_in;
                            dout << "[FRAME DETECTOR] first_peak_ind : " << first_peak_ind << std::endl;
                        }
                    }
                    else if ((nitems_read(0) + n_in - first_peak_ind) < MAX_PEAK_DISTANCE)
                    {
                        dout << "[FRAME DETECTOR] Frame detected --> Max distance between peaks: " << (nitems_read(0) + n_in - first_peak_ind) << std::endl;
                        dout << "[FRAME DETECTOR] Goes to COPY" << std::endl;

                        curr_state = COPY;
                        copied_samples = 0;
                        coarse_cfo_est = arg(in_abs[n_in]) / (d_fft_len/4.0);
                        n_peaks = 0;
                        first_peak_ind = 0;
                        insert_tag(nitems_written(0), coarse_cfo_est, nitems_read(0) + n_in); //-d_min_n_peaks
                        break;
                    }
                    else{
                        n_peaks = 0;
                        first_peak_ind = 0;
                        dout << "[FRAME DETECTOR] Detected peaks are not consecutive --> No frame detected!" << std::endl;
                    }
                }
                else if ((nitems_read(0) + n_in - first_peak_ind) > MAX_PEAK_DISTANCE)
                {
                    n_peaks = 0;
                    first_peak_ind = 0;
                }
            }

            consume_each(n_in);
            return 0;
        }

        case COPY:
        {
            int n_out = 0;
            while (n_out < ninput && n_out < noutput && copied_samples < MAX_SAMPLES)
            {
                if (in_cor[n_out] > d_threshold && in_cor[n_out] < MAX_PEAK_VALUE)
                {
                    if (n_peaks < d_min_n_peaks)
                    {
                        n_peaks++;
                        if (n_peaks == 1){
                            first_peak_ind = nitems_read(0) + n_out;
                            dout << "[FRAME DETECTOR] first_peak_ind : " << first_peak_ind << ", copied_samples: " << copied_samples << std::endl;
                        }
                    }
                    else if((nitems_read(0) + n_out - first_peak_ind) < MAX_PEAK_DISTANCE)
                    {
                        if (copied_samples > d_ignore_gap)
                        {
                            // New frame detected before MAX_SAMPLES search ends
                            dout << "[FRAME DETECTOR] Another Frame detected during COPY!!" << std::endl;
                            dout << "[FRAME DETECTOR] Max distance between peaks: " << (nitems_read(0) + n_out - first_peak_ind) << std::endl;
                            copied_samples = 0;
                            n_peaks = 0;
                            first_peak_ind = 0;
                            coarse_cfo_est = arg(in_abs[n_out]) / (d_fft_len/4.0); // Can estimate normalized CFO of +-2/fft_len
                            insert_tag(nitems_written(0) + n_out, coarse_cfo_est, nitems_read(0) + n_out); //-d_min_n_peaks
                            break;
                        }
                    }
                    else
                    {
                        dout << "[FRAME DETECTOR] Repeating consecutive peak search in COPY" << std::endl;
                        n_peaks = 0;
                        first_peak_ind = 0;
                    }
                }
                else if ((nitems_read(0) + n_out - first_peak_ind) > MAX_PEAK_DISTANCE)
                {
                    n_peaks = 0;
                    first_peak_ind = 0;
                }

                out[n_out] = in[n_out] * exp(gr_complex(0, -coarse_cfo_est * copied_samples));
                n_out++;
                copied_samples++;
            }

            if (copied_samples == MAX_SAMPLES)
            {
                curr_state = SEARCH;
                dout << "[FRAME DETECTOR] No frame detected --> Going back to SEARCH..." << std::endl;
            }

            // dout << "[FRAME DETECTOR] Copied "<< n_out <<  " samples" << std::endl;

            // std::cout << "[FRAME DETECT] Total Work Time: " <<  pc_work_time_total() << std::endl;
            // std::cout << "[FRAME DETECT] Input Buffer: " << pc_input_buffers_full(0) << ",  Output Buffer: "  << pc_output_buffers_full(0) << std::endl;
            consume_each(n_out);
            return n_out;
        }
        }

        throw std::runtime_error("[FRAME DETECTOR] Unknown state");
        return 0;
    }

    void frame_detector_impl::insert_tag(uint64_t item, double cfo_est, uint64_t input_item)
    {
        const pmt::pmt_t key = pmt::string_to_symbol("frame_start");
        const pmt::pmt_t value = pmt::from_double(cfo_est);
        const pmt::pmt_t srcid = pmt::string_to_symbol(name());
        add_item_tag(0, item, key, value, srcid);
    }

  } /* namespace mimo_ofdm_jrc */
} /* namespace gr */

