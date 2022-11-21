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
#include "fft_peak_detect_impl.h"

namespace gr {
  namespace mimo_ofdm_jrc {

    fft_peak_detect::sptr
    fft_peak_detect::make(int samp_rate, float interp_factor, float threshold, int samp_protect, std::vector<float> max_freq, bool cut_max_freq,
                                            const std::string& len_key)
    {
      return gnuradio::get_initial_sptr
        (new fft_peak_detect_impl(samp_rate, interp_factor, threshold, samp_protect, max_freq, cut_max_freq, len_key));
    }


    /*
     * The private constructor
     */
    fft_peak_detect_impl::fft_peak_detect_impl(int samp_rate, float interp_factor, float threshold, int samp_protect, std::vector<float> max_freq, bool cut_max_freq,
                                            const std::string& len_key)
      : gr::tagged_stream_block("fft_peak_detect",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make3(3, 3, sizeof(float), sizeof(float), sizeof(float)),
                              len_key)
    {
    	d_samp_rate = samp_rate;
        d_interp_factor = interp_factor;
        d_threshold = threshold;
        d_samp_protect = samp_protect;
        d_max_freq = max_freq;
        d_cut_max_freq = cut_max_freq;

        set_tag_propagation_policy(TPP_DONT);
    }

    /*
     * Our virtual destructor.
     */
    fft_peak_detect_impl::~fft_peak_detect_impl()
    {
    }

    int
    fft_peak_detect_impl::calculate_output_stream_length(const gr_vector_int &ninput_items)
    {
      return 1 ;
    }

    void fft_peak_detect_impl::set_threshold(float threshold) { d_threshold = threshold; }
    void fft_peak_detect_impl::set_samp_protect(int samp) { d_samp_protect = samp; }
    void fft_peak_detect_impl::set_max_freq(std::vector<float> freq) { d_max_freq = freq; }

    int fft_peak_detect_impl::work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
        const gr_complex* in = (const gr_complex*)input_items[0];

        float* out_freq = (float*)output_items[0];
        float* out_phase = (float*)output_items[1];
        float* out_mag = (float*)output_items[2];

        int k = -1;
        float hold = -1;
        for (int p = d_samp_protect; p < ninput_items[0] - d_samp_protect; p++) { // implementation of protected samples
            if (std::abs(in[p]) > hold && std::pow(std::abs(in[p]), 2) > std::pow(10, d_threshold / 10.0)) {
                hold = std::abs(in[p]);
                k = p;
            }
        }


        if (k != -1) {
            if (k <= ninput_items[0] / 2){
                out_freq[0] = k/(float)ninput_items[0] * (d_samp_rate * d_interp_factor);
            }
            else{
                out_freq[0] = -((float)d_samp_rate*d_interp_factor) + k * (d_samp_rate * d_interp_factor / (float)ninput_items[0]);
            }
            out_phase[0] = std::arg(in[k]);
            out_mag[0] = std::abs(in[k]);
        }

        // Tell runtime system how many output items we produced.
        return 1;
     }

  } /* namespace mimo_ofdm_jrc */
} /* namespace gr */

