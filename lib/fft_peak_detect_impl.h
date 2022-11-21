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

#ifndef INCLUDED_MIMO_OFDM_JRC_FFT_PEAK_DETECT_IMPL_H
#define INCLUDED_MIMO_OFDM_JRC_FFT_PEAK_DETECT_IMPL_H

#include <mimo_ofdm_jrc/fft_peak_detect.h>

namespace gr {
  namespace mimo_ofdm_jrc {

    class fft_peak_detect_impl : public fft_peak_detect
    {
     private:
      // Nothing to declare in this block.

     protected:
      int calculate_output_stream_length(const gr_vector_int &ninput_items);

     public:
      fft_peak_detect_impl(int samp_rate, float interp_factor, float threshold, int samp_protect, std::vector<float> max_freq, bool cut_max_freq,
                                            const std::string& len_key);
      ~fft_peak_detect_impl();
      
      void set_threshold(float threshold);
      void set_samp_protect(int samp);
      void set_max_freq(std::vector<float> freq);

      int d_samp_rate;
      float d_interp_factor;
      float d_threshold;
      int d_samp_protect;
      std::vector<float> d_max_freq;
      bool d_cut_max_freq;

      // std::vector<float> d_pks, d_freq, d_angle;

      // pmt::pmt_t d_port_id;
      // pmt::pmt_t d_ptimestamp, d_pfreq, d_ppks, d_pangle, d_value;
      std::vector<tag_t> d_tags;

      // Where all the action really happens
      int work(
              int noutput_items,
              gr_vector_int &ninput_items,
              gr_vector_const_void_star &input_items,
              gr_vector_void_star &output_items
      );
    };

  } // namespace mimo_ofdm_jrc
} // namespace gr

#endif /* INCLUDED_MIMO_OFDM_JRC_FFT_PEAK_DETECT_IMPL_H */

