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

#ifndef INCLUDED_MIMO_OFDM_JRC_OFDM_FRAME_GENERATOR_IMPL_H
#define INCLUDED_MIMO_OFDM_JRC_OFDM_FRAME_GENERATOR_IMPL_H

#include <mimo_ofdm_jrc/ofdm_frame_generator.h>

namespace gr {
  namespace mimo_ofdm_jrc {

    class ofdm_frame_generator_impl : public ofdm_frame_generator
    {
     private:
      	//! FFT length
        const int d_fft_len;
        //! Which carriers/symbols carry data
        std::vector<std::vector<int>> d_occupied_carriers;
        //! Which carriers/symbols carry pilots symbols
        std::vector<std::vector<int>> d_pilot_carriers;
        //! Value of said pilot symbols
        const std::vector<std::vector<gr_complex>> d_pilot_symbols;
        //! Synch words
        const std::vector<std::vector<gr_complex>> d_sync_words;

        const int d_ltf_len;
        int d_symbols_per_set;
        const bool d_output_is_shifted;

     protected:
      int calculate_output_stream_length(const gr_vector_int &ninput_items);

     public:
      ofdm_frame_generator_impl(int fft_len,
                                const std::vector<std::vector<int>>& occupied_carriers,
                                const std::vector<std::vector<int>>& pilot_carriers,
                                const std::vector<std::vector<gr_complex>>& pilot_symbols,
                                const std::vector<std::vector<gr_complex>>& sync_words,
                                int ltf_len,
                                const std::string& len_tag_key,
                                const bool output_is_shifted);
      ~ofdm_frame_generator_impl();

      // Where all the action really happens
      int work(
              int noutput_items,
              gr_vector_int &ninput_items,
              gr_vector_const_void_star &input_items,
              gr_vector_void_star &output_items
      );
      
        std::string len_tag_key() { return d_length_tag_key_str; };

        const int fft_len() { return d_fft_len; };
        std::vector<std::vector<int>> occupied_carriers() { return d_occupied_carriers; };
    };

  } // namespace mimo_ofdm_jrc
} // namespace gr

#endif /* INCLUDED_MIMO_OFDM_JRC_OFDM_FRAME_GENERATOR_IMPL_H */

