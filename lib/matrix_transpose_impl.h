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

#ifndef INCLUDED_MIMO_OFDM_JRC_MATRIX_TRANSPOSE_IMPL_H
#define INCLUDED_MIMO_OFDM_JRC_MATRIX_TRANSPOSE_IMPL_H

#include <mimo_ofdm_jrc/matrix_transpose.h>
#include "utils.h"

namespace gr {
  namespace mimo_ofdm_jrc {

    class matrix_transpose_impl : public matrix_transpose
    {
     private:
      	int d_input_len, d_output_len, d_interp_factor;
        std::vector< tag_t > d_tags;
        bool d_debug;

     protected:
      int calculate_output_stream_length(const gr_vector_int &ninput_items);

     public:
      matrix_transpose_impl(int input_len, int output_len, int interp_factor, bool debug, std::string len_key);
      ~matrix_transpose_impl();

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

#endif /* INCLUDED_MIMO_OFDM_JRC_MATRIX_TRANSPOSE_IMPL_H */

