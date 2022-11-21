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

#ifndef INCLUDED_MIMO_OFDM_JRC_ZERO_PAD_IMPL_H
#define INCLUDED_MIMO_OFDM_JRC_ZERO_PAD_IMPL_H

#include <mimo_ofdm_jrc/zero_pad.h>

namespace gr {
  namespace mimo_ofdm_jrc {

    class zero_pad_impl : public zero_pad
    {
     private:
        bool   d_debug;
        unsigned int d_pad_front;
        unsigned int d_pad_tail;

     protected:
      int calculate_output_stream_length(const gr_vector_int &ninput_items);

     public:
      zero_pad_impl(bool debug, unsigned int pad_front, unsigned int pad_tail);
      ~zero_pad_impl();

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

#endif /* INCLUDED_MIMO_OFDM_JRC_ZERO_PAD_IMPL_H */

