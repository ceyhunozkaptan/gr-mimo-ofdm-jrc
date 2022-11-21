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

#ifndef INCLUDED_MIMO_OFDM_JRC_MOVING_AVG_IMPL_H
#define INCLUDED_MIMO_OFDM_JRC_MOVING_AVG_IMPL_H

#include <mimo_ofdm_jrc/moving_avg.h>
#include "utils.h"

namespace gr {
  namespace mimo_ofdm_jrc {

    class moving_avg_impl : public moving_avg
    {
     private:
      	int d_length;
        float d_scale;
        int d_max_iter;
        std::vector<gr_complex> d_sum;

        int d_new_length;
        float d_new_scale;
        bool d_updated;
        bool d_debug;

     public:
      moving_avg_impl(int length, float scale, int max_iter, bool debug);
      ~moving_avg_impl();
      int length() const { return d_new_length; }
      float scale() const { return d_new_scale; }

      void set_length_and_scale(int length, float scale);
      void set_length(int length);
      void set_scale(float scale);

      // Where all the action really happens
      int work(
              int noutput_items,
              gr_vector_const_void_star &input_items,
              gr_vector_void_star &output_items
      );
    };

  } // namespace mimo_ofdm_jrc
} // namespace gr

#endif /* INCLUDED_MIMO_OFDM_JRC_MOVING_AVG_IMPL_H */

