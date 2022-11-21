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

#ifndef INCLUDED_MIMO_OFDM_JRC_FRAME_DETECTOR_IMPL_H
#define INCLUDED_MIMO_OFDM_JRC_FRAME_DETECTOR_IMPL_H

#include <mimo_ofdm_jrc/frame_detector.h>
#include "utils.h"

namespace gr {
  namespace mimo_ofdm_jrc {

    class frame_detector_impl : public frame_detector
    {
     private:
      // Nothing to declare in this block.

     public:
      frame_detector_impl(int fft_len, int cp_len, double threshold, unsigned int min_n_peaks, unsigned int ignore_gap, bool debug);
      ~frame_detector_impl();
      
      	int d_fft_len;
        int d_cp_len;
        const double d_threshold = 0.0;
        const unsigned int d_min_n_peaks = 0;
        const int d_ignore_gap;// = 480;
        bool d_debug;

        const int MAX_PEAK_DISTANCE;
        const double MAX_PEAK_VALUE;
        const int MAX_SAMPLES;// = 540 * 80;

        enum {SEARCH, COPY} curr_state;
        int copied_samples = 0;
        int n_peaks = 0;
        float coarse_cfo_est = 0.0;

        uint64_t first_peak_ind = 0;
        int distance_to_firstpeak = 0;
        
        gr::thread::mutex 		d_mutex;
        
        void insert_tag(uint64_t item, double cfo_est, uint64_t input_item);
        
      // Where all the action really happens
      int general_work(int noutput_items,
           gr_vector_int &ninput_items,
           gr_vector_const_void_star &input_items,
           gr_vector_void_star &output_items);

    };

  } // namespace mimo_ofdm_jrc
} // namespace gr

#endif /* INCLUDED_MIMO_OFDM_JRC_FRAME_DETECTOR_IMPL_H */

