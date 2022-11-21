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

#ifndef INCLUDED_MIMO_OFDM_JRC_FRAME_SYNC_IMPL_H
#define INCLUDED_MIMO_OFDM_JRC_FRAME_SYNC_IMPL_H

#include <mimo_ofdm_jrc/frame_sync.h>
#include <gnuradio/filter/fir_filter.h>
#include "utils.h"
#include <list>
#include <tuple>
#include <gnuradio/fft/fft.h>

namespace gr {
  namespace mimo_ofdm_jrc {

    class frame_sync_impl : public frame_sync
    {
     private:
      enum {SYNC, COPY, RESET} d_state;
        int d_fft_len;
        int d_cp_len;
        
        int         total_out_count;
        int         d_first;

        int         sample_offset;
        int         d_frame_start;
        float       d_freq_offset;
        double      d_cfo_coarse_est;

        int frame_tag_offset = -1;

        gr_complex *d_correlation;
        std::list<std::pair<gr_complex, int> > d_cor;
        std::vector<gr::tag_t> d_tags;

        const bool d_debug;
        const int  SYNC_LENGTH;

        static const std::vector<gr_complex> LONG;
        gr::filter::kernel::fir_filter_ccc d_ltf_fir;

        gr::thread::mutex 		d_mutex;

     public:
      frame_sync_impl(int fft_len, int cp_len, unsigned int sync_length, std::vector<gr_complex> ltf_seq_time, bool debug);
      ~frame_sync_impl();

      // Where all the action really happens
      void forecast (int noutput_items, gr_vector_int &ninput_items_required);

      int general_work(int noutput_items,
           gr_vector_int &ninput_items,
           gr_vector_const_void_star &input_items,
           gr_vector_void_star &output_items);
      void search_frame_start();
    };

  } // namespace mimo_ofdm_jrc
} // namespace gr

#endif /* INCLUDED_MIMO_OFDM_JRC_FRAME_SYNC_IMPL_H */

