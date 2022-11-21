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

#ifndef INCLUDED_MIMO_OFDM_JRC_STREAM_ENCODER_IMPL_H
#define INCLUDED_MIMO_OFDM_JRC_STREAM_ENCODER_IMPL_H

#include <mimo_ofdm_jrc/stream_encoder.h>
#include <gnuradio/digital/constellation.h>
#include "utils.h"

namespace gr {
  namespace mimo_ofdm_jrc {

    class stream_encoder_impl : public stream_encoder
    {
     private:
     	bool         d_debug;
        char*        d_symbol_values;
        gr_complex*  d_complex_symbols;

        int         d_offset;
        int         d_symbol_len;
        int 		d_data_len;
        int 		d_N_ss_radar;
        ofdm_mcs	d_ofdm_mcs;
        MCS 		d_mod_encode;
        gr::thread::mutex d_mutex;

        uint8_t      d_scrambler;

        boost::shared_ptr<gr::digital::constellation> modulator;
        digital::constellation_bpsk::sptr d_bpsk;
        digital::constellation_qpsk::sptr d_qpsk;
        digital::constellation_16qam::sptr d_16qam;

     public:
      stream_encoder_impl(MCS mod_encode, int data_len, int N_ss_radar, bool debug);
      ~stream_encoder_impl();

      // Where all the action really happens
      // void forecast (int noutput_items, gr_vector_int &ninput_items_required);

      int general_work(int noutput_items,
           gr_vector_int &ninput_items,
           gr_vector_const_void_star &input_items,
           gr_vector_void_star &output_items);
      
      void set_mcs(MCS mod_encode);
      void print_output(const char *pdu, int pdu_length);

    };

  } // namespace mimo_ofdm_jrc
} // namespace gr

#endif /* INCLUDED_MIMO_OFDM_JRC_STREAM_ENCODER_IMPL_H */

