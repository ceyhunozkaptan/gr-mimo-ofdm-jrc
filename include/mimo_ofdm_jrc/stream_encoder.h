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

#ifndef INCLUDED_MIMO_OFDM_JRC_STREAM_ENCODER_H
#define INCLUDED_MIMO_OFDM_JRC_STREAM_ENCODER_H

#include <mimo_ofdm_jrc/api.h>
#include <gnuradio/block.h>

enum MCS : uint8_t {
    BPSK_1_2  = 0,
    BPSK_3_4  = 1,
    QPSK_1_2  = 2,
    QPSK_3_4  = 3,
    QAM16_1_2 = 4,
    QAM16_3_4 = 5,
};

enum PACKET_TYPE : uint8_t {
    NDP = 1,
    DATA = 2
};

namespace gr {
  namespace mimo_ofdm_jrc {

    /*!
     * \brief <+description of block+>
     * \ingroup mimo_ofdm_jrc
     *
     */
    class MIMO_OFDM_JRC_API stream_encoder : virtual public gr::block
    {
     public:
      typedef boost::shared_ptr<stream_encoder> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of mimo_ofdm_jrc::stream_encoder.
       *
       * To avoid accidental use of raw pointers, mimo_ofdm_jrc::stream_encoder's
       * constructor is in a private implementation
       * class. mimo_ofdm_jrc::stream_encoder::make is the public interface for
       * creating new instances.
       */
      static sptr make(MCS mod_encode, int data_len, int N_ss_radar, bool debug);
      virtual void set_mcs(MCS mod_encode) = 0;
    };

  } // namespace mimo_ofdm_jrc
} // namespace gr

#endif /* INCLUDED_MIMO_OFDM_JRC_STREAM_ENCODER_H */

