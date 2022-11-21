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

#ifndef INCLUDED_MIMO_OFDM_JRC_USRP_MIMO_TRX_H
#define INCLUDED_MIMO_OFDM_JRC_USRP_MIMO_TRX_H

#include <mimo_ofdm_jrc/api.h>
#include <gnuradio/tagged_stream_block.h>

namespace gr {
  namespace mimo_ofdm_jrc {

    /*!
     * \brief <+description of block+>
     * \ingroup mimo_ofdm_jrc
     *
     */
    class MIMO_OFDM_JRC_API usrp_mimo_trx : virtual public gr::tagged_stream_block
    {
     public:
      typedef boost::shared_ptr<usrp_mimo_trx> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of mimo_ofdm_jrc::usrp_mimo_trx.
       *
       * To avoid accidental use of raw pointers, mimo_ofdm_jrc::usrp_mimo_trx's
       * constructor is in a private implementation
       * class. mimo_ofdm_jrc::usrp_mimo_trx::make is the public interface for
       * creating new instances.
       */
      static sptr make(int N_mboard, int N_tx, int N_rx, int samp_rate, float center_freq, int num_delay_samps, bool debug,  float update_period,
        std::string args, std::string clock_sources, std::string time_sources, 
        std::string antenna_tx, float gain_tx, float timeout_tx, float wait_tx, std::string wire_tx, 
        std::string antenna_rx, float gain_rx, float timeout_rx, float wait_rx, float lo_offset_rx, std::string wire_rx, 
        const std::string& len_key="packet_len");
        
        virtual void set_num_delay_samps(int num_samps) = 0;
        virtual void set_rx_gain(float gain) = 0;
        virtual void set_tx_gain(float gain) = 0;
    };

  } // namespace mimo_ofdm_jrc
} // namespace gr

#endif /* INCLUDED_MIMO_OFDM_JRC_USRP_MIMO_TRX_H */

