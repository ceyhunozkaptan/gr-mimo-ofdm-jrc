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

#ifndef INCLUDED_MIMO_OFDM_JRC_USRP_MIMO_TRX_IMPL_H
#define INCLUDED_MIMO_OFDM_JRC_USRP_MIMO_TRX_IMPL_H

#include <mimo_ofdm_jrc/usrp_mimo_trx.h>
#include <uhd/utils/thread.hpp>
#include <uhd/usrp/multi_usrp.hpp>
#include "utils.h"
#include <boost/thread/mutex.hpp>


namespace gr {
  namespace mimo_ofdm_jrc {

    class usrp_mimo_trx_impl : public usrp_mimo_trx
    {
     private:
      // Nothing to declare in this block.

     protected:
      int calculate_output_stream_length(const gr_vector_int &ninput_items);

     public:
      usrp_mimo_trx_impl(int N_mboard, int N_tx, int N_rx, int samp_rate, float center_freq, int num_delay_samps, bool debug, float update_period,
				std::string args, std::string clock_sources, std::string time_sources, 
        std::string antenna_tx, float gain_tx, float timeout_tx, float wait_tx, std::string wire_tx, 
				std::string antenna_rx, float gain_rx, float timeout_rx, float wait_rx, float lo_offset_rx, std::string wire_rx, 
				const std::string& len_key);
      ~usrp_mimo_trx_impl();
      
        void transmit();
        void receive();
        void set_num_delay_samps(int num_samps);
        void set_rx_gain(float gain);
        void set_tx_gain(float gain);
        
        bool d_debug;

        int d_N_mboard, d_N_tx, d_N_rx;

        int d_samp_rate;
        int d_fft_len;
        float d_center_freq;
        int d_num_delay_samps;
        double prev_tx_time;
        float d_update_period;

        uhd::usrp::multi_usrp::sptr d_usrp_mimo; 
        std::string d_args; 
        std::string d_wire_tx, d_wire_rx;
        std::vector<std::string> d_antenna_tx, d_antenna_rx;
        std::vector<std::string> d_clock_sources;
        std::vector<std::string> d_time_sources;
        // uhd::usrp::multi_usrp::sptr d_usrp_rx;

        uhd::tune_request_t d_tune_request_tx, d_tune_request_rx;
        uhd::tx_streamer::sptr d_tx_streamer;
        uhd::rx_streamer::sptr d_rx_streamer;
        uhd::tx_metadata_t d_metadata_tx;
        uhd::rx_metadata_t d_metadata_rx;
        
        double d_lo_offset_tx, d_lo_offset_rx;
        float d_timeout_tx, d_timeout_rx;
        float d_wait_tx, d_wait_rx;
        float d_gain_tx, d_gain_rx;
        
        uhd::time_spec_t d_time_now_tx, d_time_now_rx, d_current_time_now;
        

        gr::thread::thread d_thread_recv;
        gr::thread::thread d_thread_send;

        std::vector<std::vector<gr_complex>> d_out_buffer;

        std::vector<gr_complex *> d_out_recv_ptrs;
        std::vector<const void*>d_in_send_ptrs;

        int d_noutput_items_rx;
        int d_noutput_items_tx;    

        pmt::pmt_t d_time_key, d_time_val, d_srcid;

        boost::recursive_mutex d_mutex;

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

#endif /* INCLUDED_MIMO_OFDM_JRC_USRP_MIMO_TRX_IMPL_H */

