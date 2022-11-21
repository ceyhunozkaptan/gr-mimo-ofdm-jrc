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

#ifndef INCLUDED_MIMO_OFDM_JRC_STREAM_DECODER_IMPL_H
#define INCLUDED_MIMO_OFDM_JRC_STREAM_DECODER_IMPL_H

#include <mimo_ofdm_jrc/stream_decoder.h>
#include "viterbi_decoder.h"
#include <boost/crc.hpp>
#include "utils.h"
#include <gnuradio/digital/constellation.h>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/statistics/rolling_mean.hpp>

#define time_now std::chrono::high_resolution_clock::now()

namespace ba = boost::accumulators;
namespace bt = ba::tag;

namespace gr {
  namespace mimo_ofdm_jrc {

    class stream_decoder_impl : public stream_decoder
    {
     private:
      	bool d_debug, d_stats_record;
        bool d_new_stat_started;
        bool d_log;

        const std::string d_comm_log_file;

        packet_param d_stream_param;
        ofdm_mcs d_ofdm_mcs;

        MCS d_mcs;
        PACKET_TYPE d_packet_type;
        int d_max_ofdm_len;
        int d_n_data_carriers;

        int total_data_received;

        std::chrono::high_resolution_clock::time_point curr_time, last_update_time;

        double elapsed_time;

        ba::accumulator_set<int, ba::stats<bt::rolling_mean> > per_stats;
        ba::accumulator_set<float, ba::stats<bt::rolling_mean> > snr_data_stats;

        // double d_nom_freq;  // nominal frequency, Hz
        // double d_freq_offset;  // frequency offset, Hz
        viterbi_decoder d_viterbi_decoder;
        boost::shared_ptr<gr::digital::constellation> d_demodulator;

        digital::constellation_bpsk::sptr 				d_bpsk;
        digital::constellation_qpsk::sptr 				d_qpsk;
        digital::constellation_16qam::sptr 				d_16qam;

        uint8_t* d_rx_symbols;
        
        uint8_t* d_rx_bits;//[MAX_ENCODED_BITS];
        uint8_t d_deinterleaved_bits[MAX_ENCODED_BITS];
        uint8_t out_bytes[MAX_PAYLOAD_SIZE + 2]; // 2 for signal field
        
        float d_snr_est;  // dB        
        float d_snr_data_est;  // dB
        uint d_data_length;
        std::vector<gr_complex> chan_est_mean;
        
        const static int info_bytes = 2 + 2*sizeof(float); // success(1byte) + packet_type(1bytes) + snr(4bytes) + snr_precoded(4bytes)
        uint8_t send_out_bytes[MAX_PAYLOAD_SIZE + 2 + info_bytes]; // 3 additional bytes for info --> success(1byte) + snr(4bytes) + packet_type(1bytes)


        gr::thread::mutex d_mutex;

        float perf_display_interval;       
        std::chrono::system_clock::time_point last_perf_time;
        
        int n_copied;
        bool d_frame_rx_complete;
        bool d_start_decoding;

     public:
      stream_decoder_impl(int n_data_carriers,
                            const std::string& comm_log_file,
                            bool stats_record, 
                            bool debug);
      ~stream_decoder_impl();

      // void forecast (int noutput_items, gr_vector_int &ninput_items_required);

        int general_work(int noutput_items,
            gr_vector_int &ninput_items,
            gr_vector_const_void_star &input_items,
            gr_vector_void_star &output_items);

        void decode();
        // void deinterleave();
        void descramble (uint8_t *decoded_bits);
        void print_output(const uint8_t *pdu, int pdu_length);
        void setup_demodulator(MCS mcs_value);

        void set_stats_record(bool stats_record);

    };

  } // namespace mimo_ofdm_jrc
} // namespace gr

#endif /* INCLUDED_MIMO_OFDM_JRC_STREAM_DECODER_IMPL_H */

