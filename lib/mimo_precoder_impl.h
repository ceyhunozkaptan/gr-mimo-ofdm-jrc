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

#ifndef INCLUDED_MIMO_OFDM_JRC_MIMO_PRECODER_IMPL_H
#define INCLUDED_MIMO_OFDM_JRC_MIMO_PRECODER_IMPL_H

#include <mimo_ofdm_jrc/mimo_precoder.h>
#include <gnuradio/digital/constellation.h>
#include "utils.h"
#include <Eigen/Dense>
#include <Eigen/SVD>
// #include <Eigen/Core>

namespace gr {
  namespace mimo_ofdm_jrc {

    class mimo_precoder_impl : public mimo_precoder
    {
     private:
     	const int d_fft_len;
        const int d_N_tx;
        int d_N_ltf;
        const int d_N_ss;
        const std::vector<std::vector<gr_complex>> d_mapped_ltf_symbols;
        const bool d_debug;
        bool d_chan_est_smoothing, d_radar_aided;
        bool d_use_radar_streams;
        bool d_phased_steering;
        const std::string d_chan_est_file;
        const std::string d_radar_log_file;
        std::time_t last_chanEst_update_time;

        float angle_estimate;
        bool mean_chan_est_changed;
        int d_N_ss_radar;

        //! Which carriers/symbols carry data
        std::vector<int> d_data_carriers;
        //! Which carriers/symbols carry pilots symbols
        std::vector<int> d_pilot_carriers;
        //! Value of said pilot symbols
        const std::vector<std::vector<gr_complex>> d_pilot_symbols;
        //! Synch words
        const std::vector<std::vector<gr_complex>> d_sync_words;

        std::vector<int> d_active_carriers;
        int d_N_active_carriers;
        int d_N_data_carriers;

        gr_complex* signal_field_symbols;

        Eigen::Matrix<gr_complex,Eigen::Dynamic,Eigen::Dynamic>  mapped_ltf_mat;
        Eigen::MatrixXcf dft_matrix;
        // steering matrix for data
        std::vector<std::vector<gr_complex>> steering_matrix;
        // steering matrix averaged over subcarriers
        std::vector<gr_complex> steering_matrix_mean;
        // estimated mimo channel vector averaged over subcarriers
        std::vector<gr_complex> chan_est_vector_mean;

        digital::constellation_bpsk::sptr d_bpsk;
        digital::constellation_qpsk::sptr d_qpsk;

        float perf_display_interval;       
        std::chrono::system_clock::time_point last_perf_time;

     protected:
      int calculate_output_stream_length(const gr_vector_int &ninput_items);

     public:
      mimo_precoder_impl(int fft_len,
                            int N_tx,
                            int N_ss,
                            const std::vector<int>& data_carriers,
                            const std::vector<int>& pilot_carriers,
                            const std::vector<std::vector<gr_complex>>& pilot_symbols,
                            const std::vector<std::vector<gr_complex>>& sync_words,
                            const std::vector<std::vector<gr_complex>>& mapped_ltf_symbols,
                            const std::string& chan_est_file,
                            bool chan_est_smoothing,
                            const std::string& radar_log_file,
                            bool radar_aided,
                            bool phased_steering,
                            bool use_radar_streams,
                            const std::string& len_tag_key,
                            bool debug);
      ~mimo_precoder_impl();
      
      int get_bit(int b, int i);
        void generate_signal_field(gr_complex* out, packet_param& frame, ofdm_mcs& ofdm);
        std::vector<std::vector<gr_complex>> get_dft_matrix(int N);
        Eigen::MatrixXcf get_dft_matrix_eigen(int N);
        bool compute_steering_matrix();
        bool compute_radar_aided_steering();
        void set_chan_est_smoothing(bool chan_est_smoothing);
        void set_radar_aided(bool radar_aided);
        void set_use_radar_streams(bool use_radar_streams);
        void set_phased_steering(bool phased_steering);

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

#endif /* INCLUDED_MIMO_OFDM_JRC_MIMO_PRECODER_IMPL_H */

