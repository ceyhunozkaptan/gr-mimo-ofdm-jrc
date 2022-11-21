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

#ifndef INCLUDED_MIMO_OFDM_JRC_MIMO_OFDM_EQUALIZER_IMPL_H
#define INCLUDED_MIMO_OFDM_JRC_MIMO_OFDM_EQUALIZER_IMPL_H

#include <mimo_ofdm_jrc/mimo_ofdm_equalizer.h>
#include <gnuradio/digital/constellation.h>
#include "viterbi_decoder.h"
#include "utils.h"
#include <Eigen/Dense>

namespace gr {
  namespace mimo_ofdm_jrc {

    class mimo_ofdm_equalizer_impl : public mimo_ofdm_equalizer
    {
     private:
      	// gr::thread::mutex 		d_mutex;
        std::vector<gr::tag_t> 	tags;

        std::vector<int> 			d_pilot_carriers; 
        std::vector<int>			d_data_carriers;
        std::vector<int>			d_nondata_carriers;

        std::vector<int>			d_active_carriers;

        const std::string d_chan_est_file; 
        const std::string d_comm_log_file;

        const std::vector<std::vector<gr_complex>> d_pilot_symbols;
        const std::vector<std::vector<gr_complex>> d_mapped_ltf_symbols;
        int         d_N_mimo_ltf;
        int         d_N_tx_antenna;
        Eigen::Matrix<gr_complex,Eigen::Dynamic,Eigen::Dynamic>  mapped_ltf_mat;

        bool 		d_debug, d_stats_record;	

        bool        d_stat_record_update;

        int  		d_symbol_ind;

        int 		d_fft_len;
        int 		d_cp_len;
        int 		d_N_pilot_carr;
        int 		d_N_data_carr;
        int 		d_N_active_carr;


        int 		total_out;
        double 		d_freq;  // Hz
        double 		d_bw;  // Hz

        double 		d_freq_offset_from_synclong; 
        double 		d_er;
        double 		d_epsilon0;
        double 		d_snr_est;
        double 		d_precoded_snr_est;

        pmt::pmt_t chan_est_mean_pmt;

        bool sig_decode_success;
        bool equalize_done;

        int         data_length_SIG; 
        int         n_ofdm_symbols_SIG;
        MCS         mcs_SIG;
        PACKET_TYPE         packet_type_SIG;
        ChannelEstimator d_estimator_algo;

        double signal_power_sum;
        double noise_power_sum;
        int snr_est_count;

        gr_complex* 	d_H;
        gr_complex* 	d_H_mimo;
        gr_complex* d_H_mimo_temp;
      
        uint8_t* rx_bits;

        std::vector<gr_complex> d_lft_seq;
        std::vector<gr_complex> d_est_rx_pilots;

        gr_complex* 	d_prev_pilots;

        viterbi_decoder d_decoder;

        boost::shared_ptr<gr::digital::constellation> modulator_SIG;

        digital::constellation_bpsk::sptr d_bpsk;
        digital::constellation_qpsk::sptr d_qpsk;
        digital::constellation_16qam::sptr d_16qam;
        // digital::constellation_64qam::sptr d_64qam;

        float perf_display_interval;       
        std::chrono::system_clock::time_point last_perf_time;

     public:
      mimo_ofdm_equalizer_impl(ChannelEstimator estimator_algo, 
                                    double freq, 
                                    double bw, 
                                    int fft_len, 
                                    int cp_len, 
                                    std::vector<int> data_carriers, 
                                    std::vector<int> pilot_carriers, 
                                    const std::vector<std::vector<gr_complex>>& pilot_symbols,
                                    std::vector<gr_complex> ltf_seq, 
                                    const std::vector<std::vector<gr_complex>>& mapped_ltf_symbols,
                                    int n_mimo_ltf,
                                    const std::string& chan_est_file,
                                    const std::string& comm_log_file,
                                    bool stats_record,
                                    bool debug);
      ~mimo_ofdm_equalizer_impl();

	bool decode_signal_field(gr_complex* rx_symbols);
        bool parse_signal(uint8_t* decoded_bits);  

        void symbol_equalize (gr_complex *in, int n, gr_complex *symbols, ChannelEstimator algo); 
        double estimate_residual_cfo(gr_complex* rx_pilots, gr_complex* chan_est, const std::vector<gr_complex>& ref_pilots, std::vector<gr_complex>& est_rx_pilots);
        void set_estimator(ChannelEstimator algo);
        void set_bandwidth(double bw);
        void set_frequency(double freq);
        void set_stats_record(bool stats_record);
        // Where all the action really happens
        void forecast (int noutput_items, gr_vector_int &ninput_items_required);

        int general_work(int noutput_items,
            gr_vector_int &ninput_items,
            gr_vector_const_void_star &input_items,
            gr_vector_void_star &output_items);

    };

  } // namespace mimo_ofdm_jrc
} // namespace gr

#endif /* INCLUDED_MIMO_OFDM_JRC_MIMO_OFDM_EQUALIZER_IMPL_H */

