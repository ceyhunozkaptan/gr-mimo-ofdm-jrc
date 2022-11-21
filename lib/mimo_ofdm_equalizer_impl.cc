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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include "mimo_ofdm_equalizer_impl.h"
#include "utils.h"

#include<fstream>

namespace gr {
  namespace mimo_ofdm_jrc {

    mimo_ofdm_equalizer::sptr
    mimo_ofdm_equalizer::make(ChannelEstimator estimator_algo, 
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
                                bool debug)
    {
      return gnuradio::get_initial_sptr
        (new mimo_ofdm_equalizer_impl(estimator_algo, 
                                        freq, 
                                        bw, 
                                        fft_len, 
                                        cp_len, 
                                        data_carriers, 
                                        pilot_carriers, 
                                        pilot_symbols, 
                                        ltf_seq, 
                                        mapped_ltf_symbols,
                                        n_mimo_ltf,
                                        chan_est_file,
                                        comm_log_file,
                                        stats_record,
                                        debug));

    }


    /*
     * The private constructor
     */
    mimo_ofdm_equalizer_impl::mimo_ofdm_equalizer_impl(ChannelEstimator estimator_algo, 
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
                                                        bool debug)
      : gr::block("mimo_ofdm_equalizer",
                gr::io_signature::make(1, 1, fft_len * sizeof(gr_complex)),
                gr::io_signature::make(1, 1, data_carriers.size() * sizeof(gr_complex))),
                d_freq(freq), 
                d_bw(bw), 
                d_fft_len(fft_len), 
                d_cp_len(cp_len), 
                d_data_carriers(data_carriers), 
                d_pilot_carriers(pilot_carriers), 
                d_pilot_symbols(pilot_symbols),
                d_prev_pilots(NULL),
                d_estimator_algo(estimator_algo), 
                d_symbol_ind(0), 
                d_freq_offset_from_synclong(0.0),
                d_lft_seq(ltf_seq),
                d_mapped_ltf_symbols(mapped_ltf_symbols),
                d_N_mimo_ltf(n_mimo_ltf),
                d_chan_est_file(chan_est_file),
                d_comm_log_file(comm_log_file),
                d_stats_record(stats_record),
                d_debug(debug)
    {
            perf_display_interval = 1.0;
            last_perf_time = std::chrono::system_clock::now();

            set_tag_propagation_policy(block::TPP_DONT);

            d_N_pilot_carr = d_pilot_carriers.size();
            d_N_data_carr = d_data_carriers.size();
            dout << "[OFDM Equalizer] d_N_data_carr " << d_N_data_carr << std::endl;

          	d_prev_pilots = (gr_complex*) calloc(d_N_pilot_carr, sizeof(gr_complex));
          	d_H = (gr_complex*) calloc(d_fft_len, sizeof(gr_complex));
            d_H_mimo = (gr_complex*) calloc(d_fft_len, sizeof(gr_complex));
            d_H_mimo_temp = (gr_complex*) calloc(d_fft_len, sizeof(gr_complex));

			d_bpsk = digital::constellation_bpsk::make(); // for decoding SIG symbol 
            d_qpsk = digital::constellation_qpsk::make();
            d_16qam = digital::constellation_16qam::make();
            // d_64qam = digital::constellation_64qam::make();

			set_estimator(estimator_algo);
            
            d_est_rx_pilots.resize(d_N_pilot_carr);

	      	for (int i = 0; i < d_N_pilot_carr; i++) {
                d_pilot_carriers[i] += d_fft_len/2;
                dout << "[OFDM Equalizer] d_pilot_carriers " << d_pilot_carriers[i] << std::endl;
			}

	      	for (int i = 0; i < d_N_data_carr; i++) {
                d_data_carriers[i] += d_fft_len/2;
                dout << "[OFDM Equalizer] d_data_carriers " << d_data_carriers[i] << std::endl;
			}

	      	for (int i_sc = 0; i_sc < d_fft_len; i_sc++) 
              {
                // d_all_carriers[i_sc] = i_sc;
                if (std::find(d_data_carriers.begin(), d_data_carriers.end(), i_sc) == d_data_carriers.end())
                {
                    d_nondata_carriers.push_back(i_sc);
                    dout << "[OFDM Equalizer] d_nondata_carriers " << i_sc << std::endl;
                }
			}

            d_active_carriers.reserve(d_data_carriers.size() + d_pilot_carriers.size());
            d_active_carriers.insert(d_active_carriers.end(), d_data_carriers.begin(), d_data_carriers.end());
            d_active_carriers.insert(d_active_carriers.end(), d_pilot_carriers.begin(), d_pilot_carriers.end());
            d_N_active_carr = d_active_carriers.size();

            sort(d_active_carriers.begin(), d_active_carriers.end());


            mapped_ltf_mat.resize(d_fft_len, d_mapped_ltf_symbols[0].size());

            for (int i_row = 0; i_row < d_mapped_ltf_symbols.size(); i_row++){
                mapped_ltf_mat.row(i_row) = Eigen::VectorXcf::Map(&d_mapped_ltf_symbols[i_row][0],d_mapped_ltf_symbols[i_row].size());
            }

            d_N_tx_antenna = d_mapped_ltf_symbols[0].size() / d_N_mimo_ltf;
            dout << "[OFDM Equalizer] # of MIMO LTF Symbols: " << d_N_mimo_ltf << ", # of Tx Antennas: " <<  d_N_tx_antenna << std::endl;

            d_stat_record_update = false;

            // rx_bits = new uint8_t[d_N_data_carr];
            rx_bits = (uint8_t*) calloc(d_N_data_carr, sizeof(uint8_t));
    }

    /*
     * Our virtual destructor.
     */
    mimo_ofdm_equalizer_impl::~mimo_ofdm_equalizer_impl()
    {
        free(d_prev_pilots);
        free(d_H);
        free(d_H_mimo);
        free(d_H_mimo_temp);
        // delete [] rx_bits;
        free(rx_bits);
    }

    void
    mimo_ofdm_equalizer_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
        ninput_items_required[0] = noutput_items;
    }

    int
    mimo_ofdm_equalizer_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
        gr::block::set_thread_priority(50);
        gr::thread::scoped_lock lock(d_setlock);

        const gr_complex *in = (const gr_complex *) input_items[0];
        gr_complex *out = (gr_complex *) output_items[0];

        int n_in = 0;
        int n_out = 0;

        gr_complex eq_data_symbols_Z[d_N_data_carr];
        gr_complex rx_symbol_Y[d_fft_len];
        gr_complex rx_mimo_preamble[d_fft_len][d_N_mimo_ltf];


        // dout << "=================================================================================" << std::endl;
        // dout << "[OFDM Equalizer] input " << ninput_items[0] << "  output " << noutput_items << std::endl;

        while((n_in < ninput_items[0]) && (n_out < noutput_items)) 
        {
            get_tags_in_window(tags, 0, n_in, n_in + 1, pmt::string_to_symbol("frame_start"));

            // new frame
            if(tags.size()) 
            {
                d_symbol_ind = 0;
                total_out = 0;
                n_ofdm_symbols_SIG = 0;

                d_freq_offset_from_synclong = pmt::to_double(tags.front().value) * d_bw / (2 * M_PI);
                d_epsilon0 = pmt::to_double(tags.front().value) * d_bw / (2 * M_PI * d_freq);
                d_er = 0;

                // assume SIG field is correct for now
                sig_decode_success = true;
                equalize_done = false;

                // reset data stats 
                signal_power_sum = 0;
                noise_power_sum = 0;
                snr_est_count = 0;
                dout << "================================================================================" << std::endl;

                dout << "[OFDM Equalizer] New frame --> d_freq_offset_from_synclong = " << pmt::to_double(tags.front().value) << std::endl;
            }

            // dout << "[OFDM Equalizer] d_symbol_ind: " << d_symbol_ind << ", n_in: " << n_in << ", n_out: " << n_out << ", n_ofdm_symbols_SIG: " << n_ofdm_symbols_SIG << std::endl;

            // not interesting -> skip
            if( d_symbol_ind > (n_ofdm_symbols_SIG + 2 + d_N_mimo_ltf) || !sig_decode_success) 
            {
                n_in++;
                // dout << "[OFDM Equalizer] skipping, equalize_done: " << equalize_done << " d_symbol_ind: " << d_symbol_ind << std::endl;
                continue;
            }

            std::memcpy(rx_symbol_Y, in + n_in*d_fft_len, d_fft_len*sizeof(gr_complex));

            // compensate sampling offset 
            // E. Sourour, H. El-Ghoroury and D. McNeill, "Frequency offset estimation and correction in the IEEE 802.11a WLAN," IEEE 60th Vehicular Technology Conference, 2004. VTC2004-Fall. 2004, Los Angeles, CA, 2004, pp. 4923-4927 Vol. 7
            for(int i = 0; i < d_fft_len; i++) 
            {
                rx_symbol_Y[i] *= exp(gr_complex(0, 2*M_PI*d_symbol_ind*((d_fft_len+d_cp_len)*1.0/d_fft_len)*(d_epsilon0 + d_er)*(i-d_fft_len/2)));
            }

            gr_complex pilot_sum = gr_complex(0.0,0.0);
            gr_complex prev_pilot_sum = gr_complex(0.0,0.0);

            double residual_cfo_est;
            double er;
            // legacy-LTF-1 ====================================================================================
            if(d_symbol_ind == 0) 
            {					
                std::memcpy(d_H, rx_symbol_Y, d_fft_len * sizeof(gr_complex));
            } 
            // legacy-LTF-2 =================================================================================
            else if(d_symbol_ind == 1) 
            {
                double signal = 0;
                double noise = 0;
                for(int i_carr : d_active_carriers)
                {
                    noise += std::pow(std::abs(d_H[i_carr] - rx_symbol_Y[i_carr]), 2);
                    signal += std::pow(std::abs(d_H[i_carr] + rx_symbol_Y[i_carr]), 2);
                    d_H[i_carr] += rx_symbol_Y[i_carr];
                    d_H[i_carr] /= d_lft_seq[i_carr] * gr_complex(2, 0);
                }

                pilot_sum = gr_complex(0.0,0.0);
                for(int i_carr : d_pilot_carriers)
                {
                    // Implementation of a MIMO OFDM-based wireless LAN system --> WORKING IN SIM
                    // arg( Xest_n x Conj(P_n) ) --> exp(-1j beta)
                    pilot_sum += (rx_symbol_Y[i_carr]/d_H[i_carr]) * conj(d_lft_seq[i_carr]);

                }
                // dout << "[OFDM Equalizer] pilot_sum 2 = " << pilot_sum << std::endl;

                residual_cfo_est = arg(pilot_sum);
                for(int i = 0; i < d_fft_len; i++) 
                {
                    rx_symbol_Y[i] *= exp(gr_complex(0, -residual_cfo_est));
                }

                d_snr_est = 10 * std::log10( signal / noise / 2);
            } 
            //SIG Symbol ======================================================================================
            else if (d_symbol_ind == 2)
            {
                // Compensate Residual Frequency Offset
                residual_cfo_est = estimate_residual_cfo(rx_symbol_Y, d_H, d_pilot_symbols[0], d_est_rx_pilots);
                for(int i = 0; i < d_fft_len; i++) 
                {
                    rx_symbol_Y[i] *= exp(gr_complex(0, -residual_cfo_est));
                }

                // Equalization
                symbol_equalize(rx_symbol_Y, d_symbol_ind, eq_data_symbols_Z, d_estimator_algo);
                
                dout << "[OFDM Equalizer] Decoding SIG field..." << std::endl;
                // bool decode_success = decode_signal_field(eq_data_symbols_Z);
                sig_decode_success = decode_signal_field(eq_data_symbols_Z);

                if (sig_decode_success)
                {                    
                    dout << "[OFDM Equalizer] PACKET TYPE from SIG: " << packet_type_SIG << std::endl;
                    dout << "[OFDM Equalizer] MCS Index from SIG: " << mcs_SIG << std::endl;
                    dout << "[OFDM Equalizer] DATA LENGTH (byte) from SIG: " << data_length_SIG << std::endl;
                    dout << "[OFDM Equalizer] # of OFDM Symbol for DATA: " << n_ofdm_symbols_SIG << std::endl;

                    pmt::pmt_t dict = pmt::make_dict();
                    dict = pmt::dict_add(dict, pmt::mp("data_bytes"), pmt::from_uint64(data_length_SIG));
                    dict = pmt::dict_add(dict, pmt::mp("mcs"), pmt::from_uint64(mcs_SIG));
                    dict = pmt::dict_add(dict, pmt::mp("packet_type"), pmt::from_uint64(packet_type_SIG));
                    dict = pmt::dict_add(dict, pmt::mp("snr"), pmt::from_double(d_snr_est));
                    dict = pmt::dict_add(dict, pmt::mp("freq_offset"), pmt::from_double(d_freq_offset_from_synclong));
                    add_item_tag(0, nitems_written(0) + n_out, pmt::string_to_symbol("stream_start"),	dict, pmt::string_to_symbol(alias()));
                }
                else
                {
                    //FIXME do we need to do anything at this point?
                    // n_in = ninput_items[0]-1;
                }
            }
            // MIMO LTF Symbols ======================================================================================
            else if( 3 <= d_symbol_ind && d_symbol_ind <= 2 + d_N_mimo_ltf )  
            {
                //MIMO CHANNEL ESTIMATION
                int ltf_ind = d_symbol_ind - 3;

                // TODO Requires non-orthogonal pilots on MIMO-LTF
                // std::vector<gr_complex> ref_pilots (d_pilot_symbols[0].size());
                // pilot_sum = gr_complex(0.0,0.0);
                // for(int i_carr : d_pilot_carriers)
                // {
                //     // pilot_sum += (rx_symbol_Y[i_carr]/d_H[i_carr]) * conj(d_mapped_ltf_symbols[i_sc][0]);
                // }
                // residual_cfo_est = arg(pilot_sum);
                // for(int i = 0; i < d_fft_len; i++) 
                // {
                //     rx_symbol_Y[i] *= exp(gr_complex(0, -residual_cfo_est));
                // }
                

                // Storing MIMO preambles
                for (int i_sc = 0; i_sc < d_fft_len; i_sc++) 
                {
                    // rx_mimo_preamble[ltf_ind + i_sc*d_N_mimo_ltf] = rx_symbol_Y[i_sc];
                    rx_mimo_preamble[i_sc][ltf_ind] = rx_symbol_Y[i_sc];
                }
                    
                // Its time to estimate MIMO channel
                if (ltf_ind == d_N_mimo_ltf-1) 
                {
                    if (packet_type_SIG == PACKET_TYPE::NDP)
                    { 
                        // No precoding --> Estimate channel
                        const static Eigen::IOFormat csv_formatting(Eigen::FullPrecision, Eigen::DontAlignCols, 
                                                                ";",   //_coeffSeparator
                                                                ":",   //_rowSeparator
                                                                "",     //_rowPrefix
                                                                "",     //_rowSuffix
                                                                "",  //_matPrefix
                                                                "\n");  //_matSuffix
                        std::ofstream file_stream(d_chan_est_file, std::ofstream::trunc);
                        // TODO fix eigen3 mapping here
                        // Eigen::Matrix<gr_complex,Eigen::Dynamic,Eigen::Dynamic, Eigen::RowMajor> X_ltf;
                        
                        // std::vector<gr_complex> chan_est_vector_mean(d_N_tx_antenna);
                        Eigen::VectorXcf chan_est_vector_mean = Eigen::VectorXcf::Zero(d_N_tx_antenna);

                        for (int i_sc = 0; i_sc < d_fft_len; i_sc++)
                        {                      
                            // Channel matrix for 1 RX antenna
                            Eigen::Map<Eigen::VectorXcf> y_rx(rx_mimo_preamble[i_sc], d_N_mimo_ltf);
                            Eigen::Map< const Eigen::Matrix<gr_complex, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > X_ltf(&d_mapped_ltf_symbols[i_sc][0], d_N_tx_antenna, d_N_mimo_ltf);
                            Eigen::VectorXcf H_mimo(d_N_tx_antenna);
                            H_mimo.noalias() = X_ltf.conjugate()*y_rx;

                            if (file_stream.is_open())
                            {
                                file_stream << i_sc << ":";
                                file_stream << H_mimo.transpose().format(csv_formatting);
                                file_stream.flush();
                            }
                            else
                            {
                                throw std::runtime_error("[OFDM Equalizer] Could not open file!!");
                            }
                            
                            if (std::find(d_active_carriers.begin(), d_active_carriers.end(), i_sc) != d_active_carriers.end())
                            {
                                chan_est_vector_mean += H_mimo;
                            }
                        }
                        file_stream.close();

                        chan_est_vector_mean = chan_est_vector_mean/d_N_active_carr;

                        // std::cout << chan_est_vector_mean << std::endl;
                        chan_est_mean_pmt = pmt::init_c32vector(d_N_tx_antenna, chan_est_vector_mean.data());
                    }
                    else if (packet_type_SIG == PACKET_TYPE::DATA)
                    {
                        // TODO this may go into SIG field or block parameters
                        int rx_index = 0; //This is the actual stream index at the transmitter.

                        Eigen::Map<Eigen::VectorXcf> y_rx(rx_mimo_preamble[0], d_N_mimo_ltf);
                        Eigen::Map< const Eigen::Matrix<gr_complex, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > X_ltf( &d_mapped_ltf_symbols[0][0], d_N_tx_antenna, d_N_mimo_ltf);
                        gr_complex chan_est_vector_mean = gr_complex(0.0,0.0);

                        // MIMO Channel Estimation with MIMO Preamble
                        for(int i_data = 0; i_data < d_data_carriers.size(); i_data++) 
                        {   
                            int i_sc = d_data_carriers[i_data];

                            new (&y_rx) Eigen::Map<Eigen::VectorXcf> (rx_mimo_preamble[i_sc], d_N_mimo_ltf);
                            new (&X_ltf) Eigen::Map< const Eigen::Matrix<gr_complex, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > ( &d_mapped_ltf_symbols[i_sc][0], d_N_tx_antenna, d_N_mimo_ltf);
                            d_H_mimo[i_sc] = X_ltf.row(rx_index).dot(y_rx) / (float) d_N_mimo_ltf; 
                            chan_est_vector_mean += d_H_mimo[i_sc];
                            // dout << "(" << d_H_mimo[i_sc] << "), ";   
                        }

                        for(int i_pilot = 0; i_pilot < d_pilot_carriers.size(); i_pilot++) 
                        {   
                            int i_sc = d_pilot_carriers[i_pilot];

                            new (&y_rx) Eigen::Map<Eigen::VectorXcf> (rx_mimo_preamble[i_sc], d_N_mimo_ltf);
                            new (&X_ltf) Eigen::Map< const Eigen::Matrix<gr_complex, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > ( &d_mapped_ltf_symbols[i_sc][0], d_N_tx_antenna, d_N_mimo_ltf);                            // Eigen::VectorXcf H_mimo(1);
                            d_H_mimo[i_sc] = X_ltf.row(rx_index).dot(y_rx) / (float) d_N_mimo_ltf; 
                            chan_est_vector_mean += d_H_mimo[i_sc];
                        }

                        chan_est_vector_mean = chan_est_vector_mean / (float) d_N_active_carr;
                        chan_est_mean_pmt = pmt::init_c32vector(1, &chan_est_vector_mean);
                    }
                    else
                    {
                        // throw std::runtime_error("[OFDM Equalizer] Something is wrong --> Packet type is not correct!");
                        std::cout << "[OFDM Equalizer] Something is wrong --> Packet type is not correct: " << packet_type_SIG << std::endl;
                    }
                }
            }
            //DATA SYMBOLS
            else
            {
                std::vector<gr_complex> ref_pilots = d_pilot_symbols[(d_symbol_ind - 3 - d_N_mimo_ltf) % d_pilot_symbols.size()];

                if (packet_type_SIG == PACKET_TYPE::NDP)
                {
                    residual_cfo_est = estimate_residual_cfo(rx_symbol_Y, d_H, ref_pilots, d_est_rx_pilots);
                }
                else
                {
                    residual_cfo_est = estimate_residual_cfo(rx_symbol_Y, d_H_mimo, ref_pilots, d_est_rx_pilots);
                }

                // Compensate Residual Frequency Offset
                for(int i = 0; i < d_fft_len; i++) 
                {
                    rx_symbol_Y[i] *= exp(gr_complex(0, -residual_cfo_est));
                }

                for (int i_pilot = 0; i_pilot < d_N_pilot_carr; i_pilot++)
                {
                    // "A Novel SNR Estimation Algorithm for OFDM"
                    signal_power_sum += real( d_est_rx_pilots[i_pilot]*conj(d_est_rx_pilots[i_pilot]) );

                    gr_complex pilotErr = d_est_rx_pilots[i_pilot] - rx_symbol_Y[d_pilot_carriers[i_pilot]];
                    noise_power_sum += real(pilotErr*conj(pilotErr)); // MMSE estimate of noise
                    
                    snr_est_count++;
                }

              
                if (packet_type_SIG == PACKET_TYPE::NDP)
                {
                    symbol_equalize(rx_symbol_Y, d_symbol_ind, eq_data_symbols_Z, d_estimator_algo);

                    if (d_estimator_algo == STA)
                    {
                        uint8_t symbol_value;
                        gr_complex demod_symbol_X, H_update;
                        float symbol_mse = 0;
                        float alpha = 0.5;
                        for(int i_data = 0; i_data < d_data_carriers.size(); i_data++) 
                        {
                            int i_sc = d_data_carriers[i_data];
                            symbol_value = modulator_SIG->decision_maker(&eq_data_symbols_Z[i_data]);
                            modulator_SIG->map_to_points(symbol_value, &demod_symbol_X);
                            if(modulator_SIG->bits_per_symbol() == d_qpsk->bits_per_symbol())
                            {
                                demod_symbol_X = demod_symbol_X / (float)2.0;
                            }
                            symbol_mse += std::pow(abs(eq_data_symbols_Z[i_data] - demod_symbol_X),2);

                            // std::cout << "eq_data_symbols_Z: " << eq_data_symbols_Z[i_data] << ", demod_symbol_X: " << demod_symbol_X << std::endl;
                            H_update = rx_symbol_Y[i_sc]/demod_symbol_X;
                            // std::cout << "d_H[i_sc]: " << d_H[i_sc] << ", H_update: " << H_update << std::endl;

                            d_H[i_sc] = gr_complex(1 - alpha, 0) * d_H[i_sc] + gr_complex(alpha, 0) * H_update;
                        }
                        symbol_mse = symbol_mse / d_data_carriers.size();
                        
                        for(int i_pilot = 0; i_pilot < d_pilot_carriers.size(); i_pilot++) 
                        {
                            int i_sc = d_pilot_carriers[i_pilot];

                            // symbol_mse += std::pow(abs(eq_data_symbols_Z[i_data] - demod_symbol_X),2);
                            d_H[i_sc] = gr_complex(1 - alpha, 0) * d_H[i_sc] + gr_complex(alpha, 0) * rx_symbol_Y[i_sc]/ref_pilots[i_pilot];
                        }
                        // std::cout << "symbol_mse: " << symbol_mse << ", d_symbol_ind: " << d_symbol_ind << std::endl;

                    }
                }
                else if (packet_type_SIG == PACKET_TYPE::DATA)
                {
                    float csi_est;
                    float csi_est_mean = 0;

                    for(int i = 0; i < d_data_carriers.size(); i++) 
                    {
                        int i_sc = d_data_carriers[i];

                        csi_est = real((d_H_mimo[i_sc])*conj(d_H_mimo[i_sc])) + noise_power_sum/snr_est_count;

                        csi_est_mean += csi_est;
                        eq_data_symbols_Z[i] = rx_symbol_Y[i_sc] * conj(d_H_mimo[i_sc]) / csi_est;
                        // eq_data_symbols_Z[i] = rx_symbol_Y[i_sc] / d_H_mimo[i_sc];
                    }

                    if (d_estimator_algo == STA)
                    {
                        uint8_t symbol_value;
                        gr_complex demod_symbol_X;
                        float symbol_mse = 0;
                        float alpha = 0.4;
                        for(int i_data = 0; i_data < d_data_carriers.size(); i_data++) 
                        {
                            int i_sc = d_data_carriers[i_data];
                            symbol_value = modulator_SIG->decision_maker(&eq_data_symbols_Z[i_data]);
                            modulator_SIG->map_to_points(symbol_value, &demod_symbol_X);
                            if(modulator_SIG->bits_per_symbol() == d_qpsk->bits_per_symbol())
                            {
                                demod_symbol_X = demod_symbol_X / (float)2.0;
                            }

                            symbol_mse += std::pow(abs(eq_data_symbols_Z[i_data] - demod_symbol_X),2);

                            d_H_mimo[i_sc] = gr_complex(1 - alpha, 0) * d_H_mimo[i_sc] + gr_complex(alpha, 0) * rx_symbol_Y[i_sc]/demod_symbol_X;
                        }
                        symbol_mse = symbol_mse / d_data_carriers.size();

                        // if(symbol_mse < 0.1)
                        // {
                        //     for(int i_data = 0; i_data < d_data_carriers.size(); i_data++) 
                        //     {
                        //         int i_sc = d_data_carriers[i_data];
                        //         d_H_mimo[i_sc] = d_H_mimo_temp[i_sc];
                        //     }
                        // }
                        // std::cout << "[OFDM Equalizer] symbol_mse: " << symbol_mse << std::endl;

                        for(int i_pilot = 0; i_pilot < d_pilot_carriers.size(); i_pilot++) 
                        {
                            int i_sc = d_pilot_carriers[i_pilot];

                            // symbol_mse += std::pow(abs(eq_data_symbols_Z[i_data] - demod_symbol_X),2);
                            d_H_mimo[i_sc] = gr_complex(1 - alpha, 0) * d_H_mimo[i_sc] + gr_complex(alpha, 0) * rx_symbol_Y[i_sc]/ref_pilots[i_pilot];
                        }

                    }
                }
                else
                {
                    std::cout << "[OFDM Equalizer] Something is wrong --> Packet type is not correct!" << std::endl;
                }
                
                // pmt::pmt_t pdu = pmt::make_dict();
                // message_port_pub(pmt::mp("symbols"), pmt::cons(pmt::make_dict(), pmt::init_c32vector(d_N_data_carr, eq_data_symbols_Z)));
                
                memcpy(out + n_out*d_N_data_carr, eq_data_symbols_Z, d_N_data_carr*sizeof(gr_complex));

                n_out++;
            }
            	
            n_in++;
            d_symbol_ind++;
        }

        total_out += n_out;

        if (n_out != 0)
        {
            dout << "[OFDM Equalizer] n_in: " << n_in << ", n_out: " << n_out <<  ", total_out = " << total_out << std::endl;
        }

        if (total_out == n_ofdm_symbols_SIG && sig_decode_success && !equalize_done)
        {
            if (snr_est_count != 0)
            {
                d_precoded_snr_est = 10 * std::log10( (signal_power_sum/snr_est_count) / (noise_power_sum/snr_est_count));
                dout << "[OFDM Equalizer] Estimated Precoded SNR: " << d_precoded_snr_est << std::endl;
            }

            pmt::pmt_t dict = pmt::make_dict();
            dict = pmt::dict_add(dict, pmt::mp("snr_data"), pmt::from_double(d_precoded_snr_est));
            dict = pmt::dict_add(dict, pmt::mp("chan_mean"), chan_est_mean_pmt);
            add_item_tag(0, nitems_written(0) + n_out-1, pmt::string_to_symbol("stream_end"), dict, pmt::string_to_symbol(alias()));

            equalize_done = true;
        }

        // std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
        // std::chrono::duration<float> diff = now - last_perf_time;

        // if (diff.count() >= perf_display_interval)
        // {
        //     std::cout << "===========================================================" << std::endl;
        //     std::cout << "[OFDM Equalizer] Work Time: " << pc_work_time_avg() << std::endl;
        //     std::cout << "[OFDM Equalizer] Input Buffer: " << pc_input_buffers_full(0) << ",  Output Buffer: "  << pc_output_buffers_full(0) << std::endl;
        //     last_perf_time = std::chrono::system_clock::now();
        // }
 
        consume(0, n_in);

        return n_out;
    }

    bool mimo_ofdm_equalizer_impl::decode_signal_field(gr_complex* rx_symbols)
    {
        
        // FIXME using dynamic array instead of static fixes wrong decoding --> WHY?
        // BUT dynamic array cause memory corruption --> need to free() it 
        // uint8_t rx_bits[d_N_data_carr];
        // uint8_t* rx_bits = (uint8_t*) calloc((d_N_data_carr), sizeof(uint8_t));
        
        for(int i = 0; i < d_N_data_carr; i++) 
        {
            rx_bits[i] = d_bpsk->decision_maker(&rx_symbols[i]);	
        }

        ofdm_mcs ofdm(MCS::BPSK_1_2, d_N_data_carr);
        packet_param frame(ofdm, 0, PACKET_TYPE::NDP);


        uint8_t* decoded_bits = d_decoder.decode(&ofdm, &frame, rx_bits);
                
        int rate_bitmap = 0;
        int packet_type_bitmap = 0;
        data_length_SIG = 0;
        bool parity = false;
        for (int i = 0; i < 17; i++) 
        {
            parity ^= decoded_bits[i];

            if ((i < 4) && decoded_bits[i]) 
            {
                rate_bitmap = rate_bitmap | (1 << i);
            }

            if ((i == 4) && decoded_bits[i]) 
            {
                packet_type_bitmap = packet_type_bitmap | (1 << i-4);
            }

            if (decoded_bits[i] && (i > 4) && (i < 17)) 
            {
                data_length_SIG = data_length_SIG | (1 << (i - 5));
            }
        }

        bool trailing_zeros_correct = true;
        for (int i = 17; i < 23; i++) 
        {
            if(decoded_bits[i] != 0)
            {
                trailing_zeros_correct = false;
            }
        }

        if (parity != decoded_bits[17] && trailing_zeros_correct) 
        {
            std::cerr << "[OFDM Equalizer] SIGNAL Field Check --> FAILED!!" << std::endl;

            packet_type_bitmap = 0;
            data_length_SIG = 0;   
            n_ofdm_symbols_SIG = 0;         
            return false;
        }
        dout << "[OFDM Equalizer] SIGNAL Field Check --> SUCCESS!!" << std::endl;

        switch (packet_type_bitmap) {
                case 0:{ 
                    packet_type_SIG = PACKET_TYPE::NDP;
                    break;
                }
                case 1:{
                    packet_type_SIG = PACKET_TYPE::DATA;
                    break;
                }
                default: {
                    dout << "[OFDM Equalizer] Packet type in SIG is not correct" << std::endl;
                    return false;
                }
        }

        switch (rate_bitmap) {
        case 11:
            mcs_SIG = BPSK_1_2;
            modulator_SIG = d_bpsk;
            // dout << "Encoding: 3 Mbit/s   ";
            break;
        case 15:
            mcs_SIG = BPSK_3_4;
            modulator_SIG = d_bpsk;
            // dout << "Encoding: 4.5 Mbit/s   ";
            break;
        case 10:
            mcs_SIG = QPSK_1_2;
            modulator_SIG = d_qpsk;
            // dout << "Encoding: 6 Mbit/s   ";
            break;
        case 14:
            mcs_SIG = QPSK_3_4;
            modulator_SIG = d_qpsk;
            // dout << "Encoding: 9 Mbit/s   ";
            break;
        case 9:
            mcs_SIG = QAM16_1_2;
            modulator_SIG = d_16qam;
            // dout << "Encoding: 12 Mbit/s   ";
            break;
        case 13:
            mcs_SIG = QAM16_3_4;
            modulator_SIG = d_16qam;
            // dout << "Encoding: 18 Mbit/s   ";
            break;
        // case 8:
        //     mcs_SIG = 6;
        //     n_ofdm_symbols_SIG = (int)ceil((16 + 8 * data_length_SIG + 6) / (double)192);
        //     modulator_SIG = d_64qam;
        //     dout << "Encoding: 24 Mbit/s   ";
        //     break;
        // case 12:
        //     mcs_SIG = 7;
        //     n_ofdm_symbols_SIG = (int)ceil((16 + 8 * data_length_SIG + 6) / (double)216);
        //     modulator_SIG = d_64qam;
        //     dout << "Encoding: 27 Mbit/s   ";
        //     break;
        default:
            dout << "[OFDM Equalizer] modulation in SIG is not correct!" << std::endl;
            return false;
        }

        ofdm_mcs ofdm_param(mcs_SIG, d_N_data_carr);
        packet_param burst(ofdm_param, data_length_SIG, packet_type_SIG);
        n_ofdm_symbols_SIG = burst.n_ofdm_sym;

        return true;
    }

    bool mimo_ofdm_equalizer_impl::parse_signal(uint8_t* decoded_bits)
    {

        int rate_bitmap = 0;
        int packet_type_bitmap = 0;
        data_length_SIG = 0;
        bool parity = false;
        for (int i = 0; i < 17; i++) 
        {
            parity ^= decoded_bits[i];

            if ((i < 4) && decoded_bits[i]) 
            {
                rate_bitmap = rate_bitmap | (1 << i);
            }

            if ((i == 4) && decoded_bits[i]) 
            {
                packet_type_bitmap = packet_type_bitmap | (1 << i-4);
            }

            if (decoded_bits[i] && (i > 4) && (i < 17)) 
            {
                data_length_SIG = data_length_SIG | (1 << (i - 5));
            }
        }

        bool trailing_zeros_correct = true;
        for (int i = 17; i < 23; i++) 
        {
            if(decoded_bits[i] != 0)
            {
                trailing_zeros_correct = false;
            }
        }

        if (parity != decoded_bits[17] && trailing_zeros_correct) 
        {
            std::cerr << "[OFDM Equalizer] SIGNAL Field Check --> FAILED!!" << std::endl;

            packet_type_bitmap = 0;
            data_length_SIG = 0;   
            n_ofdm_symbols_SIG = 0;         
            return false;
        }
        dout << "[OFDM Equalizer] SIGNAL Field Check --> SUCCESS!!" << std::endl;

        switch (packet_type_bitmap) {
                case 0:{ 
                    packet_type_SIG = PACKET_TYPE::NDP;
                    break;
                }
                case 1:{
                    packet_type_SIG = PACKET_TYPE::DATA;
                    break;
                }
                default: {
                    dout << "[OFDM Equalizer] Packet type in SIG is not correct" << std::endl;
                    return false;
                }
        }

        switch (rate_bitmap) {
        case 11:
            mcs_SIG = BPSK_1_2;
            modulator_SIG = d_bpsk;
            // dout << "Encoding: 3 Mbit/s   ";
            break;
        case 15:
            mcs_SIG = BPSK_3_4;
            modulator_SIG = d_bpsk;
            // dout << "Encoding: 4.5 Mbit/s   ";
            break;
        case 10:
            mcs_SIG = QPSK_1_2;
            modulator_SIG = d_qpsk;
            // dout << "Encoding: 6 Mbit/s   ";
            break;
        case 14:
            mcs_SIG = QPSK_3_4;
            modulator_SIG = d_qpsk;
            // dout << "Encoding: 9 Mbit/s   ";
            break;
        case 9:
            mcs_SIG = QAM16_1_2;
            modulator_SIG = d_16qam;
            // dout << "Encoding: 12 Mbit/s   ";
            break;
        case 13:
            mcs_SIG = QAM16_3_4;
            modulator_SIG = d_16qam;
            // dout << "Encoding: 18 Mbit/s   ";
            break;
        // case 8:
        //     mcs_SIG = 6;
        //     n_ofdm_symbols_SIG = (int)ceil((16 + 8 * data_length_SIG + 6) / (double)192);
        //     modulator_SIG = d_64qam;
        //     dout << "Encoding: 24 Mbit/s   ";
        //     break;
        // case 12:
        //     mcs_SIG = 7;
        //     n_ofdm_symbols_SIG = (int)ceil((16 + 8 * data_length_SIG + 6) / (double)216);
        //     modulator_SIG = d_64qam;
        //     dout << "Encoding: 27 Mbit/s   ";
        //     break;
        default:
            dout << "[OFDM Equalizer] modulation in SIG is not correct!" << std::endl;
            return false;
        }

        ofdm_mcs ofdm_param(mcs_SIG, d_N_data_carr);
        packet_param burst(ofdm_param, data_length_SIG, packet_type_SIG);
        n_ofdm_symbols_SIG = burst.n_ofdm_sym;

        return true;
    }

    void mimo_ofdm_equalizer_impl::symbol_equalize (gr_complex *rx_symbols, int symbol_idx, gr_complex *eq_data_symbols, ChannelEstimator algo) 
    {
        for(int i = 0; i < d_data_carriers.size(); i++) 
        {
            eq_data_symbols[i] = rx_symbols[d_data_carriers[i]] / d_H[d_data_carriers[i]]; 
        }
    }

    double mimo_ofdm_equalizer_impl::estimate_residual_cfo(gr_complex* rx_symbols, gr_complex* chan_est, const std::vector<gr_complex>& ref_pilots, std::vector<gr_complex>& est_rx_pilots)
    {
        gr_complex pilot_sum = gr_complex(0.0,0.0);
        // prev_pilot_sum = gr_complex(0.0,0.0);
        for (int i_pilot = 0; i_pilot < d_N_pilot_carr; i_pilot++)
        {
            //https://www.albany.edu/faculty/dsaha/teach/2019Spring_CEN574/slides/08_WLAN.pdf 

            est_rx_pilots[i_pilot] = chan_est[d_pilot_carriers[i_pilot]] * ref_pilots[i_pilot];
            // commonPhaseErrorEstimate --> cpe = rxPilots.*conj(chanEstPilotsR.*refPilotsR) --> exp(-1j*cpe)
            pilot_sum += rx_symbols[d_pilot_carriers[i_pilot]] * conj(est_rx_pilots[i_pilot]);

        }
        return arg(pilot_sum);;
    }

    void
    mimo_ofdm_equalizer_impl::set_estimator(ChannelEstimator algo) 
    {
        gr::thread::scoped_lock lock(d_setlock);
        d_estimator_algo = algo;

        switch(algo) {
            case LS:
                std::cout << "[OFDM Equalizer] Channel Estimator Algorithm: LS " << std::endl;
                break;
            case STA:
                std::cout << "[OFDM Equalizer] Channel Estimator Algorithm: STA " << std::endl;
                break;
            default:
                throw std::runtime_error("[OFDM Equalizer] Estimator not implemented");
        }
    }

    void
    mimo_ofdm_equalizer_impl::set_bandwidth(double bw) 
    {
        gr::thread::scoped_lock lock(d_setlock);
        d_bw = bw;
    }

    void
    mimo_ofdm_equalizer_impl::set_frequency(double freq) 
    {
        gr::thread::scoped_lock lock(d_setlock);
        d_freq = freq;
    }

    void mimo_ofdm_equalizer_impl::set_stats_record(bool stats_record)
    {
        gr::thread::scoped_lock lock(d_setlock);
        d_stats_record = stats_record;
    }

  } /* namespace mimo_ofdm_jrc */
} /* namespace gr */

