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
#include "mimo_precoder_impl.h"
#include <gnuradio/math.h>
#include <mimo_ofdm_jrc/stream_encoder.h>
#include <random> 
#include <volk/volk.h>
#include "viterbi_decoder.h"
#include <math.h>

namespace gr {
  namespace mimo_ofdm_jrc {

    mimo_precoder::sptr
    mimo_precoder::make(int fft_len,
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
                        bool debug)
    {
      return gnuradio::get_initial_sptr
        (new mimo_precoder_impl(fft_len, 
                                N_tx, 
                                N_ss, 
                                data_carriers, 
                                pilot_carriers, 
                                pilot_symbols, 
                                sync_words, 
                                mapped_ltf_symbols,
                                chan_est_file, 
                                chan_est_smoothing,
                                radar_log_file,
                                radar_aided,
                                phased_steering,
                                use_radar_streams,
                                len_tag_key, 
                                debug));
    }


    /*
     * The private constructor
     */
    mimo_precoder_impl::mimo_precoder_impl(int fft_len,
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
                                            bool debug)
      : gr::tagged_stream_block("mimo_precoder",
                gr::io_signature::make(N_ss, N_ss, sizeof(gr_complex)),
                gr::io_signature::make(N_tx, N_tx, fft_len*sizeof(gr_complex)), len_tag_key),
                d_fft_len(fft_len),
                d_N_tx(N_tx),
                d_N_ss(N_ss),
                d_data_carriers(data_carriers),
                d_pilot_carriers(pilot_carriers),
                d_pilot_symbols(pilot_symbols),
                d_sync_words(sync_words),
                d_mapped_ltf_symbols(mapped_ltf_symbols),
                d_chan_est_file(chan_est_file),
                d_radar_log_file(radar_log_file),
                d_chan_est_smoothing(chan_est_smoothing),
                d_radar_aided(radar_aided),
                d_phased_steering(phased_steering),
                d_use_radar_streams(use_radar_streams),
                d_debug(debug)
    {
        perf_display_interval = 1.0;
        last_perf_time = std::chrono::system_clock::now();

        // for mimo-ofdm radar we always use N_tx LTF symbols
        d_N_ltf = d_N_tx;

        // Sanity checks
        // If that is is null, the input is wrong -> force user to use () in python
        if (d_data_carriers.empty()) {
            throw std::invalid_argument(
                "Data carriers must be of type vector of vector i.e. ().");
        }

        for (unsigned j = 0; j < d_data_carriers.size(); j++) {
            if (data_carriers[j] < 0) 
            {
                d_data_carriers[j] += d_fft_len;
            }
            if (d_data_carriers[j] > d_fft_len || d_data_carriers[j] < 0) 
            {
                throw std::invalid_argument("data carrier index out of bounds");
            }
            d_data_carriers[j] = (d_data_carriers[j] + fft_len / 2) % fft_len;
            
        }
        
        if (d_pilot_carriers.empty()) {
            throw std::invalid_argument("Pilot carriers must be of type vector of vector i.e. ((),).");
        }

        for (unsigned j = 0; j < d_pilot_carriers.size(); j++) {
            if (d_pilot_carriers[j] < 0) {
                d_pilot_carriers[j] += d_fft_len;
            }
            if (d_pilot_carriers[j] > d_fft_len || d_pilot_carriers[j] < 0) {
                throw std::invalid_argument("pilot carrier index out of bounds");
            }
            
            d_pilot_carriers[j] = (d_pilot_carriers[j] + fft_len / 2) % fft_len;
        }

        if (d_pilot_symbols.empty()) {
            throw std::invalid_argument("Pilot symbols must be of type vector of vector i.e. ((),).");
        }
        for (unsigned i = 0; i < d_pilot_symbols.size(); i++) 
        {
            if (d_pilot_carriers.size() != d_pilot_symbols[i].size()) 
            {
                throw std::invalid_argument("pilot_carriers do not match pilot_symbols");
            }
        }

        dout << "[MIMO PRECODER] d_sync_words.size: " << sync_words.size() << std::endl;
        dout << "[MIMO PRECODER] mapped_ltf_symbols.size: " << mapped_ltf_symbols.size() << std::endl;
        dout << "[MIMO PRECODER] mapped_ltf_symbols[0].size: " << mapped_ltf_symbols[0].size() << std::endl;

        for (unsigned i = 0; i < d_sync_words.size(); i++) 
        {
            if (d_sync_words[i].size() != (unsigned)d_fft_len) 
            {
                throw std::invalid_argument("[MIMO PRECODER] sync words must be fft length");
            }
        }
        
        
        if (d_mapped_ltf_symbols.size() != d_fft_len) //d_N_tx*d_fft_len
        {
            throw std::invalid_argument("[MIMO PRECODER] MIMO LTF symbols should have (fft length x Ntx) rows!!");
        }

        for (int i_row = 0; i_row < d_mapped_ltf_symbols.size(); i_row++)
        {
            if (d_mapped_ltf_symbols[i_row].size() != d_N_tx*d_N_ltf) //d_N_tx
            {
                throw std::invalid_argument("[MIMO PRECODER] MIMO LTF symbols should have (Ntx) columns!!");
            }
        }
        mapped_ltf_mat.resize(d_fft_len, d_N_tx*d_N_ltf);

        for (int i_row = 0; i_row < d_mapped_ltf_symbols.size(); i_row++){
            mapped_ltf_mat.row(i_row) = Eigen::VectorXcf::Map(&d_mapped_ltf_symbols[i_row][0],d_mapped_ltf_symbols[i_row].size());
        }
        mapped_ltf_mat.resize(d_fft_len, d_N_tx*d_N_ltf);

        for (int i_row = 0; i_row < d_mapped_ltf_symbols.size(); i_row++){
            mapped_ltf_mat.row(i_row) = Eigen::VectorXcf::Map(&d_mapped_ltf_symbols[i_row][0],d_mapped_ltf_symbols[i_row].size());
        }

        d_N_ss_radar = d_N_tx - d_N_ss;
        dout << "[MIMO PRECODER] Number of Data Streams: " << d_N_ss << std::endl;
        dout << "[MIMO PRECODER] Number of Radar Streams: " << d_N_ss_radar << std::endl;

        d_N_data_carriers = d_data_carriers.size();
        dout << "[MIMO PRECODER] Number of Data Subcarriers: " << d_N_data_carriers << std::endl;
        dout << "[MIMO PRECODER] Number of Pilot Subcarriers: " << d_pilot_carriers.size() << std::endl;

        d_active_carriers.reserve(d_data_carriers.size() + d_pilot_carriers.size());
        d_active_carriers.insert(d_active_carriers.end(), d_data_carriers.begin(), d_data_carriers.end());
        d_active_carriers.insert(d_active_carriers.end(), d_pilot_carriers.begin(), d_pilot_carriers.end());
        d_N_active_carriers = d_active_carriers.size();

        signal_field_symbols = new gr_complex[d_N_data_carriers];

        d_qpsk = digital::constellation_qpsk::make();
        d_bpsk = digital::constellation_bpsk::make();

        steering_matrix.resize(d_fft_len);
        for (auto element : steering_matrix)
        {
            element.resize(d_N_tx*d_N_tx);
        }

        steering_matrix_mean.resize(d_N_tx*d_N_tx);
        chan_est_vector_mean.resize(d_N_tx);
        mean_chan_est_changed = true;

        dft_matrix = get_dft_matrix_eigen(d_N_tx);

        std::cout << "[MIMO PRECODER] Channel Estimate File: " << d_chan_est_file << std::endl;

        std::ifstream file_stream(d_chan_est_file);
        if (!file_stream.is_open())
        {
            std::cout << "[MIMO PRECODER] Could not open channel estimate file!" << std::endl;
        }
        file_stream.close();
        
        std::cout << "[MIMO PRECODER] Radar Log File: " << d_radar_log_file << std::endl;
        file_stream.open(d_radar_log_file);

        if (!file_stream.is_open())
        {
            std::cout << "[MIMO PRECODER] Could not open radar log file!" << std::endl;
        }
        file_stream.close();

        last_chanEst_update_time = 0;

        set_tag_propagation_policy(TPP_DONT);

        //FIXME set_relative_rate
        set_relative_rate(1, (uint64_t)d_data_carriers.size());
    }

    /*
     * Our virtual destructor.
     */
    mimo_precoder_impl::~mimo_precoder_impl()
    {
        delete signal_field_symbols;
    }

    int
    mimo_precoder_impl::calculate_output_stream_length(const gr_vector_int &ninput_items)
    {
        int nin = ninput_items[0];
        int nout = (nin / d_data_carriers.size());

        return d_sync_words.size() + 1 + d_N_tx + nout; //SYNC + MIMO_LTF + DATA
    }

    int                                                         // The number of items produced is returned, this can be less than noutput_items
    mimo_precoder_impl::work (int noutput_items,                // total number of items in each output buffer
                       gr_vector_int &ninput_items,             // vector describing the length of each input buffer
                       gr_vector_const_void_star &input_items,  // vector of input buffers, where each element corresponds to an input port
                       gr_vector_void_star &output_items)       // vector of output buffers, where each element corresponds to an output port
    {
        gr::block::set_thread_priority(30);
        gr::thread::scoped_lock lock(d_setlock);

        dout << "================================================================================" << std::endl;        
        dout << "[MIMO PRECODER] ninput_items[0]: " << ninput_items[0] << " noutput_items: " << noutput_items << std::endl;        

        if (output_items.size() != d_N_tx){
            throw std::runtime_error("[MIMO PRECODER] output_items should contain {N_tx} buffers!!");
        }

        int N_ofdm_symbol = ninput_items[0] / d_N_data_carriers;
        dout << "[MIMO PRECODER] OFDM Symbols for Data: " << N_ofdm_symbol << std::endl;

        int N_total_symbols = N_ofdm_symbol + d_sync_words.size() + d_N_tx + 1; //+1 for SIG field
        dout << "[MIMO PRECODER] Total OFDM symbols with preambles and SIG: " << N_total_symbols << std::endl;

        // SETTING POINTERS FOR INPUT/OUTPUT BUFFERS
        const gr_complex* in = (const gr_complex*)input_items[0];

        std::vector<gr_complex *> out_ptrs(output_items.size());
        for (int i_out = 0; i_out < output_items.size(); i_out++){
            out_ptrs[i_out] = (gr_complex *) output_items[i_out];
        }
    
        // PROCESS TAGS (is this fast enough?)
        std::vector<tag_t> tags;
        get_tags_in_range(tags, 0, nitems_read(0), nitems_read(0) + ninput_items[0], pmt::mp("mcs"));
        if (tags.size() != 1) {
            throw std::runtime_error("no mcs tag in input stream!");
        }
        MCS mcs = (MCS)pmt::to_long(tags[0].value);
        dout << "[MIMO PRECODER] MCS from encoder: " << mcs << std::endl;

        get_tags_in_range(tags, 0, nitems_read(0), nitems_read(0) + ninput_items[0], pmt::mp("packet_type"));
        if (tags.size() != 1) {
            throw std::runtime_error("no packet_type tag in input stream!");
        }
        PACKET_TYPE packet_type = (PACKET_TYPE) pmt::to_uint64(tags[0].value);
        dout << "[MIMO PRECODER] Packet Type from encoder: " << packet_type << std::endl;

        get_tags_in_range(tags, 0, nitems_read(0), nitems_read(0) + ninput_items[0], pmt::mp("pdu_len"));
        if (tags.size() != 1) {
            throw std::runtime_error("no pdu_len tag in input stream!");
        }
        int data_size_crc = pmt::to_long(tags[0].value);
        dout << "[MIMO PRECODER] Data size with CRC from encoder: " << data_size_crc << " byte" << std::endl;

        ofdm_mcs ofdm_param(mcs, d_N_data_carriers);
        packet_param frame_param(ofdm_param, data_size_crc, packet_type);

        if (frame_param.n_ofdm_sym != N_ofdm_symbol)
        {
            throw std::runtime_error("[MIMO PRECODER] something is wrong!!");
        }

        // COPYING LEGACY (SYNC) PREAMBLE TO ALL OUTPUTS (NO PRECODING)
        for (int i_out = 0; i_out < out_ptrs.size(); i_out++){
            memset(out_ptrs[i_out], 0, sizeof(gr_complex) * d_fft_len * N_total_symbols);

                // Copy Sync word
                for (unsigned i = 0; i < d_sync_words.size(); i++) 
                {
                    if (i_out < 2)
                    {
                        memcpy(out_ptrs[i_out], &d_sync_words[i][0], sizeof(gr_complex) * d_fft_len);
                    }
                    out_ptrs[i_out] += d_fft_len;
                }


        }

        // gr_complex* signal_field_symbols = (gr_complex *) calloc((d_N_data_carriers), sizeof(gr_complex));
        memset(signal_field_symbols, 0, sizeof(gr_complex) * d_N_data_carriers);
        generate_signal_field(signal_field_symbols, frame_param, ofdm_param);

        // COPYING SIG FIELD TO ALL OUTPUTS (NO PRECODING)
        for (int i_out = 0; i_out < out_ptrs.size(); i_out++)
        {
            if (i_out < 2)
            {
                for (int i_data = 0; i_data < d_N_data_carriers; i_data++)
                {
                    out_ptrs[i_out][d_data_carriers[i_data]] = signal_field_symbols[i_data];
                }
                for (int i_pilot = 0; i_pilot < d_pilot_carriers.size(); i_pilot++)
                {
                    out_ptrs[i_out][d_pilot_carriers[i_pilot]] = d_pilot_symbols[0][i_pilot];
                }
            }
            out_ptrs[i_out] += d_fft_len;
        }
        
        long i_ofdm_symbol = 0; // Number of output items
        if (packet_type == PACKET_TYPE::NDP)
        {
            dout << "[MIMO PRECODER] Generating an NDP frame..." << std::endl;
            
            // COPYING MIMO PREAMBLE 
            for (int i_tx = 0; i_tx < out_ptrs.size(); i_tx++)
            {
                for (int i_sc = 0; i_sc < d_fft_len; i_sc++){
                    for (int i_ltf = 0; i_ltf < d_N_ltf; i_ltf++)
                    {
                        out_ptrs[i_tx][i_ltf*d_fft_len + i_sc] = d_mapped_ltf_symbols[i_sc][i_ltf + i_tx*d_N_ltf];
                    }
                }
                out_ptrs[i_tx] += d_N_ltf*d_fft_len;
            }
            

            //TODO Precode NDP? (fourier?)
            // COPYING NDP 
            int symbols_placed = 0;
            for (int i_in = 0; i_in < ninput_items[0]; i_in++) 
            {
                if (symbols_placed == 0) 
                {
                    i_ofdm_symbol++;
                }

                for (int i_out = 0; i_out < out_ptrs.size(); i_out++)
                {
                    if (i_out < 2)
                    {
                        out_ptrs[i_out][(i_ofdm_symbol-1)*d_fft_len + d_data_carriers[symbols_placed]] = in[i_in];
                    }
                }

                symbols_placed++;
                if (symbols_placed == d_data_carriers.size()) {
                    symbols_placed = 0;
                }
            }

            // Copy pilot symbols
            for (int i_sym = 0; i_sym < i_ofdm_symbol; i_sym++) 
            {
                for (unsigned i_pilot = 0; i_pilot < d_pilot_carriers.size(); i_pilot++) 
                {
                    for (int i_out = 0; i_out < out_ptrs.size(); i_out++)
                    {
                        if (i_out < 2)
                        {
                        out_ptrs[i_out][i_sym * d_fft_len + d_pilot_carriers[i_pilot]] = d_pilot_symbols[i_sym % d_pilot_symbols.size()][i_pilot];
                        }
                    }
                }
            }
        }
        else if (packet_type == PACKET_TYPE::DATA)
        {
            dout << "[MIMO PRECODER] Generating an DATA frame..." << std::endl;
            
            // GENERATE RADAR STREAMS
            std::random_device rd;
            std::default_random_engine e1(rd());
            std::uniform_int_distribution<unsigned int> uniform_dist( 0, 3 );

            int symbols_generated = 0;

            // gr_complex stream_symbols [N_ofdm_symbol*d_fft_len] [d_N_tx];
            // stream_symbols

            Eigen::MatrixXcf stream_symbols;
            int N_streams_touse ;
            if (d_use_radar_streams)
            {
                N_streams_touse = d_N_tx;
            }
            else
            {
                N_streams_touse = 1;
            }

            stream_symbols.resize(N_streams_touse, N_ofdm_symbol*d_fft_len);

            for (int i_symbol = 0; i_symbol < N_ofdm_symbol; i_symbol++)
            {   
                for (int i_ss = 0; i_ss < N_streams_touse; i_ss++)
                {
                    for (int i_data= 0; i_data < d_N_data_carriers; i_data++)
                    {
                        int i_sc = d_data_carriers[i_data];

                        if (i_ss == 0)
                        {
                            stream_symbols(i_ss, i_symbol*d_fft_len + i_sc)  = in[i_symbol*d_N_data_carriers + i_data];
                        }
                        else
                        {
                            d_qpsk->map_to_points(uniform_dist(e1), &stream_symbols(i_ss, i_symbol*d_fft_len + i_sc));  
                            stream_symbols(i_ss, i_symbol*d_fft_len + i_sc) = stream_symbols(i_ss, i_symbol*d_fft_len + i_sc) / (float)2.0;                  
                        }
                    }

                    for (int i_pilot = 0; i_pilot < d_pilot_carriers.size(); i_pilot++)
                    {
                        int i_sc = d_pilot_carriers[i_pilot];
                        // stream_symbols[i_symbol*d_N_active_carriers + i_carrier] = new gr_complex[d_N_ss_radar] ();
                        if (i_ss == 0)
                        {
                            // stream_symbols[i_symbol*d_fft_len + i_sc][i_ss] =   d_pilot_symbols[i_symbol % d_pilot_symbols.size()][i_pilot];
                            stream_symbols(i_ss, i_symbol*d_fft_len + i_sc) = d_pilot_symbols[i_symbol % d_pilot_symbols.size()][i_pilot];
                        }
                        else
                        {
                            // d_bpsk->map_to_points(uniform_dist(e1), &stream_symbols[i_symbol*d_fft_len + i_sc][i_ss]);
                            d_qpsk->map_to_points(uniform_dist(e1), &stream_symbols(i_ss, i_symbol*d_fft_len + i_sc));  
                            stream_symbols(i_ss, i_symbol*d_fft_len + i_sc) = stream_symbols(i_ss, i_symbol*d_fft_len + i_sc) / (float)2.0;          
                        }
                    }
                }
            }
        
            std::vector<gr_complex> chan_est_vector_mean_norm(d_N_tx);

            bool use_fourier_precoding = false;
            bool use_custom_precoding = false;

            if (d_radar_aided)
            {
                use_custom_precoding = compute_radar_aided_steering();
            }
            else
            {
                use_custom_precoding = compute_steering_matrix();
            }

            if (use_custom_precoding)
            {
                double angle_estimate;
                dout << "[MIMO PRECODER] Mean channel vector: " << std::endl;
                for (int i_tx = 0; i_tx < d_N_tx; i_tx++ )
                {
                    // dout << "[MIMO PRECODER] Estimated mean channel vector: " << chan_est_vector_mean[i_tx] << std::endl;
                    chan_est_vector_mean_norm[i_tx] = chan_est_vector_mean[i_tx] / chan_est_vector_mean[0];
                    
                    dout << chan_est_vector_mean[i_tx] << ", ";

                    if (i_tx == 1)
                    {
                        angle_estimate = asin( arg( chan_est_vector_mean_norm[i_tx])/M_PI )*180/M_PI;
                    }
                }
                dout << std::endl;
                dout << "[MIMO PRECODER] Angle estimate: " << angle_estimate << std::endl;
            }
            else
            {
                use_fourier_precoding = true;
                dout << "[MIMO PRECODER] Custom Precoding Failed --> USING FOURIER PRECODING!" << std::endl;
            }

            //= PRECODING MIMO PREAMBLE ===============================================================
            // TODO reduce for loops
            Eigen::MatrixXcf mimo_ltf_precoded(d_N_tx, d_N_ltf);
            Eigen::Map<const Eigen::Matrix <gr_complex, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > X_ltf( &d_mapped_ltf_symbols[0][0], d_N_tx, d_N_ltf );
            Eigen::Map<Eigen::MatrixXcf > Q_steer(&steering_matrix_mean[0], d_N_tx, d_N_tx);

            for (int i_sc = 0; i_sc < d_fft_len; i_sc++)
            {
                new (&X_ltf) Eigen::Map<const Eigen::Matrix <gr_complex, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> ( &d_mapped_ltf_symbols[i_sc][0], d_N_tx, d_N_ltf );

                if (X_ltf.isZero())
                {
                    mimo_ltf_precoded = Eigen::MatrixXcf::Zero(d_N_tx, d_N_ltf);
                }
                else
                {
                    if (use_fourier_precoding)
                    {
                        mimo_ltf_precoded = dft_matrix*X_ltf;
                    }
                    else
                    {
                        // TODO improve performance of smoothed steering!! https://eigen.tuxfamily.org/dox-devel/group__TutorialMapClass.html
                        if(!d_chan_est_smoothing && !d_radar_aided)
                        {
                            new (&Q_steer) Eigen::Map<Eigen::MatrixXcf> (&steering_matrix[i_sc][0], d_N_tx, d_N_tx);
                        }
                        mimo_ltf_precoded = Q_steer*X_ltf;
                    }
                }



                for (int i_tx = 0; i_tx < out_ptrs.size(); i_tx++)
                {
                    for (int i_ltf = 0; i_ltf < d_N_ltf; i_ltf++)
                    {
                        // memcpy(out_ptrs[i_out], &d_mapped_ltf_symbols[i_ltf][i_out], sizeof(gr_complex));
                        // out_ptrs[i_tx][i_ltf*d_fft_len + i_sc] = d_mapped_ltf_symbols[i_sc][i_ltf + i_tx*d_N_ltf];
                        out_ptrs[i_tx][i_ltf*d_fft_len + i_sc] = mimo_ltf_precoded(i_tx, i_ltf);
                    }
                }
            }

            for (int i_tx = 0; i_tx < out_ptrs.size(); i_tx++)
            {
                out_ptrs[i_tx] += d_N_ltf*d_fft_len;
            }

            //= PRECODING DATA ==================================================================
            // Eigen::MatrixXcf symbols_precoded(d_N_tx, ninput_items[0]);
            // Eigen::Map<const Eigen::VectorXcf> x_tx(in, ninput_items[0]);     
            // // symbols_precoded = dft_matrix*(x_tx.transpose().replicate( dft_matrix.rows(), 1));
            // symbols_precoded = dft_matrix.col(0)*x_tx.transpose();

            Eigen::MatrixXcf symbols_precoded;
            Eigen::VectorXcf pilots_precoded(d_N_tx);
            // Eigen::VectorXcf in_streams(d_N_tx);
            // Eigen::Map<Eigen::VectorXcf > in_streams(&stream_symbols[0][0], d_N_tx);

            bool use_subcarrier_precoding = !use_fourier_precoding && !d_chan_est_smoothing && !d_radar_aided;

            if (use_subcarrier_precoding)
            {
                // each symbol in subcarriers are precoded with separate steering matrices
                symbols_precoded.resize(d_N_tx,1);
            }
            else
            {
                // symbols_precoded.resize(d_N_tx, ninput_items[0]);
                symbols_precoded.resize(d_N_tx, N_ofdm_symbol*d_fft_len);

                if (use_fourier_precoding)
                {
                    // Eigen::Map<const Eigen::VectorXcf> x_tx(in, ninput_items[0]); 
                    // symbols_precoded = dft_matrix*(x_tx.transpose().replicate( dft_matrix.rows(), 1));
                    // symbols_precoded = dft_matrix.col(0)*x_tx.transpose();
                    
                    if (d_use_radar_streams)
                    {
                        symbols_precoded = dft_matrix*stream_symbols;
                    }
                    else
                    {
                        symbols_precoded = dft_matrix.col(0)*stream_symbols;
                    }
                }
                else
                {
                    new (&Q_steer) Eigen::Map<Eigen::MatrixXcf> (&steering_matrix_mean[0], d_N_tx, d_N_tx);
                    
                    if (d_use_radar_streams)
                    {
                        symbols_precoded = Q_steer*stream_symbols;
                    }
                    else
                    {
                        symbols_precoded = Q_steer.col(0)*stream_symbols;
                    }
                }
            }


            // COPYING PRECODED DATA & PILOTS
            int i_sc_symbol = 0;

            for (int i_in = 0; i_in < ninput_items[0]; i_in++) 
            {
                int i_sc = d_data_carriers[i_sc_symbol];

                if (!use_subcarrier_precoding)
                { 
                    for (int i_out = 0; i_out < out_ptrs.size(); i_out++)
                    {
                        out_ptrs[i_out][i_ofdm_symbol*d_fft_len + i_sc] = symbols_precoded(i_out, i_ofdm_symbol*d_fft_len + i_sc);
                    }
                }
                else
                {
                    new (&Q_steer) Eigen::Map<Eigen::MatrixXcf >(&steering_matrix[i_sc][0], d_N_tx, d_N_tx);                    

                    if (d_use_radar_streams)
                    {
                        // new (&in_streams) Eigen::Map<Eigen::VectorXcf > (&stream_symbols[i_ofdm_symbol*d_fft_len + i_sc][0], d_N_tx);
                        symbols_precoded = Q_steer*stream_symbols.col(i_ofdm_symbol*d_fft_len + i_sc);
                    }
                    else
                    {
                        symbols_precoded = Q_steer.col(0)*in[i_in];
                    }

                    for (int i_out = 0; i_out < out_ptrs.size(); i_out++)
                    {
                        out_ptrs[i_out][i_ofdm_symbol*d_fft_len + i_sc] = symbols_precoded(i_out);//symbols_precoded(i_out, i_in); //in[i_in];
                    }
                }

                i_sc_symbol++;

                if (i_sc_symbol == d_N_data_carriers) 
                {
                    //Now it is time to COPY PILOT symbols to current ofdm symbol
                    for (int i_pilot = 0; i_pilot < d_pilot_carriers.size(); i_pilot++) 
                    {
                        int i_sc = d_pilot_carriers[i_pilot];

                        if (!use_subcarrier_precoding)
                        { 
                            for (int i_out = 0; i_out < out_ptrs.size(); i_out++)
                            {
                                out_ptrs[i_out][(i_ofdm_symbol)*d_fft_len + i_sc] = symbols_precoded(i_out, i_ofdm_symbol*d_fft_len + i_sc);
                            }
                        }
                        else // precoding subcarriers separately
                        {
                            new (&Q_steer) Eigen::Map<Eigen::MatrixXcf >(&steering_matrix[i_sc][0], d_N_tx, d_N_tx);                    

                            if (d_use_radar_streams)
                            {
                                // new (&in_streams) Eigen::Map<Eigen::VectorXcf > (&stream_symbols[i_ofdm_symbol*d_fft_len + i_sc][0], d_N_tx);
                                symbols_precoded = Q_steer*stream_symbols.col(i_ofdm_symbol*d_fft_len + i_sc);
                            }
                            else
                            {
                                symbols_precoded = Q_steer.col(0)*stream_symbols(0, i_ofdm_symbol*d_fft_len + i_sc);
                            }

                            for (int i_out = 0; i_out < out_ptrs.size(); i_out++)
                            {
                                out_ptrs[i_out][i_ofdm_symbol*d_fft_len + i_sc] = symbols_precoded(i_out);//symbols_precoded(i_out, i_in); //in[i_in];
                            }
                        }

                    }
                    i_sc_symbol = 0;
                    i_ofdm_symbol++;

                }
            }
            
        }
        else
        {
            std::cout << "[MIMO PRECODER] packet_type: " << packet_type << std::endl;
            throw std::invalid_argument("[MIMO PRECODER] packet type is not defined!");
        }
        
        dout << "[MIMO PRECODER] OFDM symbols generated: " << i_ofdm_symbol << std::endl;

        if (i_ofdm_symbol != N_ofdm_symbol)
        {
            throw std::runtime_error("[MIMO PRECODER] something is wrong!!");
        }

        // std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
        // std::chrono::duration<float> diff = now - last_perf_time;

        // if (diff.count() >= perf_display_interval)
        // {
        //     std::cout << "===========================================================" << std::endl;
        //     std::cout << "[MIMO PRECODER] Work Time: " << pc_work_time_avg() <<  ", Throughput: " <<  pc_throughput_avg() << std::endl;
        //     std::cout << "[MIMO PRECODER] Input Buffer: " << pc_input_buffers_full_avg(0) << ",  Output Buffer: "  << pc_output_buffers_full_avg(0) << std::endl;
        //     last_perf_time = std::chrono::system_clock::now();
        // }

        // Don't call consume!
        return i_ofdm_symbol + d_sync_words.size() + d_N_tx + 1;
    }

    std::vector<std::vector<gr_complex>> mimo_precoder_impl::get_dft_matrix(int N)
    {
        std::vector<std::vector<gr_complex>> dft_mat;
        dft_mat.resize(N);

        for (int i_row = 0; i_row < N; i_row++)
        {
            dft_mat[i_row].resize(N);
            for(int i_col = 0; i_col < N; i_col++)
            {
                gr_complex exp_part = gr_complex(0, -2 * GR_M_PI * float(i_row*i_col)/float(N));
                dft_mat[i_row][i_col] = std::exp(exp_part)/(gr_complex)std::sqrt(N);
            }
        }

        return dft_mat;
    }

    Eigen::MatrixXcf mimo_precoder_impl::get_dft_matrix_eigen(int N)
    {
        Eigen::MatrixXcf dft_mat(N, N);

        for (int i_row = 0; i_row < N; i_row++){
            for(int i_col = 0; i_col < N; i_col++){
                gr_complex exp_part = gr_complex(0, -2 * GR_M_PI * float(i_row*i_col)/float(N));
                dft_mat(i_row,i_col) = std::exp(exp_part)/(gr_complex)std::sqrt(N);
            }
        }
        return dft_mat;
    }

    // this function may fail if there is an error with reading or parsing the file 
    bool mimo_precoder_impl::compute_steering_matrix()
    {
        // https://stackoverflow.com/questions/34247057/how-to-read-csv-file-and-assign-to-eigen-matrix

        std::ifstream file_stream(d_chan_est_file);
        std::string line;

        if (!file_stream.is_open())
        {
            if (!d_chan_est_file.empty())
            {
                std::cerr << "[MIMO PRECODER] Could not open channel estimate file at " << d_chan_est_file << std::endl;
            }
            return false;
        }

        std::time_t curr_write_time = boost::filesystem::last_write_time(d_chan_est_file);
        dout << "[MIMO PRECODER] Last write time of the chanEst file: " << last_chanEst_update_time << ", Current write time of the file: " << curr_write_time << std::endl;
        
        if(curr_write_time <= last_chanEst_update_time && !mean_chan_est_changed)
        {
            dout << "[MIMO PRECODER] Steering matrix is up-to-date --> No need to compute steering matrix again " << std::endl;
            return true;
        }


        std::fill(chan_est_vector_mean.begin(), chan_est_vector_mean.end(), gr_complex(0,0));

        Eigen::Map<Eigen::VectorXcf> H_est_n(nullptr, d_N_tx);

        int n_line_read = 0;
        while (std::getline(file_stream, line))
        {
            std::stringstream lineStream(line);
            std::vector<gr_complex> chan_est_vector(d_N_tx);

            std::string line_entry;

            getline(lineStream, line_entry, ':');
            int sc_idx = std::stoi(line_entry);
            // dout << sc_idx << ": ";

            float real_part, imag_part;
            int n_col_read = 0;
            while (getline(lineStream, line_entry, ';')) 
            {
                sscanf(line_entry.c_str(), "(%f,%f)", &real_part, &imag_part);
                
                if (std::find(d_active_carriers.begin(), d_active_carriers.end(), n_line_read) != d_active_carriers.end() )
                {
                    chan_est_vector_mean[n_col_read] += gr_complex(real_part, imag_part);
                }

                // chan_est_vector.push_back(gr_complex(real_part, imag_part));
                chan_est_vector[n_col_read] = gr_complex(real_part, imag_part);

                n_col_read++;
                // std::cout << real_part << ", " << imag_part << "; ";
            }
            // std::cout << std::endl;

            n_line_read++;
            
            if(n_col_read != d_N_tx)
            {
                std::cerr << "[MIMO PRECODER] Steering matrix computation FAILED! --> Line is not correct: " << n_line_read << std::endl;
                std::cerr << "[MIMO PRECODER] Channel estimate file at " << d_chan_est_file << std::endl;

                return false;
            }

            new (&H_est_n) Eigen::Map<Eigen::VectorXcf >(&chan_est_vector[0], chan_est_vector.size());

            if(d_phased_steering)
            {
                Eigen::MatrixXcf Q = Eigen::MatrixXcf::Zero(d_N_tx, d_N_tx);
                Q.col(0) = H_est_n.conjugate();
                Q = Q * sqrt(d_N_tx) / Q.norm();
                steering_matrix[sc_idx] = std::vector<gr_complex> (Q.data(), Q.data() + Q.size());
            }
            else
            {
                Eigen::JacobiSVD<Eigen::MatrixXcf> SVD_out(H_est_n.transpose(), Eigen::ComputeFullV);
                // normQ = Q * sqrt(Ntx)/norm(Q, 'fro'); % Normalization
                Eigen::MatrixXcf Q = SVD_out.matrixV() * sqrt(d_N_tx) / SVD_out.matrixV().norm();
                steering_matrix[sc_idx] = std::vector<gr_complex> (Q.data(), Q.data() + Q.size());
            }
        }

        if(n_line_read < d_fft_len)
        {
            std::cerr << "[MIMO PRECODER] Steering matrix computation FAILED! --> Number of parsed lines not correct: " << n_line_read << std::endl;
            std::cerr << "[MIMO PRECODER] Channel estimate file at " << d_chan_est_file << std::endl;

            return false;
        }

        for (int i_tx = 0; i_tx < d_N_tx; i_tx++ )
        {
            chan_est_vector_mean[i_tx] = chan_est_vector_mean[i_tx] / (gr_complex) d_N_active_carriers;
        }

        // Eigen::Map<Eigen::VectorXcf> H_est_mean(&chan_est_vector_mean[0], chan_est_vector_mean.size());
        new (&H_est_n) Eigen::Map<Eigen::VectorXcf >(&chan_est_vector_mean[0], chan_est_vector_mean.size());

        if(d_phased_steering)
        {
            Eigen::MatrixXcf Q = Eigen::MatrixXcf::Zero(d_N_tx, d_N_tx);
            Q.col(0) = H_est_n.conjugate();
            Q = Q * sqrt(d_N_tx) / Q.norm();
            steering_matrix_mean = std::vector<gr_complex> (Q.data(), Q.data() + Q.size());
        }
        else
        {
            Eigen::JacobiSVD<Eigen::MatrixXcf> SVD_out(H_est_n.transpose(), Eigen::ComputeFullV);
            // normQ = Q * sqrt(Ntx)/norm(Q, 'fro'); % Normalization
            Eigen::MatrixXcf Q = SVD_out.matrixV() * sqrt(d_N_tx) / SVD_out.matrixV().norm();
            steering_matrix_mean = std::vector<gr_complex> (Q.data(), Q.data() + Q.size());
        }

        mean_chan_est_changed = false;
        last_chanEst_update_time = curr_write_time;
        return true;
    }



    // this function may fail if there is an error with reading or parsing the file 
    bool mimo_precoder_impl::compute_radar_aided_steering()
    {
        std::ifstream file_stream(d_radar_log_file);

        if (!file_stream.is_open())
        {
            std::cerr << "[MIMO PRECODER] Could not open radar log file at " << d_radar_log_file << std::endl;
            return false;
        }

        file_stream.seekg(-1, std::ios_base::end);
        if(file_stream.peek() == '\n')
        {
            file_stream.seekg(-1, std::ios_base::cur);
            int i = file_stream.tellg();
            for(i; i > 0; i--)
            {
                if(file_stream.peek() == '\n')
                {
                    //Found
                    file_stream.get();
                    break;
                }
                file_stream.seekg(i, std::ios_base::beg);
            }
        }
        std::string lastline;
        std::getline(file_stream, lastline);

        dout << "[MIMO PRECODER] lastline: " << lastline << std::endl;

        if (lastline.empty())
        {
            std::cerr << "[MIMO PRECODER] Radar log file is empty at " << d_radar_log_file << std::endl;
            std::cerr << "[MIMO PRECODER] Radar log file's lastline: " << lastline << std::endl;
            std::cerr << "[MIMO PRECODER] will use channel estimate file" << std::endl;

            compute_steering_matrix();
            return false;
        }
        std::stringstream lineStream(lastline);
        std::string line_entry;

        getline(lineStream, line_entry, ',');
        getline(lineStream, line_entry, ',');
        getline(lineStream, line_entry, ',');
        getline(lineStream, line_entry, ',');
        getline(lineStream, line_entry, '\n');

        angle_estimate = std::stof(line_entry);
        dout << "[MIMO PRECODER] line_entry: " << line_entry << ", angle_estimate:" << angle_estimate <<  std::endl;
        
        // std::vector<gr_complex> chan_est_vector(d_N_tx);
        for (int i_tx = 0; i_tx < d_N_tx; i_tx++)
        {
            chan_est_vector_mean[i_tx] = std::exp(gr_complex(0, GR_M_PI * sin(angle_estimate/180.0*GR_M_PI) * i_tx));
            // dout << chan_est_vector_mean[i_tx] << ", ";
        }
        // dout <<  std::endl;   

        Eigen::Map<Eigen::VectorXcf> H_est_n(&chan_est_vector_mean[0], chan_est_vector_mean.size());

        if(d_phased_steering)
        {
            Eigen::MatrixXcf Q = Eigen::MatrixXcf::Zero(d_N_tx, d_N_tx);
            Q.col(0) = H_est_n.conjugate();
            Q = Q * sqrt(d_N_tx) / Q.norm();
            steering_matrix_mean = std::vector<gr_complex> (Q.data(), Q.data() + Q.size());
        }
        else
        {
            Eigen::JacobiSVD<Eigen::MatrixXcf> SVD_out(H_est_n.transpose(), Eigen::ComputeFullV);
            // normQ = Q * sqrt(Ntx)/norm(Q, 'fro'); % Normalization
            Eigen::MatrixXcf Q = SVD_out.matrixV() * sqrt(d_N_tx) / SVD_out.matrixV().norm();
            steering_matrix_mean = std::vector<gr_complex> (Q.data(), Q.data() + Q.size());
        }

        mean_chan_est_changed = true;

        return true;
    }

    void mimo_precoder_impl::generate_signal_field(gr_complex* out, packet_param& frame, ofdm_mcs& ofdm)
    {
        // data bits of the signal header
        char* signal_header = (char*)calloc((ofdm.d_n_data_carriers/2), sizeof(char));
        // signal header after...
        // convolutional encoding
        char* encoded_signal_header = (char*)calloc((ofdm.d_n_data_carriers), sizeof(char));

        // interleaving
        // char* interleaved_signal_header = (char*)calloc((ofdm.d_n_data_carriers), sizeof(char));

        // splitted symbols per subcarrier
        char* symbols = (char*)calloc((ofdm.d_n_data_carriers), sizeof(char));
        // char *symbols = (char*)calloc((frame.n_encoded_bits / ofdm.n_bpsc), sizeof(char));

        int length = frame.data_size_byte;
        
        // first 4 bits represent the modulation and coding scheme
        signal_header[0] = get_bit(ofdm.rate_field, 3);
        signal_header[1] = get_bit(ofdm.rate_field, 2);
        signal_header[2] = get_bit(ofdm.rate_field, 1);
        signal_header[3] = get_bit(ofdm.rate_field, 0);
        // 5th bit is reserved for packet_type --> two packet types --> 1-bit
        signal_header[4] = get_bit(frame.packet_type_field, 0);
        // then 12 bits represent the length
        signal_header[5] = get_bit(length, 0);
        signal_header[6] = get_bit(length, 1);
        signal_header[7] = get_bit(length, 2);
        signal_header[8] = get_bit(length, 3);
        signal_header[9] = get_bit(length, 4);
        signal_header[10] = get_bit(length, 5);
        signal_header[11] = get_bit(length, 6);
        signal_header[12] = get_bit(length, 7);
        signal_header[13] = get_bit(length, 8);
        signal_header[14] = get_bit(length, 9);
        signal_header[15] = get_bit(length, 10);
        signal_header[16] = get_bit(length, 11);
        // 18-th bit is the parity bit for the first 17 bits
        int sum = 0;
        for (int i = 0; i < 17; i++) {
            if (signal_header[i]) {
                sum++;
            }
        }
        signal_header[17] = sum % 2;
        
        // last 6 bits must be set to 0 for convolutional encoder
        for (int i = 0; i < 6; i++) {
            signal_header[18 + i] = 0;
        }
        
        // for (int i = 0; i < (ofdm.d_n_data_carriers/2); i++) {
        //     dout << (int) signal_header[i] << ", ";
        // }
        // dout<<std::endl;
        
        ofdm_mcs signal_ofdm(MCS::BPSK_1_2, d_N_data_carriers);
        packet_param signal_param(signal_ofdm, 0, frame.packet_type);

        // convolutional encoding (scrambling is not needed)
        convolutional_encoding(signal_header, encoded_signal_header, signal_param);

        // interleaving
        // interleave(encoded_signal_header, out, signal_param, signal_ofdm);
        split_symbols(encoded_signal_header, symbols, signal_param, signal_ofdm);

        for (int i = 0; i < ofdm.d_n_data_carriers; i++)
        {
            d_bpsk->map_to_points(symbols[i], &out[i]);
        }
    
        free(signal_header);
        free(encoded_signal_header);
        // free(interleaved_signal_header);
        free(symbols);
    }       

    int mimo_precoder_impl::get_bit(int b, int i) 
    { 
        return (b & (1 << i) ? 1 : 0); 
    }

    void mimo_precoder_impl::set_chan_est_smoothing(bool chan_est_smoothing)
    {
        gr::thread::scoped_lock lock(d_setlock);
        
        d_chan_est_smoothing = chan_est_smoothing;
    }

    void mimo_precoder_impl::set_radar_aided(bool radar_aided)
    {
        gr::thread::scoped_lock lock(d_setlock);
        
        d_radar_aided = radar_aided;
    }


    void mimo_precoder_impl::set_use_radar_streams(bool use_radar_streams)
    {
        gr::thread::scoped_lock lock(d_setlock);
        
        d_use_radar_streams = use_radar_streams;
    }

    void mimo_precoder_impl::set_phased_steering(bool phased_steering)
    {
        gr::thread::scoped_lock lock(d_setlock);
        
        mean_chan_est_changed = true;
        d_phased_steering = phased_steering;
    }

  } /* namespace mimo_ofdm_jrc */
} /* namespace gr */

