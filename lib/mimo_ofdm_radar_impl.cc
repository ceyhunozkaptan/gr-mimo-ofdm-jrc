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
#include "mimo_ofdm_radar_impl.h"

namespace gr {
namespace mimo_ofdm_jrc {

    mimo_ofdm_radar::sptr
    mimo_ofdm_radar::make(int fft_len,
                            int N_tx,
                            int N_rx,
                            int N_sym,
                            int N_pre,
                            bool background_removal,
                            bool background_recording,
                            int record_len,
                            int interp_factor,
                            bool enable_tx_interleave,
                            const std::string& radar_chan_file,
                            const std::string& len_tag_key,
                            bool debug)
    {
      return gnuradio::get_initial_sptr
        (new mimo_ofdm_radar_impl(fft_len,
                                            N_tx,
                                            N_rx,
                                            N_sym,
                                            N_pre,
                                            background_removal,
                                            background_recording,
                                            record_len,
                                            interp_factor,
                                            enable_tx_interleave,
                                            radar_chan_file,
                                            len_tag_key,
                                            debug));
    }


    /*
     * The private constructor
     */
    mimo_ofdm_radar_impl::mimo_ofdm_radar_impl(
                int fft_len,
				int N_tx,
				int N_rx,
				int N_sym,
                int N_pre,
                bool background_removal,
                bool background_recording,
                int record_len,
                int interp_factor,
                bool enable_tx_interleave,
                const std::string& radar_chan_file,
				const std::string& len_tag_key,
                bool debug
    )
      : gr::block("mimo_ofdm_radar",
                gr::io_signature::make(N_tx+N_rx, N_tx+N_rx, sizeof(gr_complex) * fft_len),
                gr::io_signature::make(1, 1, sizeof(gr_complex) * fft_len * interp_factor)),
                    d_fft_len(fft_len),
                    d_N_tx(N_tx),
                    d_N_rx(N_rx),
                    d_N_sym(N_sym),
                    d_N_pre(N_pre),
                    d_background_removal(background_removal),
                    d_background_recording(background_recording),
                    d_record_len(record_len),
                    d_interp_factor(interp_factor),
                    d_enable_tx_interleave(enable_tx_interleave),
                    d_radar_chan_file(radar_chan_file),
                    d_debug(debug)
    {
        //TODO Add sanity checks!

        // radar_chan_est =[ [fft_len], [fft_len], [fft_len], ...]
        radar_chan_est = new gr_complex[d_N_tx*d_N_rx*d_fft_len];
        
        rx_radar_frame.resize(d_N_rx);
        for (size_t i_rx = 0; i_rx < d_N_rx; i_rx++)
        {
            rx_radar_frame[i_rx].resize(d_fft_len*d_N_sym);
        }


        tx_radar_frame.resize(d_N_tx);
        for (size_t i_tx = 0; i_tx < d_N_tx; i_tx++)
        {
            tx_radar_frame[i_tx].resize(d_fft_len*d_N_sym);
        }

        radar_chan_est_buffer.set_capacity(d_record_len);
        radar_chan_est_temp.resize(d_N_tx*d_N_rx*d_fft_len);

        captured = false;
        set_tag_propagation_policy(TPP_DONT); // does not apply on stream tags!
    }

    /*
     * Our virtual destructor.
     */
    mimo_ofdm_radar_impl::~mimo_ofdm_radar_impl()
    {
        
    }


    int
    mimo_ofdm_radar_impl::general_work(int noutput_items,
                            gr_vector_int &ninput_items,
                            gr_vector_const_void_star &input_items,
                            gr_vector_void_star &output_items)
    {
        const gr_complex *in_tx;
        const gr_complex *in_rx;

        // const gr_complex* in = reinterpret_cast<const gr_complex*>(input_items[0]);
        gr_complex* out = (gr_complex *) output_items[0];

        
        dout << "================================================================================" << std::endl;        
        dout << "[MIMO OFDM RADAR]: input_items.size: " << input_items.size() << std::endl;
        dout << "[MIMO OFDM RADAR]: ninput_items tx0: " << ninput_items[0] << std::endl;
        dout << "[MIMO OFDM RADAR]: ninput_items rx0: " << ninput_items[d_N_tx] << std::endl;
        dout << "[MIMO OFDM RADAR]: nitems_read tx0: " << nitems_read(0) << std::endl;
        dout << "[MIMO OFDM RADAR]: nitems_read rx0: " << nitems_read(d_N_tx) << std::endl;
        dout << "[MIMO OFDM RADAR]: noutput_items: " << noutput_items << std::endl;

        
        std::vector<gr::tag_t> 	rx_tags;
        std::vector<gr::tag_t> 	tx_tags;

        uint start_tag_offset;
        uint rx_packet_len;
        uint tx_packet_len;

        // TODO get tags in range (should be at 0 offet since USRP block generates tagged_streams ) --> start copying all until enough sample collected --> perform matrix computations

        int n_tx_tags;
        uint n_tx_samples_discard = 0;

        int tx_frame_offset = 0;

        get_tags_in_range(rx_tags, d_N_tx, nitems_read(d_N_tx), nitems_read(d_N_tx) + ninput_items[d_N_tx], pmt::mp("packet_len"));
        dout << "[MIMO OFDM RADAR] number of rx_tags:" << rx_tags.size() << std::endl;

        if (rx_tags.size() == 0) 
        {
            dout << "[MIMO OFDM RADAR] no packet_len tag on RX input"<< std::endl;
        }
        else
        {
            get_tags_in_range(tx_tags, 0, nitems_read(0), nitems_read(0) + ninput_items[0], pmt::mp("packet_len"));
            n_tx_tags = tx_tags.size();

            dout << "[MIMO OFDM RADAR] number of tx_tags:" << tx_tags.size() << std::endl;
            new_radar_frame = true;
            radar_frame_received = false;
            start_tag_offset = rx_tags[0].offset - nitems_read(d_N_tx); // offset relative to nitems_read

            if (start_tag_offset != 0)
            {
                // TODO 
            }
            
            if(tx_tags.size() > rx_tags.size())
            {
                tx_frame_offset = tx_tags.size() - rx_tags.size();

                for (int i_tx_frame = 0; i_tx_frame < tx_frame_offset; i_tx_frame++)
                {
                    n_tx_samples_discard += pmt::to_uint64(tx_tags[i_tx_frame].value);
                }
            }

            rx_packet_len = pmt::to_uint64(rx_tags[0].value);
            tx_packet_len = pmt::to_uint64(tx_tags[tx_frame_offset].value);
            dout << "[MIMO OFDM RADAR] start_tag_offset: " << start_tag_offset << ", rx_packet_len: " << rx_packet_len << ", tx_packet_len: " << tx_packet_len << std::endl;
            dout << "[MIMO OFDM RADAR] tx_frame_offset: " << tx_frame_offset << ", n_tx_samples_discard: " << n_tx_samples_discard << std::endl;
        }   

        if (new_radar_frame)
        {
            if (!radar_frame_received)
            {
                radar_frame_received = true;
            }
            else
            {
                //whole frame received --> process
            }
        }
        else
        {
            dout << "[MIMO RADAR]: No tags on RX input --> Consume all ninput_items tx0: " << ninput_items[0] <<  ", ninput_items rx0: " << ninput_items[d_N_tx] << std::endl;
            // No tags --> consume tx inputs
            for (int i_tx = 0; i_tx < d_N_tx; i_tx++)
            {
                consume(i_tx, ninput_items[i_tx]);
            }

            for (int i_rx = 0; i_rx < d_N_rx; i_rx++)
            {
                consume(i_rx+d_N_tx, ninput_items[i_rx+d_N_tx]);
            }

            return 0;
        }

        get_tags_in_range(rx_tags, d_N_tx, nitems_read(d_N_tx), nitems_read(d_N_tx) + ninput_items[d_N_tx], pmt::mp("rx_time"));
        if (rx_tags.size() != 1) 
        {
            dout << "[MIMO OFDM RADAR]: no rx_time tag found! " << std::endl;
        }
        else
        {
            dout << "[MIMO OFDM RADAR]: rx_time: " << pmt::to_double(rx_tags[0].value) << std::endl;
        }

        memset(out, 0, d_fft_len*d_N_tx*d_N_rx*d_interp_factor*sizeof(gr_complex));
        memset(radar_chan_est, 0, d_fft_len*d_N_tx*d_N_rx*sizeof(gr_complex) );       
        
        int i_chan_indx;

        gr_complex buffer_elem_mean(0.0, 0.0);

        for (int i_sc = 0; i_sc < d_fft_len; i_sc++)
        {
            for (int i_rx = 0; i_rx < d_N_rx; i_rx++)
            {
                in_rx = reinterpret_cast<const gr_complex*>(input_items[i_rx + d_N_tx]);
                in_rx += d_fft_len*d_N_pre;
                for (int i_tx = 0; i_tx < d_N_tx; i_tx++)
                {
                    in_tx = reinterpret_cast<const gr_complex*>(input_items[i_tx]);
                    in_tx += d_fft_len*d_N_pre;
                    in_tx += d_fft_len*n_tx_samples_discard;

                    if (d_enable_tx_interleave)
                    {
                        i_chan_indx = i_sc + d_fft_len * (i_tx*d_N_rx+i_rx); // col + width * (row)
                    }
                    else
                    {
                        i_chan_indx = i_sc + d_fft_len * (i_rx*d_N_tx+i_tx); // col + width * (row)
                    }

                    for (int i_sym = 0; i_sym < d_N_sym; i_sym++)
                    {
                        radar_chan_est[i_chan_indx] = radar_chan_est[i_chan_indx] + in_rx[i_sc+i_sym*d_fft_len]*std::conj(in_tx[i_sc+i_sym*d_fft_len]);
                    }

                    if (d_background_recording)
                    {
                        memcpy(&radar_chan_est_temp[i_chan_indx], radar_chan_est+i_chan_indx, sizeof(gr_complex));
                    }

                    if (d_background_removal)
                    {
                        int curr_buffer_size = radar_chan_est_buffer.size();
                        buffer_elem_mean = 0.0;

                        for (int i_buffer = 0; i_buffer < curr_buffer_size; i_buffer++)
                        {
                            buffer_elem_mean += radar_chan_est_buffer[i_buffer][i_chan_indx]/(float)curr_buffer_size;
                        }

                        radar_chan_est[i_chan_indx] = radar_chan_est[i_chan_indx] - buffer_elem_mean;
                    }
                }
            }
        }

        if (d_background_removal)
        {
            radar_chan_est_buffer.push_back(radar_chan_est_temp);
        }
        dout << "[MIMO OFDM RADAR]: curr_buffer_size = " << radar_chan_est_buffer.size() << std::endl;

        noutput_items = d_N_tx*d_N_rx;

        // Update len key tag
        pmt::pmt_t key = pmt::string_to_symbol("packet_len");
        pmt::pmt_t value = pmt::from_long(noutput_items);
        pmt::pmt_t srcid = pmt::string_to_symbol(alias());
        add_item_tag(0, nitems_written(0), key, value, srcid);    


        for(int i_pair = 0; i_pair < d_N_tx*d_N_rx; i_pair++)
        {
            memcpy(out+i_pair*d_fft_len*d_interp_factor, radar_chan_est+i_pair*d_fft_len, d_fft_len*sizeof(gr_complex));
        }
        // dout << "[PILOT PROCESSOR] pc_work_time_total: " << this->pc_work_time_total()  << std::endl;
        // dout << "[PILOT PROCESSOR] pc_input_buffers_full: " << pc_input_buffers_full(0) << "  pc_output_buffers_full: "  << pc_output_buffers_full(0) << std::endl;

        // int N;
        // std::vector<gr_complex> x(N);
        // vector<gr_complex> y(N);
        // vector<gr_complex> z(N);
        // // fill x and y with stuff
        // volk_32fc_x2_dot_prod_32fc(&z[0], &x[0], &y[0], N);

        for (int i_rx = 0; i_rx < d_N_rx; i_rx++)
        {
            consume(i_rx + d_N_tx, rx_packet_len);
        }

        for (int i_tx = 0; i_tx < d_N_tx; i_tx++)
        {
            consume(i_tx, n_tx_samples_discard + tx_packet_len);
        }

        new_radar_frame = false;

        // Tell runtime system how many output items we produced.
        return noutput_items;
    }

    void mimo_ofdm_radar_impl::set_background_record(bool background_recording)
    {
        std::cout << "[MIMO OFDM RADAR] Background recording set to  " << background_recording << std::endl;
        d_background_recording = background_recording;
    }

    void mimo_ofdm_radar_impl::capture_radar_data(bool capture_sig)
    {
        if (capture_sig)
        {
            const static Eigen::IOFormat csv_formatting(Eigen::FullPrecision, Eigen::DontAlignCols, 
                                                                    ";",   //_coeffSeparator
                                                                    ":",   //_rowSeparator
                                                                    "",     //_rowPrefix
                                                                    "",     //_rowSuffix
                                                                    "",  //_matPrefix
                                                                    ";\n");  //_matSuffix

            // radar_chan_est --> [ [fft_len], [fft_len], [fft_len], ...]
            std::ofstream file_stream(d_radar_chan_file, std::ofstream::app);
            Eigen::Map<Eigen::VectorXcf> H_radar(radar_chan_est, d_N_tx*d_N_rx*d_fft_len);

            if (file_stream.is_open())
            {
                file_stream << current_date_time2() << ", " << d_N_tx << ", " << d_N_rx << ", " << d_fft_len << ":";
                file_stream << H_radar.transpose().format(csv_formatting);
                file_stream << "\n";
                file_stream.flush();
            }
            else
            {
                throw std::runtime_error("[MIMO OFDM RADAR] Could not open file!!");
            }
            
            file_stream.close();

            std::cout << "[MIMO OFDM RADAR] Radar image captured!" << std::endl;

            // captured = true;
        }
        // else
        // {
        //     captured = false;
        // }
        
    }

} /* namespace mimo_ofdm_jrc */
} /* namespace gr */

