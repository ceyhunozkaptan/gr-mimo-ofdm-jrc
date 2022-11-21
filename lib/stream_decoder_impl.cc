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
#include "stream_decoder_impl.h"
// #include <chrono>

namespace gr {
  namespace mimo_ofdm_jrc {

    stream_decoder::sptr
    stream_decoder::make(int n_data_carriers,
                            const std::string& comm_log_file,
                            bool stats_record,  
                            bool debug)
    {
      return gnuradio::get_initial_sptr
        (new stream_decoder_impl(n_data_carriers,                            
                                    comm_log_file,
                                    stats_record,  
                                    debug));
    }


    /*
     * The private constructor
     */
    stream_decoder_impl::stream_decoder_impl(int n_data_carriers, 
                                                const std::string& comm_log_file,
                                                bool stats_record, 
                                                bool debug)
      : gr::block("stream_decoder",
              gr::io_signature::make(1, 1, n_data_carriers*sizeof(gr_complex)),
              gr::io_signature::make(0, 1, sizeof(float))),
            d_debug(debug),
            d_n_data_carriers(n_data_carriers),
            d_ofdm_mcs(BPSK_1_2, n_data_carriers),
            d_stream_param(d_ofdm_mcs, 0, PACKET_TYPE::NDP), 
            d_comm_log_file(comm_log_file),
            d_stats_record(stats_record),
            d_frame_rx_complete(true),
            per_stats(bt::rolling_window::window_size = 25),
            snr_data_stats(bt::rolling_window::window_size = 1)
    {
            message_port_register_out(pmt::mp("sym"));
            message_port_register_out(pmt::mp("stats"));

            dout << "[STREAM DECODER] Starting... d_n_data_carriers: " << d_n_data_carriers << std::endl;
            d_rx_symbols = new uint8_t[MAX_SYM*d_n_data_carriers];
            
            // d_rx_bits = new uint8_t[MAX_ENCODED_BITS];
            d_rx_bits = (uint8_t*) calloc((MAX_ENCODED_BITS), sizeof(uint8_t));

            d_bpsk = digital::constellation_bpsk::make(); 
			d_qpsk = digital::constellation_qpsk::make();
			d_16qam = digital::constellation_16qam::make();

            perf_display_interval = 1.0;
            last_perf_time = std::chrono::system_clock::now();
            d_new_stat_started = false;

            set_tag_propagation_policy(TPP_DONT);

            dout << "[STREAM DECODER] Log file: " << d_comm_log_file << std::endl;
            std::ifstream file_stream(d_comm_log_file);
            if (!file_stream.is_open())
            {
                std::cerr << "[MIMO PRECODER] Could not open log file!" << std::endl;
            }
            file_stream.close();
    }

    /*
     * Our virtual destructor.
     */
    stream_decoder_impl::~stream_decoder_impl()
    {
        delete[] d_rx_symbols;
        free(d_rx_bits);
    }

    int
    stream_decoder_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
        int n_input_items = ninput_items[0];

        const gr_complex *in = (const gr_complex*) input_items[0];
        // dout << "[STREAM DECODER] n_input_items: " << n_input_items << "  noutput_items " << noutput_items << ", output_items.size(): " << output_items.size() <<  std::endl;
        
        float* per_out = (float*) output_items[0];

        std::vector<float> per_vector;

        std::vector<gr::tag_t> tags;
        const uint64_t nread = this->nitems_read(0);

        gr::block::set_thread_priority(60);
        gr::thread::scoped_lock lock(d_setlock);

        int n_in = 0;
        int n_out = 0;

        while(n_in < n_input_items) 
        {
            get_tags_in_range(tags, 0, nread + n_in, nread + n_in + 1, pmt::string_to_symbol("stream_start"));

            if(tags.size()) 
            {
                dout << "================================================================================" << std::endl;

                if (d_frame_rx_complete == false) 
                {
                    dout << "[STREAM DECODER] Warning: starting to receive new frame before old frame was complete" << std::endl;
                    dout << "[STREAM DECODER] Already n_copied " << n_copied << " out of " << d_stream_param.n_ofdm_sym << " symbols of last frame" << std::endl;
                }
                d_frame_rx_complete = false;

                pmt::pmt_t dict = tags[0].value;
                d_snr_est = pmt::to_double(pmt::dict_ref(dict, pmt::mp("snr"), pmt::from_double(0)));
                dout << "[STREAM DECODER] SNR estimate from equalizer: " << d_snr_est << std::endl;
                d_data_length = pmt::to_uint64(pmt::dict_ref(dict, pmt::mp("data_bytes"), pmt::from_double(0)));
                dout << "[STREAM DECODER] DATA length from equalizer: " << d_data_length << std::endl;
                d_mcs = (MCS) pmt::to_uint64(pmt::dict_ref(dict, pmt::mp("mcs"), pmt::from_double(0)));
                dout << "[STREAM DECODER] MCS from equalizer: " << d_mcs << std::endl;
                d_packet_type = (PACKET_TYPE) pmt::to_uint64(pmt::dict_ref(dict, pmt::mp("packet_type"), pmt::from_double(0)));
                dout << "[STREAM DECODER] PacketType from equalizer: " << d_packet_type << std::endl;

                // check for maximum frame size
                ofdm_mcs ofdm = ofdm_mcs(d_mcs, d_n_data_carriers);
                packet_param stream = packet_param(ofdm, d_data_length, d_packet_type); 

                if (stream.n_ofdm_sym <= MAX_SYM && stream.data_size_byte <= MAX_PAYLOAD_SIZE) 
                {
                    d_ofdm_mcs = ofdm;
                    d_stream_param = stream;
                    n_copied = 0;
                    d_start_decoding = true;
                    setup_demodulator(d_mcs);
                    dout << "[STREAM DECODER] Decoding starts --> " << d_stream_param.n_ofdm_sym  << "  symbols,  " << d_stream_param.data_size_byte << "  bytes" << std::endl;
                    dout << "[STREAM DECODER] d_demodulator:" << d_demodulator->bits_per_symbol() << std::endl;

                    // dout << "[STREAM DECODER] Bits per Subcarrier:" << d_ofdm_mcs.n_bpsc << ", Total Data Bits:" << d_stream_param.n_data_bits << ", Total OFDM Symbols:" << d_stream_param.n_ofdm_sym << std::endl;
                } else 
                {
                    dout << "[STREAM DECODER] Dropping packet is too large!!" << std::endl;
                    d_start_decoding = false;
                }
            }

            get_tags_in_range(tags, 0, nread + n_in, nread + n_in + 1, pmt::mp("stream_end"));
            if (tags.size()) 
            {
                pmt::pmt_t dict = tags[0].value;
                d_snr_data_est = pmt::to_float(pmt::dict_ref(dict, pmt::mp("snr_data"), pmt::from_double(0)));
                chan_est_mean = pmt::c32vector_elements( (pmt::dict_ref(dict, pmt::mp("chan_mean"), pmt::from_double(0))));

                // d_snr_data_est = pmt::to_float(tags[0].value);
                dout << "[STREAM DECODER] Data SNR: " << d_snr_data_est <<" ,chan_est_mean.size: " << chan_est_mean.size() << std::endl;

            }

            if( n_copied < d_stream_param.n_ofdm_sym && d_start_decoding ) 
            {
                for (int i_data = 0; i_data < d_n_data_carriers; i_data++)
                {
                    d_rx_symbols[n_copied * d_n_data_carriers + i_data] = d_demodulator->decision_maker(&in[i_data]);
                }
                n_copied++;

                if(n_copied == d_stream_param.n_ofdm_sym) 
                {
                    dout << "[STREAM DECODER] Complete Frame Received -> Decoding STARTS..." << std::endl;
                    dout << "[STREAM DECODER] d_comm_log_file:" << d_comm_log_file << std::endl;
                    decode();
                    dout << "[STREAM DECODER] Decode done!" << std::endl;
                    in += d_n_data_carriers;
                    n_in++;
                    d_frame_rx_complete = true;
                    n_out++;
                    per_vector.push_back(100.0*boost::accumulators::rolling_mean(per_stats));
                    n_copied = 0;
                    break;
                }
            }

            in += d_n_data_carriers;
            n_in++;
        }

        // std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
        // std::chrono::duration<float> diff = now - last_perf_time;

        // if (diff.count() >= perf_display_interval)
        // {
        //     std::cout << "===========================================================" << std::endl;
        //     std::cout << "[STREAM DECODER] Work Time: " << pc_work_time_avg() <<  ", Throughput: " <<  pc_throughput_avg() << std::endl;
        //     std::cout << "[STREAM DECODER] Input Buffer: " << pc_input_buffers_full(0) << ",  Output Buffer: "  << pc_output_buffers_full(0) << std::endl;
        //     last_perf_time = std::chrono::system_clock::now();
        // }

        memcpy(per_out, &per_vector[0], n_out * sizeof(float));

        consume(0, n_in);
        return n_out;
    }

    void stream_decoder_impl::decode() 
    {

            // dout << "[STREAM DECODER] d_comm_log_file:" << d_comm_log_file << std::endl;
            std::ofstream file_stream(d_comm_log_file, std::ofstream::app);

            if(d_stats_record)
            {
                dout << "[STREAM DECODER] d_comm_log_file:" << d_comm_log_file << ", " << file_stream.is_open() << std::endl;

                if (file_stream.is_open())
                {
                    if(!d_new_stat_started)
                    {
                        file_stream << "\n NEW RECORD - " << current_date_time() << "\n";;
                        file_stream.flush();
                        d_new_stat_started = true;
                    }
                }
                else
                {
                    dout << "[STREAM DECODER] d_comm_log_file:" << d_comm_log_file << ", " << file_stream.is_open() << std::endl;

                    throw std::runtime_error("[STREAM DECODER] Could not open file!!");
                }
            }

           
            for(int i = 0; i < d_stream_param.n_ofdm_sym * d_n_data_carriers; i++) 
            {
                for(int k = 0; k < d_ofdm_mcs.n_bpsc; k++) 
                {
                    d_rx_bits[i*d_ofdm_mcs.n_bpsc + k] = !!(d_rx_symbols[i] & (1 << k));
                }
            }

            // deinterleave();

            // TODO change decoder return type 
            uint8_t* decoded = d_viterbi_decoder.decode(&d_ofdm_mcs, &d_stream_param, d_rx_bits);
            descramble(decoded);

            // Printing decoded data bytes
            dout << "[STREAM DECODER] Decoded bytes with CRC:" << d_stream_param.data_size_byte << std::endl;
            // print_output(out_bytes+2, d_stream_param.data_size_byte); //2(SERVICE) + PDU + 4(CRC32) 

            //CHECKSUM
            boost::crc_32_type result;
            result.process_bytes(out_bytes+2, d_stream_param.data_size_byte); // +2 for 16-bit zeros prepended to reset the scrambler
            if(result.checksum() != 558161692) 
            {
                std::cerr << "[STREAM DECODER] Data Checksum is WRONG!!! --> Dropping Packet, bytes:" << d_stream_param.data_size_byte << std::endl;

                // create PDU
                send_out_bytes[0] = (uint8_t) 0; 
                send_out_bytes[1] = d_packet_type;
                memcpy(send_out_bytes+2, &d_snr_est, sizeof(float));
                memcpy(send_out_bytes+2+sizeof(float), &d_snr_data_est, sizeof(float));
                // copy data bytes
                memcpy(send_out_bytes+info_bytes, out_bytes+2, d_stream_param.data_size_byte-4);
                // create PDU
                pmt::pmt_t blob = pmt::make_blob(send_out_bytes, info_bytes+d_stream_param.data_size_byte-4); // -4 for CRC
                pmt::pmt_t dict = pmt::make_dict();
                dict = pmt::dict_add(dict, pmt::mp("SNR"), pmt::from_double(d_snr_est));
                message_port_pub(pmt::mp("sym"), pmt::cons(dict, blob));

                float per_val = 100.0*boost::accumulators::rolling_mean(per_stats);
                pmt::pmt_t d_per_key = pmt::string_to_symbol("per");          // identifier 
                pmt::pmt_t d_per_value = pmt::init_f32vector(1,&per_val); // pmt
                pmt::pmt_t d_per_pack = pmt::list2(d_per_key, d_per_value); // make list 

                float snr_val = boost::accumulators::rolling_mean(snr_data_stats);
                pmt::pmt_t d_snr_key = pmt::string_to_symbol("snr");          // identifier 
                pmt::pmt_t d_snr_value = pmt::init_f32vector(1, &snr_val); // pmt
                pmt::pmt_t d_snr_pack = pmt::list2(d_snr_key, d_snr_value); // make list 

                pmt::pmt_t stats_msg = pmt::list2(d_per_pack, d_snr_pack);
                message_port_pub(pmt::mp("stats"), stats_msg); // publish message

                per_stats(1);
                snr_data_stats(d_snr_data_est);

                if(d_stats_record)
                {       
                    dout << "[STREAM DECODER] d_comm_log_file:" << d_comm_log_file << ", " << file_stream.is_open() << std::endl;

                    if (file_stream.is_open())
                    {
                        file_stream << current_date_time2() << ", \t" << 0 << ", \t" << (int) d_packet_type << ", \t" << (int) d_mcs << ", \t" << d_snr_est << ", \t" << d_snr_data_est << ", \t" << d_stream_param.data_size_byte << ", \t";
                        for (int i = 0; i < chan_est_mean.size(); i++)
                        {
                            file_stream << chan_est_mean[i] << ";";
                        }
                        file_stream << "\n";
                        file_stream.flush();
                    }
                    else
                    {
                        dout << "[STREAM DECODER] d_comm_log_file:" << d_comm_log_file << ", " << file_stream.is_open() << std::endl;
                        throw std::runtime_error("[OFDM Equalizer] Could not open file!!");
                    }
                    // file_stream.close();
                }

                return;
            }

            dout << "[STREAM DECODER] Checksum CORRECT!" << std::endl;
            total_data_received += d_stream_param.data_size_byte;
            per_stats(0);
            snr_data_stats(d_snr_data_est);
            dout << "[STREAM DECODER] Estimated SNR:" << d_snr_est << std::endl;
            
            // time (&curr_time);

            // curr_time = time_now;
            // double diff_ms = std::chrono::duration_cast<std::chrono::microseconds>(curr_time - last_update_time).count()/1e3;

            // if (diff_ms >= perf_display_interval*1e3)
            // {
            //     std::cout << "[STREAM DECODER] Average Data Throughput:" << total_data_received*8*1e3/diff_ms << std::endl;
            //     total_data_received = 0;
            //     last_update_time = time_now;
            // }

            //copy info bytes
            send_out_bytes[0] = (uint8_t) 1; 
            send_out_bytes[1] = d_packet_type;
            memcpy(send_out_bytes+2, &d_snr_est, sizeof(float));
            memcpy(send_out_bytes+2+sizeof(float), &d_snr_data_est, sizeof(float));
            // copy data bytes
            memcpy(send_out_bytes+info_bytes, out_bytes+2, d_stream_param.data_size_byte-4);
            // create PDU
            pmt::pmt_t blob = pmt::make_blob(send_out_bytes, info_bytes+d_stream_param.data_size_byte-4); // -4 for CRC
            // pmt::pmt_t enc = pmt::from_uint64(d_ofdm_mcs.d_mcs);
            pmt::pmt_t dict = pmt::make_dict();
            dict = pmt::dict_add(dict, pmt::mp("SNR"), pmt::from_double(d_snr_est));
            message_port_pub(pmt::mp("sym"), pmt::cons(dict, blob));

            pmt::pmt_t d_per_key = pmt::string_to_symbol("per");          // identifier 
            float per_val = 100.0*boost::accumulators::rolling_mean(per_stats);
            pmt::pmt_t d_per_value = pmt::init_f32vector(1,&per_val); // pmt
            pmt::pmt_t d_per_pack = pmt::list2(d_per_key, d_per_value); // make list 


            float snr_val = boost::accumulators::rolling_mean(snr_data_stats);
            pmt::pmt_t d_snr_key = pmt::string_to_symbol("snr");          // identifier 
            pmt::pmt_t d_snr_value = pmt::init_f32vector(1, &snr_val); // pmt
            pmt::pmt_t d_snr_pack = pmt::list2(d_snr_key, d_snr_value); // make list 

            pmt::pmt_t stats_msg = pmt::list2(d_per_pack, d_snr_pack);
            message_port_pub(pmt::mp("stats"), stats_msg); // publish message

            if(d_stats_record)
            {
                dout << "[STREAM DECODER] d_comm_log_file:" << d_comm_log_file << ", " << file_stream.is_open() << std::endl;
                if (file_stream.is_open())
                {
                    file_stream << current_date_time2() << ", \t" << 1 << ", \t" << (int) d_packet_type << ", \t" << (int) d_mcs << ", \t" << d_snr_est << ", \t" << d_snr_data_est << ", \t" << d_stream_param.data_size_byte << ", \t" ;
                    for (int i = 0; i < chan_est_mean.size(); i++)
                    {
                        file_stream << chan_est_mean[i] << ";";
                    }
                    file_stream << "\n";
                    file_stream.flush();
                }
                else
                {
                    dout << "[STREAM DECODER] d_comm_log_file:" << d_comm_log_file << ", " << file_stream.is_open() << std::endl;
                    throw std::runtime_error("[STREAM DECODER] Could not open file!!");
                }
                // file_stream.close();
            }

            file_stream.close();
            
        }

        void stream_decoder_impl::descramble (uint8_t *decoded_bits) 
        {
            int state = 0;
            std::memset(out_bytes, 0, d_stream_param.data_size_byte+2);

            for(int i = 0; i < 7; i++) 
            {
                if(decoded_bits[i]) 
                {
                    state |= 1 << (6 - i);
                }
            }
            out_bytes[0] = state;

            int feedback;
            int bit;

            for(int i = 7; i < d_stream_param.data_size_byte*8+16; i++) 
            {
                feedback = ((!!(state & 64))) ^ (!!(state & 8));
                bit = feedback ^ (decoded_bits[i] & 0x1);
                out_bytes[i/8] |= bit << (i%8);
                state = ((state << 1) & 0x7e) | feedback;
            }
        }

        void stream_decoder_impl::print_output(const uint8_t *psdu, int psdu_length) 
        {
            dout << "[STREAM DECODER] PSDU Length: " << psdu_length << std::endl;
            dout << std::setfill('0') << std::setw(4) << 0 << ": ";;
            for(int i = 0; i < psdu_length; i++) 
            {
                dout << std::setfill('0') << std::setw(2) << std::hex << (unsigned int) psdu[i] << std::dec << " " ;
                if(i % 16 == 15) 
                {
                    dout << std::endl;
                    dout << std::setfill('0') << std::setw(4) << i/15 << ": ";
                }
            }
            dout << std::dec << std::endl;
        }
        
        void stream_decoder_impl::setup_demodulator(MCS mcs_value)
        {
            switch (mcs_value) {
            case BPSK_1_2:
            case BPSK_3_4:
                d_demodulator = d_bpsk;
                break;
            case QPSK_1_2:
            case QPSK_3_4:
                d_demodulator = d_qpsk;
                break;
            case QAM16_1_2:
            case QAM16_3_4:
                d_demodulator = d_16qam;
                break;
            default:
                std::cout << "[STREAM DECODER] Something is wrong! --> mcs_value is not correct: " << (int) mcs_value << std::endl;
            }
        }

        void stream_decoder_impl::set_stats_record(bool stats_record)
        {
            gr::thread::scoped_lock lock(d_mutex);
            d_stats_record = stats_record;
            d_new_stat_started = false;
            
            if(stats_record)
            {
                dout << "[STREAM DECODER] Recording Stats --> Enabled" << std::endl;
            }
            else{
                dout << "[STREAM DECODER] Recording Stats --> Disabled" << std::endl;
            }

        }
  } /* namespace mimo_ofdm_jrc */
} /* namespace gr */

