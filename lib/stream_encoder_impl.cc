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
#include "stream_encoder_impl.h"
#include <boost/crc.hpp>

namespace gr {
  namespace mimo_ofdm_jrc {

    stream_encoder::sptr
    stream_encoder::make(MCS mod_encode, int data_len, int N_ss_radar, bool debug)
    {
      return gnuradio::get_initial_sptr
        (new stream_encoder_impl(mod_encode, data_len, N_ss_radar, debug));
    }

    /*
     * The private constructor
     */
    stream_encoder_impl::stream_encoder_impl(MCS mod_encode, int data_len, int N_ss_radar, bool debug)
      : gr::block("stream_encoder",
                gr::io_signature::make(0, 0, 0),
                gr::io_signature::make(1, 1, sizeof(gr_complex))),
            d_offset(0),
			d_symbol_values(NULL),
			d_complex_symbols(NULL),
			d_debug(debug),
			d_mod_encode(mod_encode),
			d_data_len(data_len),
			d_ofdm_mcs(mod_encode,data_len),
			d_scrambler(1),
			d_N_ss_radar(N_ss_radar)
    {
        message_port_register_in(pmt::mp("pdu_in"));
		
		d_bpsk = digital::constellation_bpsk::make();
		d_qpsk = digital::constellation_qpsk::make();
		d_16qam = digital::constellation_16qam::make();

		modulator = d_bpsk;
    }

    /*
     * Our virtual destructor.
     */
    stream_encoder_impl::~stream_encoder_impl()
    {
		free(d_symbol_values);
		free(d_complex_symbols);        
    }

    int
    stream_encoder_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
        gr_complex *out = (gr_complex *)output_items[0];

		std::string str;
		const char *data_packet;
		int packet_size_byte;

        PACKET_TYPE packet_type; 

		while (!d_offset)
		{
			pmt::pmt_t msg(delete_head_nowait(pmt::intern("pdu_in")));

			if (!msg.get())
			{
				return 0;
			}
			else
			{
				gr::thread::scoped_lock lock(d_mutex);


                dout << "====================================================================================" << std::endl;
				dout << "[STREAM ENCODER] New PDU Received -> Encoding STARTS..." << std::endl;

				if(pmt::is_symbol(msg)) 
				{
					str = pmt::symbol_to_string(msg);
					packet_size_byte = str.length();
					data_packet = str.data();

                    packet_type = (PACKET_TYPE) str.at(0);
                    dout << "[STREAM ENCODER] PDU is a symbol " << std::endl;

				} 
				else if(pmt::is_pair(msg)) 
				{
					packet_size_byte = pmt::blob_length(pmt::cdr(msg));
					data_packet = reinterpret_cast<const char *>(pmt::blob_data(pmt::cdr(msg)));
                    packet_type = (PACKET_TYPE) data_packet[0];
                    
                    dout << "[STREAM ENCODER] PDU is a pair " << std::endl;
				} 
				else 
				{
					throw std::invalid_argument("[STREAM ENCODER] Encoder expects PDUs or strings");
					return 0;
				}


				d_ofdm_mcs = ofdm_mcs(d_mod_encode, d_data_len);
				packet_param burst(d_ofdm_mcs, packet_size_byte+4, packet_type); //+4 added for 32bit CRC checksum

                // int max_ofdm_sym = (((16 + 8 * MAX_PAYLOAD_SIZE + 6) / ((double) d_ofdm_mcs.n_dbps)) + 1); // 16 zeros for scrambler + psdu + 6-bits to terminate convolutional encoder
                // std::cout << "[STREAM ENCODER] PDU too Large -> Maximun PDU Length (byte): " << max_ofdm_sym << std::endl;

                //FIXME make sure this check is correct
				if(packet_size_byte + 4 > MAX_PAYLOAD_SIZE) // 4 bytes for CRC32, 3 bytes for SERVICE and 1 bytes for Padding Bits (>6 bits)
				{
					std::cout << "[STREAM ENCODER] Data Packet too Large -> Maximun Packet Size (byte): " << MAX_PAYLOAD_SIZE << std::endl;
					return 0;
				}

                dout << "[STREAM ENCODER] Packet Type: " << packet_type << std::endl;
				dout << "[STREAM ENCODER] Data Packet Size with CRC (byte): " << (packet_size_byte + 4) << std::endl;
				dout << "[STREAM ENCODER] OFDM Symbols: " << burst.n_ofdm_sym << ", Bits per OFDM Symbol: " << d_ofdm_mcs.n_dbps << ", Total Data Bits: " << burst.n_data_bits << std::endl;
				dout << "[STREAM ENCODER] Padding Bits: " << burst.n_pad_bits << ", Encoded Bits: " << burst.n_encoded_bits <<  std::endl;

                // 4 bytes (32-bit) CRC (FCS) at the end
				char *data_packet_crc =(char*)calloc(packet_size_byte + 4, sizeof(char));
				//compute and store fcs
				boost::crc_32_type result;
				result.process_bytes(data_packet, packet_size_byte);

				uint32_t fcs = result.checksum();
				memcpy(data_packet_crc, data_packet, packet_size_byte*sizeof(char));
				memcpy(data_packet_crc + packet_size_byte, &fcs, sizeof(uint32_t));
				
                dout << "[STREAM ENCODER] Encoding following bytes with CRC: " << burst.n_pad_bits << ", Encoded Bits: " << burst.n_encoded_bits <<  std::endl;
                print_output(data_packet_crc, packet_size_byte+4);

				char *data_bits 		= (char*)calloc(burst.n_data_bits, sizeof(char));
				char *scrambled_data   	= (char*)calloc(burst.n_data_bits, sizeof(char));
				char *encoded_data     	= (char*)calloc(burst.n_data_bits * 2, sizeof(char));
				char *punctured_data   	= (char*)calloc(burst.n_encoded_bits, sizeof(char));
				char *interleaved_data 	= (char*)calloc(burst.n_encoded_bits, sizeof(char));
				char *symbols          = (char*)calloc((burst.n_encoded_bits / d_ofdm_mcs.n_bpsc), sizeof(char));

				//generate the data field, adding service field and pad bits
				generate_bits(data_packet_crc, data_bits, burst, d_ofdm_mcs);

				// scrambling
				scramble(data_bits, scrambled_data, burst, d_scrambler++);
				if(d_scrambler > 127) 
				{
					d_scrambler = 1;
				}

				// reset tail bits
				reset_tail_bits(scrambled_data, burst);
				// encoding
				convolutional_encoding(scrambled_data, encoded_data, burst);
				// puncturing
				puncturing(encoded_data, punctured_data, burst, d_ofdm_mcs);
				// no interleaving applied
				// interleave(punctured_data, interleaved_data, burst, d_ofdm_mcs);

				// one byte per symbol
				split_symbols(punctured_data, symbols, burst, d_ofdm_mcs);

				d_symbol_len = burst.n_ofdm_sym * d_data_len;

				d_symbol_values = (char *)calloc((d_symbol_len), sizeof(char));
				std::memcpy(d_symbol_values, symbols, d_symbol_len);

				switch (d_ofdm_mcs.d_mcs) {
					case BPSK_1_2:
					case BPSK_3_4:
						modulator = d_bpsk;
						break;

					case QPSK_1_2:
					case QPSK_3_4:
						modulator = d_qpsk;
						break;

					case QAM16_1_2:
					case QAM16_3_4:
						modulator = d_16qam;
						break;

					default:
						throw std::invalid_argument("wrong encoding");
						break;
				}
				d_complex_symbols = (gr_complex *)calloc((d_symbol_len), sizeof(gr_complex));
				for (int i = 0; i < d_symbol_len; i++)
				{
					modulator->map_to_points(d_symbol_values[i], &d_complex_symbols[i]);

                    if(modulator->bits_per_symbol() == d_qpsk->bits_per_symbol()){
                        d_complex_symbols[i] = d_complex_symbols[i] / (float)2.0;
                    }
				}

                // add tags
                pmt::pmt_t key = pmt::string_to_symbol("packet_len");
                pmt::pmt_t value = pmt::from_long(d_symbol_len);
                pmt::pmt_t srcid = pmt::string_to_symbol(alias());
                add_item_tag(0, nitems_written(0), key, value, srcid);
                
                pmt::pmt_t type_tag = pmt::from_long(packet_type);
                add_item_tag(0, nitems_written(0), pmt::mp("packet_type"), type_tag, srcid);

                pmt::pmt_t mcs_tag = pmt::from_long(d_ofdm_mcs.d_mcs);
                add_item_tag(0, nitems_written(0), pmt::mp("mcs"), mcs_tag, srcid);

                int total_pdu_length = packet_size_byte+4;
                pmt::pmt_t pdu_bytes = pmt::from_long(total_pdu_length); //+4 added for 32bit CRC checksum
                add_item_tag(0, nitems_written(0), pmt::mp("pdu_len"), pdu_bytes, srcid);


				free(data_packet_crc);
				free(data_bits);
				free(scrambled_data);
				free(encoded_data);
				free(punctured_data);
				free(interleaved_data);
				free(symbols);
                
				break;
			}
		} 

		int n_out = std::min(noutput_items, d_symbol_len - d_offset);


		std::memcpy(out, d_complex_symbols + d_offset, n_out * sizeof(gr_complex));
		d_offset += n_out;

		if (d_offset == d_symbol_len)
		{
			d_offset = 0;
			free(d_symbol_values);
			free(d_complex_symbols);
			d_symbol_values = 0;
			d_complex_symbols = 0;
		}

		dout << "[STREAM ENCODER] Encoding COMPLETED --> Total Data Symbols: " << n_out << std::endl;
		dout << "========================================================================================" << std::endl;
		return n_out;
    }
    void stream_encoder_impl::set_mcs(MCS mod_encode) 
    {
		std::cout << "[STREAM ENCODER] MCS Changed to  " << mod_encode << std::endl;
		gr::thread::scoped_lock lock(d_mutex);
		d_mod_encode = mod_encode;
		d_ofdm_mcs = ofdm_mcs(mod_encode, d_data_len);
	}

	void stream_encoder_impl::print_output(const char *data_packet, int packet_size_byte) 
	{
		dout << std::setfill('0') << std::setw(4) << 0 << ": ";;
		for(int i = 0; i < packet_size_byte; i++) {
			dout << std::setfill('0') << std::setw(2) << std::hex << (unsigned int) (unsigned char) data_packet[i] << " " ;
			if(i % 16 == 15) {
				dout << std::endl;
				dout << std::setfill('0') << std::setw(4) << i/15 << ": ";
			}
		}
		dout << std::dec << std::endl;
	}
  } /* namespace mimo_ofdm_jrc */
} /* namespace gr */

