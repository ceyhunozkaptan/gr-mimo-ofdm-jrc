/* -*- c++ -*- */
/* 
 * Copyright 2019 gr-ofdm_radar author.
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

#ifndef INCLUDED_OFDM_RADAR_UTILS_H
#define INCLUDED_OFDM_RADAR_UTILS_H

#include <mimo_ofdm_jrc/api.h>
#include <gnuradio/block.h>
#include <gnuradio/config.h>
#include <mimo_ofdm_jrc/stream_encoder.h>
#include <cinttypes>
#include <iostream>


//FIXME
#define MAX_PAYLOAD_SIZE 3100
#define MAX_SYM (((16 + 8 * MAX_PAYLOAD_SIZE + 6) / 24) + 1) //16 zeros for scrambler + psdu + 6-bits to terminate convolutional encoder
#define MAX_ENCODED_BITS ((16 + 8 * MAX_PAYLOAD_SIZE + 6) * 2 + 288)

#define dout d_debug&& std::cout

class ofdm_mcs {
	public:
		ofdm_mcs(MCS mod_encode, int data_len);
		// Modulation and Coding Scheme
		MCS 	d_mcs;
		// rate field of the SIGNAL header
		char     rate_field; 
		// number of coded bits per sub carrier
		int      n_bpsc; 
		// number of coded bits per OFDM symbol
		int      n_cbps; 
		// number of data bits per OFDM symbol
		int      n_dbps; 

		int		d_n_data_carriers; 	

        int     n_constellations; //e.g., 3 for QPSK, 15 for 16-QAM

		void print();
};

/**
 * packet specific parameters
 */
class packet_param {
	public:
        // mcs, data_size, packet_type
		packet_param(ofdm_mcs &packet_mcs, int data_size_byte, PACKET_TYPE packet_type);
		// data size in bytes
		int data_size_byte; 
		// number of OFDM symbols (17-11)
		int n_ofdm_sym; 
		// number of padding bits in the DATA field (17-13)
		int n_pad_bits; 
		int n_encoded_bits;
		// number of data bits, including service and padding (17-12)
		int n_data_bits; 

        PACKET_TYPE packet_type;
        char packet_type_field;
        
		void print();
};

void scramble(const char *input, char *out, packet_param &frame, char initial_state);

void reset_tail_bits(char *scrambled_data, packet_param &frame);

void convolutional_encoding(const char *input, char *out, packet_param &frame);

void puncturing(const char *input, char *out, packet_param &frame, ofdm_mcs &ofdm);

void interleave(const char *input, char *out, packet_param &frame, ofdm_mcs &ofdm, bool reverse = false);

void split_symbols(const char *input, char *out, packet_param &frame, ofdm_mcs &ofdm);

void generate_bits(const char *psdu, char *data_bits, packet_param &frame, ofdm_mcs &ofdm);

int get_bit(int b, int i);

std::string current_date_time();
std::string current_date_time2();

#endif /* INCLUDED_OFDM_RADAR_UTILS_H */
