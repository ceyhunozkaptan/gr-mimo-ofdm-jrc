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

#include "utils.h"
#include <cassert>
#include <cstring>
#include <math.h>

packet_param::packet_param(ofdm_mcs &ofdm, int data_size_byte, PACKET_TYPE packet_type) 
{
	this->data_size_byte = data_size_byte;
    
    //16 zeros for scrambler + data packet bytes + 6 (at least) bits to terminate convolutional encoder
	n_ofdm_sym = (int) ceil((16 + 8 * data_size_byte + 6) / (double) ofdm.n_dbps); 
    
	n_data_bits = n_ofdm_sym * ofdm.n_dbps;

	// number of padding bits (17-13)
	n_pad_bits = n_data_bits - (16 + 8 * data_size_byte + 6); //16 bit header (16 zeros) [, 16 info moved to SIG]

	n_encoded_bits = n_ofdm_sym * ofdm.n_cbps;
    
    this->packet_type = packet_type;

    switch(packet_type) {
		case NDP:
            packet_type_field = 0b00000000;
            break;
        case DATA:
            packet_type_field = 0b00000001;
            break;
        defaut:
			assert(false);
			break;
    }
}

ofdm_mcs::ofdm_mcs(MCS mod_encode, int n_data_carriers) {
	d_mcs = mod_encode;
	d_n_data_carriers = n_data_carriers;
	switch(d_mcs) {
		case BPSK_1_2:
			n_bpsc = 1;
			n_cbps = d_n_data_carriers*n_bpsc;
			n_dbps = n_cbps/2;
			rate_field = 0x0D; // 0b00001101
            n_constellations = 2;
			break;

		case BPSK_3_4:
			n_bpsc = 1;
			n_cbps = d_n_data_carriers*n_bpsc;
			n_dbps = n_cbps*3/4;
			rate_field = 0x0F; // 0b00001111
            n_constellations = 2;
			break;
		
		case QPSK_1_2:
			n_bpsc = 2;
			n_cbps = d_n_data_carriers*n_bpsc;
			n_dbps = n_cbps/2;
			rate_field = 0x05; // 0b00000101
            n_constellations = 4;
			break;

		case QPSK_3_4:
			n_bpsc = 2;
			n_cbps = d_n_data_carriers*n_bpsc;
			n_dbps = n_cbps*3/4;
			rate_field = 0x07; // 0b00000111
            n_constellations = 4;
			break;

		case QAM16_1_2:
			n_bpsc = 4;
			n_cbps = d_n_data_carriers*n_bpsc;
			n_dbps = n_cbps/2;
			rate_field = 0x09; // 0b00001001
            n_constellations = 16;
			break;

		case QAM16_3_4:
			n_bpsc = 4;
			n_cbps = d_n_data_carriers*n_bpsc;
			n_dbps = n_cbps*3/4;
			rate_field = 0x0B; // 0b00001011
            n_constellations = 16;
			break;

		defaut:
			assert(false);
			break;
	}
}

void
ofdm_mcs::print() {
	std::cout << "OFDM Parameters:" << std::endl;
	std::cout << "endcoding :" << d_mcs << std::endl;
	std::cout << "rate_field :" << (int)rate_field << std::endl;
	std::cout << "n_bpsc :" << n_bpsc << std::endl;
	std::cout << "n_cbps :" << n_cbps << std::endl;
	std::cout << "n_dbps :" << n_dbps << std::endl;
}

void
packet_param::print() {
	std::cout << "FRAME Parameters:" << std::endl;
	std::cout << "data_size_byte: " << data_size_byte << std::endl;
	std::cout << "n_ofdm_sym: " << n_ofdm_sym << std::endl;
	std::cout << "n_pad_bits: " << n_pad_bits << std::endl;
	std::cout << "n_encoded_bits: " << n_encoded_bits << std::endl;
	std::cout << "n_data_bits: " << n_data_bits << std::endl;
}

int get_bit(int b, int i) {
	return (b & (1 << i) ? 1 : 0);
}

void generate_bits(const char *psdu, char *data_bits, packet_param &frame, ofdm_mcs &ofdm) 
{
    //TODO IF CHANGED, n_ofdm_sym and crc checks should be updated

	// first 16 bits are zero to reset scrambler (SERVICE field)
	memset(data_bits, 0, 16);
	data_bits += 16;
        
    // // next 16 bits are reserved for packet info 
	// data_bits[ 0] = get_bit(ofdm.rate_field, 0);
	// data_bits[ 1] = get_bit(ofdm.rate_field, 1);
	// data_bits[ 2] = get_bit(ofdm.rate_field, 2);
	// data_bits[ 3] = get_bit(ofdm.rate_field, 3);
	// data_bits[ 4] = get_bit(length,  0);
	// data_bits[ 5] = get_bit(length,  1);
	// data_bits[ 6] = get_bit(length,  2);
	// data_bits[ 7] = get_bit(length,  3);
	// data_bits[ 8] = get_bit(length,  4);
	// data_bits[ 9] = get_bit(length,  5);
	// data_bits[10] = get_bit(length,  6);
	// data_bits[11] = get_bit(length,  7);
	// data_bits[12] = get_bit(length,  8);
	// data_bits[13] = get_bit(length,  9);
	// data_bits[14] = get_bit(length,  10);
	// data_bits[15] = get_bit(length,  11);
	// data_bits += 16;

    //TODO: CRC check at the decoder should ignore first 4 bytes 
    // CRC is computed with psdu before these field

	for(int i = 0; i < frame.data_size_byte; i++) {
		for(int b = 0; b < 8; b++) {
			data_bits[i * 8 + b] = !!(psdu[i] & (1 << b));
			// std::cout << "[generate_bits]" << i * 8 + b << "th data_bit: " << (int) data_bits[i * 8 + b] << std::endl;
		}
	}
}

void scramble(const char *in, char *out, packet_param &frame, char initial_state) 
{
    int state = initial_state;
    int feedback;

    for (int i = 0; i < frame.n_data_bits; i++) 
	{
		feedback = (!!(state & 64)) ^ (!!(state & 8));
		out[i] = feedback ^ in[i];
		state = ((state << 1) & 0x7e) | feedback;
    }
}


void reset_tail_bits(char *scrambled_data, packet_param &frame) 
{
	memset(scrambled_data + frame.n_data_bits - frame.n_pad_bits - 6, 0, 6 * sizeof(char));
}


int ones(int n) 
{
	int sum = 0;
	for(int i = 0; i < 8; i++) {
		if(n & (1 << i)) {
			sum++;
		}
	}
	return sum;
}


void convolutional_encoding(const char *in, char *out, packet_param &frame) 
{
	int state = 0;

	for(int i = 0; i < frame.n_data_bits; i++) {
		assert(in[i] == 0 || in[i] == 1);
		state = ((state << 1) & 0x7e) | in[i];
		out[i * 2]     = ones(state & 0155) % 2;
		out[i * 2 + 1] = ones(state & 0117) % 2;
	}
}


void puncturing(const char *in, char *out, packet_param &frame, ofdm_mcs &ofdm) 
{
	int mod;

	for (int i = 0; i < frame.n_data_bits * 2; i++) {
		switch(ofdm.d_mcs) {
			case BPSK_1_2:
			case QPSK_1_2:
			case QAM16_1_2:
				*out = in[i];
				out++;
				break;

			case BPSK_3_4:
			case QPSK_3_4:
			case QAM16_3_4:
			// case QAM64_3_4:
				mod = i % 6;
				if (!(mod == 3 || mod == 4)) {
					*out = in[i];
					out++;
				}
				break;
			defaut:
				assert(false);
				break;
		}
	}
}


void interleave(const char *in, char *out, packet_param &frame, ofdm_mcs &ofdm, bool reverse) 
{
	int n_cbps = ofdm.n_cbps;
	int first[n_cbps];
	int second[n_cbps];
	int s = std::max(ofdm.n_bpsc / 2, 1);

	for(int j = 0; j < n_cbps; j++) {
		first[j] = s * (j / s) + ((j + int(floor(16.0 * j / n_cbps))) % s);
	}

	for(int i = 0; i < n_cbps; i++) {
		second[i] = 16 * i - (n_cbps - 1) * int(floor(16.0 * i / n_cbps));
	}

	for(int i = 0; i < frame.n_ofdm_sym; i++) {
		for(int k = 0; k < n_cbps; k++) {
			if(reverse) {
				out[i * n_cbps + second[first[k]]] = in[i * n_cbps + k];
			} else {
				out[i * n_cbps + k] = in[i * n_cbps + second[first[k]]];
			}
		}
	}
}


void split_symbols(const char *in, char *out, packet_param &frame, ofdm_mcs &ofdm) 
{
	int symbols = frame.n_ofdm_sym * ofdm.d_n_data_carriers;

	for (int i = 0; i < symbols; i++) 
	{
		out[i] = 0;
		for(int k = 0; k < ofdm.n_bpsc; k++) 
		{
			// std::cout << "[split_symbols] i: " << i << "  k:" << k << " in:" << (uint) *in  << std::endl;
			assert(*in == 1 || *in == 0);
			out[i] |= (*in << k);
			in++;
		}
	}
}

std::string current_date_time()
{
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);

    strftime(buf, sizeof(buf), "%m-%d-%Y %H:%M:%S", &tstruct);

    return buf;
}

std::string current_date_time2()
{
	auto now = std::chrono::system_clock::now();
	auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;
	auto timer = std::chrono::system_clock::to_time_t(now);
	std::tm bt = *std::localtime(&timer);
	std::ostringstream oss;

	oss << std::put_time(&bt, "%H:%M:%S"); // HH:MM:SS
    oss << '.' << std::setfill('0') << std::setw(3) << ms.count();

    return oss.str();
}


