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
#include "matrix_transpose_impl.h"

namespace gr {
  namespace mimo_ofdm_jrc {

    matrix_transpose::sptr
    matrix_transpose::make(int input_len, int output_len, int interp_factor,bool debug, std::string len_key)
    {
      return gnuradio::get_initial_sptr
        (new matrix_transpose_impl(input_len, output_len, interp_factor, debug, len_key));
    }

    /*
     * The private constructor
     */
    matrix_transpose_impl::matrix_transpose_impl(int input_len, int output_len, int interp_factor, bool debug, std::string len_key)
      : gr::tagged_stream_block("matrix_transpose",
				gr::io_signature::make(1, 1, sizeof(gr_complex)*input_len),
				gr::io_signature::make(1, 1, sizeof(gr_complex)*output_len*interp_factor), len_key),
				d_input_len(input_len),
				d_output_len(output_len),
				d_interp_factor(interp_factor),
				d_debug(debug)
    {
        set_relative_rate((double)input_len/(double)output_len);
        set_tag_propagation_policy(TPP_DONT);

    }

    /*
     * Our virtual destructor.
     */
    matrix_transpose_impl::~matrix_transpose_impl()
    {
    }

    int
    matrix_transpose_impl::calculate_output_stream_length(const gr_vector_int &ninput_items)
    {
        int noutput_items = d_input_len;
        return noutput_items ;
    }

    int
    matrix_transpose_impl::work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
        const gr_complex *in = (const gr_complex *) input_items[0];
        gr_complex *out = (gr_complex *) output_items[0];
        
        dout << "[MATRIX TRANSPOSE] ninput_items " << ninput_items[0] << " noutput_items " << noutput_items << std::endl;


        // Error handling
        if(ninput_items[0]*float(d_input_len)/float(d_output_len)-ninput_items[0]*d_input_len/d_output_len!=0) 
            throw std::runtime_error("[MATRIX TRANSPOSE] input_len and output_len do not match to packet length");
        
        
        if (pc_output_buffers_full(0) > 0.001)
        {
            return 0;
        }

        // Set noutput items
        noutput_items = d_input_len;
        
        // Update len key tag
        update_length_tags(noutput_items,0);
        
        memset(out, 0, sizeof(gr_complex)*d_interp_factor*d_output_len*d_input_len);

        // Reorganize samples
        for(int l = 0; l < d_input_len; l++){ // go through single input vector
            for(int k=0; k < ninput_items[0]; k++){ // go through all input vectors
                out[l*d_output_len*d_interp_factor + k] = in[k*d_input_len + l];
            }
        }
        dout << "[MATRIX TRANSPOSE] pc_work_time_total: " << this->pc_work_time_total()  << std::endl;
        dout << "[MATRIX TRANSPOSE] pc_input_buffers_full: " << pc_input_buffers_full(0) << "  pc_output_buffers_full: "  << pc_output_buffers_full(0) << std::endl;

        // Tell runtime system how many output items we produced.
        return noutput_items;
    }

  } /* namespace mimo_ofdm_jrc */
} /* namespace gr */

