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
#include "zero_pad_impl.h"
#include <random> 

namespace gr {
  namespace mimo_ofdm_jrc {

    zero_pad::sptr
    zero_pad::make(bool debug, unsigned int pad_front, unsigned int pad_tail)
    {
      return gnuradio::get_initial_sptr
        (new zero_pad_impl(debug, pad_front, pad_tail ));
    }


    /*
     * The private constructor
     */
    zero_pad_impl::zero_pad_impl(bool debug, unsigned int pad_front, unsigned int pad_tail)
      : gr::tagged_stream_block("zero_pad",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(1, 1, sizeof(gr_complex)), "packet_len"),
            d_debug(debug),
			d_pad_front(pad_front),
			d_pad_tail(pad_tail)
    {
        set_tag_propagation_policy(block::TPP_DONT);
    }

    /*
     * Our virtual destructor.
     */
    zero_pad_impl::~zero_pad_impl()
    {
    }

    int
    zero_pad_impl::calculate_output_stream_length(const gr_vector_int &ninput_items)
    {
      return ninput_items[0] + d_pad_front + d_pad_tail;
    }

    int
    zero_pad_impl::work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
        const gr_complex *in = (const gr_complex*)input_items[0];
        gr_complex *out = (gr_complex*)output_items[0];

        std::random_device rd;
        std::default_random_engine rand_generator(rd());
        std::normal_distribution<float> normal_distribution( 0.0, 1e-2 );
        
        std::vector<gr_complex> random_vector(d_pad_front + d_pad_tail);

        std::generate(random_vector.begin(), random_vector.end(),
        [&rand_generator, &normal_distribution]()-> gr_complex
        { 
          return gr_complex(normal_distribution(rand_generator), normal_distribution(rand_generator));
        });

        std::memcpy(out, &random_vector[0], sizeof(gr_complex)*d_pad_front);
        std::memcpy(out+d_pad_front, in, sizeof(gr_complex)*ninput_items[0]);
        std::memcpy(out+d_pad_front+ninput_items[0], &random_vector[d_pad_front], sizeof(gr_complex)*d_pad_tail);

        int produced = ninput_items[0] + d_pad_front + d_pad_tail;
        return produced;
    }

  } /* namespace mimo_ofdm_jrc */
} /* namespace gr */

