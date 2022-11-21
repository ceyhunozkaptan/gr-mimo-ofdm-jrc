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
#include "moving_avg_impl.h"

namespace gr {
  namespace mimo_ofdm_jrc {

    moving_avg::sptr
    moving_avg::make(int length, float scale, int max_iter, bool debug)
    {
      return gnuradio::get_initial_sptr
        (new moving_avg_impl(length, scale, max_iter, debug));
    }

    /*
     * The private constructor
     */
    moving_avg_impl::moving_avg_impl(int length, float scale, int max_iter, bool debug)
      : gr::sync_block("moving_avg",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(1, 1, sizeof(gr_complex))),
        d_length(length),
        d_scale(scale),
        d_max_iter(max_iter),
        d_new_length(length),
        d_new_scale(scale),
        d_updated(false),
        d_debug(debug)
    {
        this->set_history(length);
    }

    /*
     * Our virtual destructor.
     */
    moving_avg_impl::~moving_avg_impl()
    {
    }

    int
    moving_avg_impl::work(int noutput_items,
        gr_vector_const_void_star &input_items,
        gr_vector_void_star &output_items)
    {

        gr::block::set_thread_priority(50);

        if (d_updated) 
        {
            d_length = d_new_length;
            d_scale = d_new_scale;
            this->set_history(d_length);
            d_updated = false;
            return 0; 
        }

        dout << "[MOVING AVG]: noutput_items: " << noutput_items << std::endl;

        const gr_complex* in = (const gr_complex*)input_items[0];
        gr_complex* out = (gr_complex*) output_items[0];

        unsigned int num_iter = (unsigned int)((noutput_items > d_max_iter) ? d_max_iter : noutput_items);

        gr_complex sum = in[0];
        for (int i = 1; i < d_length - 1; i++) 
        {
            sum += in[i];
        }

        for (unsigned int i = 0; i < num_iter; i++) 
        {
            sum += in[i + d_length - 1];
            out[i] = sum * d_scale;
            sum -= in[i];
        }

        return num_iter;
    }

    void moving_avg_impl::set_length_and_scale(int length, float scale)
    {
        d_new_length = length;
        d_new_scale = scale;
        d_updated = true;
    }

    void moving_avg_impl::set_length(int length)
    {
        d_new_length = length;
        d_updated = true;
    }

    void moving_avg_impl::set_scale(float scale)
    {
        d_new_scale = scale;
        d_updated = true;
    }

  } /* namespace mimo_ofdm_jrc */
} /* namespace gr */

