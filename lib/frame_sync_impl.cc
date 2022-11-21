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
#include "frame_sync_impl.h"

bool compare_abs2(const std::pair<gr_complex, int>& first, const std::pair<gr_complex, int>& second) 
{
	return abs(std::get<0>(first)) > abs(std::get<0>(second));
}

namespace gr {
  namespace mimo_ofdm_jrc {

    frame_sync::sptr
    frame_sync::make(int fft_len, int cp_len, unsigned int sync_length, std::vector<gr_complex> ltf_seq_time, bool debug)
    {
      return gnuradio::get_initial_sptr
        (new frame_sync_impl(fft_len, cp_len, sync_length, ltf_seq_time, debug));
    }


    /*
     * The private constructor
     */
    frame_sync_impl::frame_sync_impl(int fft_len, int cp_len, unsigned int sync_length, std::vector<gr_complex> ltf_seq_time, bool debug)
      : gr::block("frame_sync",
              gr::io_signature::make2(2, 2, sizeof(gr_complex), sizeof(gr_complex)),
              gr::io_signature::make(1, 1, sizeof(gr_complex))),
                d_fft_len(fft_len),
                d_cp_len(cp_len),
                d_ltf_fir(gr::filter::kernel::fir_filter_ccc(1, ltf_seq_time)),
                d_debug(debug),
                sample_offset(0),
                d_state(SYNC),
                SYNC_LENGTH(sync_length)
    {
        set_tag_propagation_policy(block::TPP_DONT);
		d_correlation = gr::fft::malloc_complex(8192);
        frame_tag_offset = -1;
    }

    /*
     * Our virtual destructor.
     */
    frame_sync_impl::~frame_sync_impl()
    {
        gr::fft::free(d_correlation);
    }

    void
    frame_sync_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
		if(d_state == SYNC) {
			ninput_items_required[0] = d_fft_len+d_cp_len;
			ninput_items_required[1] = d_fft_len+d_cp_len;

		} else {
			ninput_items_required[0] = noutput_items;
			ninput_items_required[1] = noutput_items;

            // ninput_items_required[0] = d_fft_len+d_cp_len;
			// ninput_items_required[1] = d_fft_len+d_cp_len;
		}        
    }

    int
    frame_sync_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
		gr::thread::scoped_lock lock(d_mutex);

	    const gr_complex *in = (const gr_complex*)input_items[0];
		const gr_complex *in_delayed = (const gr_complex*)input_items[1];
		gr_complex *out = (gr_complex*)output_items[0];

		gr::block::set_thread_priority(49);

		int ninput = std::min(std::min(ninput_items[0], ninput_items[1]), 8192);

		const uint64_t nread = nitems_read(0);
        
        dout << "===================================================================" << std::endl; 
        dout << "[FRAME SYNC] ninput_items[0]: " << ninput_items[0] << ", ninput_items[1]: " <<
				ninput_items[1] << ", nitems_read(0): " << nitems_read(0) << ",  noutput_items: " << noutput_items <<
				",  state: " << d_state << ", sample_offset: " << sample_offset << std::endl;    


		get_tags_in_range(d_tags, 0, nread, nread + ninput);
		if (d_tags.size()) 
        {
			std::sort(d_tags.begin(), d_tags.end(), gr::tag_t::offset_compare);

			frame_tag_offset = d_tags.front().offset; //const uint64_t 
			if(frame_tag_offset > nread) 
			{
				ninput = frame_tag_offset - nread;
                dout << "[FRAME SYNC] New frame tag received --> Still copying old samples: " << ninput << std::endl;
			} 
			else 
			{
				if(sample_offset && (d_state == SYNC)) 
				{
					throw std::runtime_error("[FRAME SYNC] Something is wrong!");
				}
				if(d_state == COPY) 
				{
					d_state = RESET;
				}
                else{
                    dout << std::endl;
                    dout << "[FRAME SYNC] New frame tag received --> Starting frame_tag_offset : " << frame_tag_offset << ", ninput: " << ninput << ", d_cfo_coarse_est: " << d_cfo_coarse_est << std::endl;
                }
				d_cfo_coarse_est = pmt::to_double(d_tags.front().value);
			}
		}

		int n_in = 0;
		int n_out = 0;

		switch(d_state) 
        {
		case SYNC:
			d_ltf_fir.filterN(d_correlation, in, std::min(SYNC_LENGTH, std::max(ninput - d_fft_len - 1, 0)));

            dout << "[FRAME SYNC] Current correlator length : " << std::min(SYNC_LENGTH, std::max(ninput - d_fft_len - 1, 0)) << std::endl;

			while(n_in + d_fft_len-1 < ninput) 
			{
				d_cor.push_back(std::pair<gr_complex, int>(d_correlation[n_in], sample_offset));
                // dout << "[FRAME SYNC] d_cor: " << d_cor.size() << ", n_in : " << n_in << " sample_offset: " << sample_offset << std::endl;

				n_in++;
				sample_offset++;

				if(sample_offset == SYNC_LENGTH) {
                    dout << "[FRAME SYNC] Enough samples collected --> Searching frame start..." << std::endl;

					search_frame_start();
                    dout << "[FRAME SYNC] d_frame_start: " << d_frame_start << std::endl;

					sample_offset = 0;
					total_out_count = 0;
					d_state = COPY;

					break;
				}
			}

			break;

		case COPY:
			while(n_in < ninput && n_out < noutput_items) 
			{
				int rel = sample_offset - d_frame_start;

				if(!rel)  
				{
					add_item_tag(0, nitems_written(0),
							pmt::string_to_symbol("frame_start"),
							pmt::from_double(d_cfo_coarse_est - d_freq_offset), //d_cfo_coarse_est - d_freq_offset
							pmt::string_to_symbol(name()));
                    dout << "[FRAME SYNC] Frame start tag added at " << nitems_written(0) << ", freq_offset tag: "<< d_cfo_coarse_est - d_freq_offset << std::endl;
				}

				if(rel >= 0 && (rel < d_fft_len*2 || ((rel - d_fft_len*2) % (d_fft_len+d_cp_len)) > d_cp_len-1)) 
				{
					out[n_out] = in_delayed[n_in] * exp(gr_complex(0, sample_offset * d_freq_offset));
					n_out++;
				}
				n_in++;
				sample_offset++;
			}
			break;

		case RESET: { //Removes some samples from STF
			while(n_out < noutput_items) 
            {
				if( ((total_out_count + n_out) % d_fft_len) == 0 ) 
                {
					sample_offset = 0;
					d_state = SYNC;
					break;
				} 
                else 
                {
					out[n_out] = 0;
					n_out++;
				}
			}
			break;
		}
		}

		dout << "[FRAME SYNC] Total Work Time: " <<  pc_work_time_total() << std::endl;
		dout << "[FRAME SYNC] Input Buffer: " << pc_input_buffers_full(0) << ",  Output Buffer: "  << pc_output_buffers_full(0) << std::endl;

        total_out_count += n_out;

        consume(0, n_in);
        consume(1, n_in);
        dout << "[FRAME SYNC] produced: " << n_out << ", consumed: " << n_in << ", total_out_count: " << total_out_count << std::endl;

        return n_out;
    }
    
    void 
	frame_sync_impl::search_frame_start() 
	{
		// sort list (highest correlation first)
		assert(d_cor.size() == SYNC_LENGTH);
		d_cor.sort(compare_abs2);

		// copy list in vector for nicer access
		std::vector<std::pair<gr_complex, int> > vec(d_cor.begin(), d_cor.end());
		d_cor.clear();

		// in case we don't find anything use SYNC_LENGTH
		d_frame_start = SYNC_LENGTH;

		for(int i = 0; i < 3; i++) 
        {
			for(int k = i + 1; k < 4; k++) 
            {
				gr_complex first_corr_val;
				gr_complex second_corr_val;

				if( std::get<1>(vec[i]) > std::get<1>(vec[k]) ) 
                {
					first_corr_val = std::get<0>(vec[k]);
					second_corr_val = std::get<0>(vec[i]);
				} 
                else 
                {
					first_corr_val = std::get<0>(vec[i]);
					second_corr_val = std::get<0>(vec[k]);
				}

				int diff  = abs(std::get<1>(vec[i]) - std::get<1>(vec[k]));
				if( diff == d_fft_len ) 
                {
					d_frame_start = std::min(std::get<1>(vec[i]), std::get<1>(vec[k]));
					d_freq_offset = arg(first_corr_val * conj(second_corr_val)) / d_fft_len;
					// nice match found, return immediately
                    dout << "[FRAME SYNC] Perfect peaks found! --> d_freq_offset: " << d_freq_offset << ", peak diff= " << abs(std::get<0>(vec[i])-std::get<0>(vec[k])) << std::endl;
					return;
				} 
                else if( diff == d_fft_len-1 ) 
                {
					d_frame_start = std::min(std::get<1>(vec[i]), std::get<1>(vec[k]));
					d_freq_offset = std::arg(first_corr_val * conj(second_corr_val)) / (d_fft_len-1);
                    dout << "[FRAME SYNC] Peaks found --> d_freq_offset: " << d_freq_offset << " --> diff:" << diff << ", peak diff= " << abs(std::get<0>(vec[i])-std::get<0>(vec[k])) << std::endl;
				} 
                else if( diff == d_fft_len+1 ) 
                {
					d_frame_start = std::min(std::get<1>(vec[i]), std::get<1>(vec[k]));
					d_freq_offset = std::arg(first_corr_val * conj(second_corr_val)) / (d_fft_len+1);
                    dout << "[FRAME SYNC] Peaks found --> d_freq_offset: " << d_freq_offset << " --> diff:" << diff  << ", peak diff= " << abs(std::get<0>(vec[i])-std::get<0>(vec[k])) << std::endl;
				}
			}
		}
	}

  } /* namespace mimo_ofdm_jrc */
} /* namespace gr */

