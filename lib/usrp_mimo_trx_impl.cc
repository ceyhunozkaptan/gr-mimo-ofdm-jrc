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
#include "usrp_mimo_trx_impl.h"
#include <iostream>
#include <boost/algorithm/string.hpp>
// #include <thread>

namespace gr {
  namespace mimo_ofdm_jrc {

	usrp_mimo_trx::sptr
	usrp_mimo_trx::make(int N_mboard, int N_tx, int N_rx, int samp_rate, float center_freq, int num_delay_samps, bool debug, float update_period,
				std::string args, std::string clock_sources, std::string time_sources, 
                std::string antenna_tx, float gain_tx, float timeout_tx, float wait_tx, std::string wire_tx, 
				std::string antenna_rx, float gain_rx, float timeout_rx, float wait_rx, float lo_offset_rx, std::string wire_rx, 
				const std::string& len_key)
	{
		return gnuradio::get_initial_sptr
		(new usrp_mimo_trx_impl(N_mboard, N_tx, N_rx, samp_rate, center_freq, num_delay_samps, debug, update_period,
				args, clock_sources, time_sources, 
                antenna_tx, gain_tx, timeout_tx, wait_tx, wire_tx, 
				antenna_rx, gain_rx, timeout_rx, wait_rx, lo_offset_rx, wire_rx, 
				len_key));
	}

	/*
		* The private constructor
		*/
	usrp_mimo_trx_impl::usrp_mimo_trx_impl(int N_mboard, int N_tx, int N_rx, int samp_rate, float center_freq, int num_delay_samps, bool debug, float update_period,
				std::string args, std::string clock_sources, std::string time_sources, 
                std::string antenna_tx, float gain_tx, float timeout_tx, float wait_tx, std::string wire_tx, 
				std::string antenna_rx, float gain_rx, float timeout_rx, float wait_rx, float lo_offset_rx, std::string wire_rx, 
				const std::string& len_key)
		: gr::tagged_stream_block("usrp_mimo_trx",
				gr::io_signature::make(N_tx, N_tx, sizeof(gr_complex)),
				gr::io_signature::make(N_rx, N_rx, sizeof(gr_complex)), len_key),
				d_debug(debug)
	{

        d_N_mboard = N_mboard;
        d_N_tx = N_tx;
        d_N_rx = N_rx;

		d_samp_rate = samp_rate;
		d_center_freq = center_freq;
		d_num_delay_samps = num_delay_samps;

        d_out_buffer.resize(d_N_rx);
        d_out_recv_ptrs.resize(d_N_rx);
        for (size_t i_rx = 0; i_rx < d_N_rx; i_rx++)
        {
            d_out_buffer[i_rx].resize(0);
            d_out_recv_ptrs[i_rx] = &d_out_buffer[i_rx].front();
        }
        // d_out_buffer.resize(0);

        d_in_send_ptrs.resize(d_N_tx);

		prev_tx_time = 0.0;
		d_update_period = update_period;

		//===========================================================================================================
		//= Setup USRP TX ===========================================================================================
		//===========================================================================================================
		d_args = args;
		d_wire_tx = wire_tx;
        d_wire_rx = wire_rx;

        clock_sources.erase(std::remove_if(clock_sources.begin(), clock_sources.end(), isspace), clock_sources.end());
		std::cout << "clock_source_tx: " << clock_sources << std::endl;

        boost::split(d_clock_sources, clock_sources, boost::is_any_of(","));
        // for (auto i : d_clock_sources)
        // {
        //     std::cout << "clock_source_tx: " << i << std::endl;
        // }

        time_sources.erase(std::remove_if(time_sources.begin(), time_sources.end(), isspace), time_sources.end());
        boost::split(d_time_sources, time_sources, boost::is_any_of(","));

        antenna_tx.erase(std::remove_if(antenna_tx.begin(), antenna_tx.end(), isspace), antenna_tx.end());
        boost::split(d_antenna_tx, antenna_tx, boost::is_any_of(","));

        antenna_rx.erase(std::remove_if(antenna_rx.begin(), antenna_rx.end(), isspace), antenna_rx.end());
        boost::split(d_antenna_rx, antenna_rx, boost::is_any_of(","));

		d_gain_tx = gain_tx;
		d_timeout_tx = timeout_tx; // timeout for sending
		d_wait_tx = wait_tx; // secs to wait befor sending

		d_lo_offset_rx = lo_offset_rx;
		d_gain_rx = gain_rx;
		d_timeout_rx = timeout_rx; // timeout for receiving
		d_wait_rx = wait_rx; // secs to wait befor receiving

		// Setup USRP TX: args (addr,...)
		d_usrp_mimo = uhd::usrp::multi_usrp::make(d_args);

		std::cout << "Summary of USRP Device : " << std::endl << d_usrp_mimo->get_pp_string() << std::endl;
		std::cout << "Number of MBoards : " << std::endl << d_usrp_mimo->get_num_mboards() << std::endl;

        std::cout << "Using Master Clock Rate (TX): " << std::endl << d_usrp_mimo->get_master_clock_rate() << std::endl;
        // Setup USRP TX: clock adn time source
        for (int i_mboard = 0; i_mboard < d_N_mboard; i_mboard++ )
        {
            d_usrp_mimo->set_clock_source(d_clock_sources[i_mboard], i_mboard);  // Set TX clock, TX is master
            d_usrp_mimo->set_time_source(d_time_sources[i_mboard], i_mboard);  // Set TX time, TX is master 
            std::cout << "USRP Clock Source (Mboard-" << i_mboard << "): " << d_usrp_mimo->get_clock_source(i_mboard) << std::endl;
            std::cout << "USRP Time Source (Mboard-" << i_mboard << "): " << d_usrp_mimo->get_time_source(i_mboard) << std::endl;
        }

        d_usrp_mimo->set_tx_lo_export_enabled(true, "lo1", 0);
        d_usrp_mimo->set_rx_lo_export_enabled(true, "lo1", 0);

        for (int i_tx = 0; i_tx < d_N_tx; i_tx++ )
        {
            d_usrp_mimo->set_tx_lo_source("external", "lo1", i_tx);
        }

        for (int i_rx = 0; i_rx < d_N_rx; i_rx++ )
        {
            d_usrp_mimo->set_rx_lo_source("external", "lo1", i_rx);
        }

        // TODO Add parameters to block for UI control
        // LO Output Switch to loop back
        d_usrp_mimo->get_device()->get_tree()->access<bool>("blocks/0/Radio#0/dboard/tx_frontends/0/los/lo1/lo_distribution/LO_OUT_1/export").set(true);
        d_usrp_mimo->get_device()->get_tree()->access<bool>("blocks/0/Radio#0/dboard/rx_frontends/0/los/lo1/lo_distribution/LO_OUT_1/export").set(true);
        // LO Output Switch for other USRP N320
        d_usrp_mimo->get_device()->get_tree()->access<bool>("blocks/0/Radio#0/dboard/tx_frontends/0/los/lo1/lo_distribution/LO_OUT_0/export").set(true);
        d_usrp_mimo->get_device()->get_tree()->access<bool>("blocks/0/Radio#0/dboard/rx_frontends/0/los/lo1/lo_distribution/LO_OUT_0/export").set(true);

		// Setup USRP TX: sample rate
		std::cout << "Setting TX/RX Rate: " << d_samp_rate << std::endl;
		d_usrp_mimo->set_tx_rate(d_samp_rate);
        d_usrp_mimo->set_rx_rate(d_samp_rate);
		std::cout << "Actual RX Rate: " << d_usrp_mimo->get_rx_rate() << std::endl;
		std::cout << "Actual TX Rate: " << d_usrp_mimo->get_tx_rate() << std::endl;

		// Setup USRP TX: time sync
        // d_usrp_mimo->set_time_now(uhd::time_spec_t(0.0)); // Do set time on startup if not gpsdo is activated.
        
        d_usrp_mimo->set_time_next_pps(uhd::time_spec_t(0.0));
        boost::this_thread::sleep(boost::posix_time::milliseconds(500));
        
        d_usrp_mimo->clear_command_time();
        d_usrp_mimo->set_command_time(d_usrp_mimo->get_time_now() + uhd::time_spec_t(0.1)); //set cmd time for .1s in the future

		// Setup USRP TX: tune request
		d_tune_request_tx = uhd::tune_request_t(d_center_freq); 
		// d_tune_request_tx = uhd::tune_request_t(d_center_freq, 0, dsp_policy=uhd.tune_request.POLICY_MANUAL, dsp_freq=LO-d_center_freq, lo_freq_policy=uhd.tune_request.POLICY_MANUAL, lo_freq=LO);
		
        for (int i_tx = 0; i_tx < d_N_tx; i_tx++ )
        {
            d_usrp_mimo->set_tx_freq(d_tune_request_tx, i_tx);
        }
        boost::this_thread::sleep(boost::posix_time::milliseconds(150));
        d_usrp_mimo->clear_command_time();
        std::cout << "TX Frequency is set!" << std::endl;

        d_usrp_mimo->set_command_time(d_usrp_mimo->get_time_now() + uhd::time_spec_t(0.1)); //set cmd time for .1s in the future
            
		// Setup USRP RX: tune request
		d_tune_request_rx = uhd::tune_request_t(d_center_freq, d_lo_offset_rx); 

        for (int i_rx = 0; i_rx < d_N_rx; i_rx++ )
        {
            d_usrp_mimo->set_rx_freq(d_tune_request_rx, i_rx);
        }
        boost::this_thread::sleep(boost::posix_time::milliseconds(150));
        d_usrp_mimo->clear_command_time();
        std::cout << "RX Frequency is set!" << std::endl;

		// Setup USRP TX: gain
		set_tx_gain(d_gain_tx);
        set_rx_gain(d_gain_rx);
        std::cout << "TX/RX Gains are set!" << std::endl;

        // Setup USRP TX: antenna
        for (int i_tx = 0; i_tx < d_N_tx; i_tx++ )
        {
            std::cout << "TX Channel " << i_tx << " using " << d_antenna_tx[i_tx] << std::endl;
            d_usrp_mimo->set_tx_antenna(d_antenna_tx[i_tx], i_tx);
        }

        for (int i_rx = 0; i_rx < d_N_rx; i_rx++ )
        {
            std::cout << "RX Channel " << i_rx << " using " << d_antenna_rx[i_rx] << std::endl;
            d_usrp_mimo->set_rx_antenna(d_antenna_rx[i_rx], i_rx);
        }

        // d_usrp_mimo->set_tx_dc_offset(0.02, 0);
        // d_usrp_mimo->set_tx_dc_offset(0.05, 1);

		// Setup transmit streamer
		uhd::stream_args_t stream_args_tx("fc32", d_wire_tx); // complex floats
        for (int i_tx = 0; i_tx < d_N_tx; i_tx++ )
        {
            stream_args_tx.channels.push_back(i_tx);
        }
		d_tx_streamer = d_usrp_mimo->get_tx_stream(stream_args_tx);
        std::cout << "Total TX Channels: " << d_usrp_mimo->get_tx_num_channels() << std::endl;
        std::cout << "Used TX Channels: " << d_tx_streamer->get_num_channels() << std::endl;

		// Setup receive streamer
		uhd::stream_args_t stream_args_rx("fc32", d_wire_rx); // complex floats
        for (int i_rx = 0; i_rx < d_N_rx; i_rx++ )
        {
            stream_args_rx.channels.push_back(i_rx);
        }
		// std::vector<size_t> channel_nums; channel_nums.push_back(0); // define channel!
		// stream_args_rx.channels = channel_nums;
		d_rx_streamer = d_usrp_mimo->get_rx_stream(stream_args_rx);
        std::cout << "Total RX Channels: " << d_usrp_mimo->get_rx_num_channels() << std::endl;
        std::cout << "Used RX Channels: " << d_rx_streamer->get_num_channels() << std::endl;


		//===========================================================================================================
		//= Other Setup =============================================================================================
		//===========================================================================================================

        std::vector<std::string> tree_list = d_usrp_mimo->get_tree()->list("blocks/0/Radio#0/dboard/tx_frontends/0/los/lo1/lo_distribution");
        for (auto i : tree_list)
        {
            std::cout << "USRP get_tree: " << i << std::endl;
        }
        
		// Setup rx_time pmt
		d_time_key = pmt::string_to_symbol("rx_time");
		d_srcid = pmt::string_to_symbol("usrp_mimo_trx");

		// Setup thread priority
		// uhd::set_thread_priority_safe(); // necessary? doesnt work...

		// Sleep to get sync done
		boost::this_thread::sleep(boost::posix_time::milliseconds(200)); // FIXME: necessary?

        std::vector<uhd::time_spec_t> current_time_specs_tx(d_N_mboard);
        std::vector<uhd::time_spec_t> current_time_specs_rx(d_N_mboard);

        for (int i_mboard = 0; i_mboard < d_N_mboard; i_mboard++ )
        {
            current_time_specs_tx[i_mboard] = d_usrp_mimo->get_time_now(i_mboard);
            // current_time_specs_rx[i_mboard] = d_usrp_rx->get_time_now(i_mboard);
        }

        for (int i_mboard = 0; i_mboard < d_N_mboard; i_mboard++ )
        {
            double current_time_tx = current_time_specs_tx[i_mboard].get_full_secs()+current_time_specs_tx[i_mboard].get_frac_secs();
            // double current_time_rx = current_time_specs_rx[i_mboard].get_full_secs()+current_time_specs_rx[i_mboard].get_frac_secs();

            std::cout << "USRP Current Time Source (Mboard-" << i_mboard << "): " << current_time_tx << std::endl;
            // std::cout << "USRP Current Time Source (RX-" << i_mboard << "): " << current_time_rx << std::endl;
        }
	}

	/*
		* Our virtual destructor.
		*/
	usrp_mimo_trx_impl::~usrp_mimo_trx_impl()
	{
	}

	int
	usrp_mimo_trx_impl::work (int noutput_items,
						gr_vector_int &ninput_items,
						gr_vector_const_void_star &input_items, // hold vector of pointers for multiple channels
						gr_vector_void_star &output_items) // hold vector of pointers for multiple channels
	{
		// gr_complex *in = (gr_complex *) input_items[0]; // remove const
		
		// std::cout  << "[USRP] active_thread_priority : " << gr::block::active_thread_priority()	<< std::endl;
		// std::cout  << "[USRP] thread_priority : " << gr::block::thread_priority()	<< std::endl;
		gr::block::set_thread_priority(70);
        boost::recursive_mutex::scoped_lock lock(d_mutex);

		// Set output items on packet length
		noutput_items = ninput_items[0];
        dout << "[USRP] ninput_items[0]: " << ninput_items[0] << std::endl;
        // dout << "[USRP] ninput_items[1]: " << ninput_items[1] << std::endl;
        // dout << "[USRP] ninput_items[2]: " << ninput_items[2] << std::endl;
        // dout << "[USRP] ninput_items[3]: " << ninput_items[3] << std::endl;
        dout << "[USRP] noutput_items: " << noutput_items << std::endl;

        dout << "[USRP] input_items.size(): " << input_items.size() << std::endl;
        dout << "[USRP] output_items.size(): " << output_items.size() << std::endl;

        for (size_t i_rx = 0; i_rx < d_N_rx; i_rx++)
        {
            if(d_out_buffer[i_rx].size() != noutput_items)
            {
                d_out_buffer[i_rx].resize(noutput_items);
                dout << "[USRP] Receive Buffer " << i_rx << " resized to " << noutput_items << std::endl;
                d_out_recv_ptrs[i_rx] = &d_out_buffer[i_rx].front();
            }
        }

        dout << "[USRP] d_out_recv_ptrs.size: " << d_out_recv_ptrs.size() << std::endl;

        for (size_t i_tx = 0; i_tx < d_N_tx; i_tx++)
        {
            // d_in_send_ptrs.push_back(input_items[i_tx]);
            d_in_send_ptrs[i_tx] = input_items[i_tx];
        }


		// Get time from USRP TX
		// d_time_now_tx = d_usrp_mimo->get_time_now();
		uhd::time_spec_t current_time_spec = d_usrp_mimo->get_time_now();
		double current_time = current_time_spec.get_full_secs()+current_time_spec.get_frac_secs();

		if (current_time >= prev_tx_time + d_update_period) // TX/RX MODE
		{
            dout << "[USRP] Performing TX/RX " << std::endl;
            prev_tx_time = current_time_spec.get_full_secs()+current_time_spec.get_frac_secs();
                
            d_noutput_items_tx = noutput_items;
            d_noutput_items_rx = noutput_items;
                
            d_time_now_tx = d_usrp_mimo->get_time_now();
            d_time_now_rx = d_time_now_tx;

            // Trasnmit thread
            d_thread_send = gr::thread::thread(boost::bind(&usrp_mimo_trx_impl::transmit, this));
            // Receive thread
            d_thread_recv = gr::thread::thread(boost::bind(&usrp_mimo_trx_impl::receive, this));
            
            // Wait for threads to complete
            d_thread_send.join();
            d_thread_recv.join();
            dout << "[USRP] TX/RX threads are joined!" << std::endl;

        }
        else
        {
            dout << "[USRP] Performing only TX " << std::endl;

            d_noutput_items_tx = noutput_items;
			d_noutput_items_rx = 0;
            
            // Send thread
            d_thread_send = gr::thread::thread(boost::bind(&usrp_mimo_trx_impl::transmit, this));

            d_thread_send.join();
			
			return 0;
        }
        
        gr_complex *out;

        for (int i_rx = 0; i_rx < d_N_rx; i_rx++) 
        {
            out = (gr_complex*) output_items[i_rx];
            
            memcpy(out, &d_out_buffer[i_rx][0] + d_num_delay_samps, (noutput_items - d_num_delay_samps)*sizeof(gr_complex)); // push buffer to output
            memset(out + (noutput_items-d_num_delay_samps), 0, d_num_delay_samps*sizeof(gr_complex)); // set zeros
            
            // Setup rx_time tag
            add_item_tag(i_rx, nitems_written(0), d_time_key, d_time_val, d_srcid);
        }
        // dout << "[USRP] Receive buffers are copied to output!" << std::endl;

		// Tell runtime system how many output items we produced.
		return noutput_items;
	}

	void
	usrp_mimo_trx_impl::transmit()
	{
	// Setup metadata for first package
		d_metadata_tx.start_of_burst = true;
		d_metadata_tx.end_of_burst = false;
		d_metadata_tx.has_time_spec = true;
		// std::cout << "[SEND] d_noutput_items_rx: " << d_noutput_items_rx << std::endl;

					
		if (d_noutput_items_rx == 0) // TX Only -> No RX Scheduled
		{
            d_time_now_tx = d_usrp_mimo->get_time_now();
            d_metadata_tx.time_spec = d_time_now_tx + uhd::time_spec_t(d_wait_tx); 
        }
        else
        {
            d_metadata_tx.time_spec = d_time_now_tx + uhd::time_spec_t(d_wait_tx);
        }

		// Send input buffer
		size_t num_tx_samps, total_num_samps;
		total_num_samps = d_noutput_items_tx;
        // dout << "[USRP] Transmitting samples: " << total_num_samps << std::endl;

        // Data to USRP
		num_tx_samps = d_tx_streamer->send(d_in_send_ptrs, total_num_samps, d_metadata_tx, total_num_samps/(float)d_samp_rate+d_timeout_tx);
        // num_tx_samps = d_tx_streamer->send(d_in_send, total_num_samps, d_metadata_tx, total_num_samps/(float)d_samp_rate+d_timeout_tx);
        
        dout << "[USRP] Transmission Done!: " << num_tx_samps << std::endl;


		// Get timeout
		if (num_tx_samps < total_num_samps){
            std::cerr << "Send timeout..." << std::endl;
        }

		//send a mini EOB packet
		d_metadata_tx.start_of_burst = false;
		d_metadata_tx.end_of_burst = true;
		d_metadata_tx.has_time_spec = false;
		d_tx_streamer->send("", 0, d_metadata_tx);
	}

	void
	usrp_mimo_trx_impl::receive()
	{
		// Setup RX streaming
		size_t total_num_samps = d_noutput_items_rx;
		uhd::stream_cmd_t stream_cmd(uhd::stream_cmd_t::STREAM_MODE_NUM_SAMPS_AND_DONE);
		// std::cout << "[USRP] Total Number of Samples for TX/RX : " << total_num_samps << std::endl;
		stream_cmd.num_samps = total_num_samps;
		stream_cmd.stream_now = false;
		stream_cmd.time_spec = d_time_now_rx + uhd::time_spec_t(d_wait_rx);

        // dout << "[USRP] Will receive samples: " << total_num_samps << std::endl;

		d_current_time_now = d_usrp_mimo->get_time_now();
		double current_time = d_current_time_now.get_full_secs() + d_current_time_now.get_frac_secs();
		double planned_rx_time = stream_cmd.time_spec.get_full_secs() + stream_cmd.time_spec.get_frac_secs();
        // dout << "[USRP] RX scheduled:" << planned_rx_time << std::endl;
        // dout << "[USRP] current_time: " << current_time << std::endl;

		if (planned_rx_time >= current_time)
		{
			d_rx_streamer->issue_stream_cmd(stream_cmd);
			size_t num_rx_samps;

            // dout << "[USRP] Receiving now: " << total_num_samps << std::endl;

			// Receive a packet
			num_rx_samps = d_rx_streamer->recv(d_out_recv_ptrs, total_num_samps, d_metadata_rx, total_num_samps/(float)d_samp_rate+d_timeout_rx);
            // num_rx_samps = d_rx_streamer->recv(d_out_recv, total_num_samps, d_metadata_rx, total_num_samps/(float)d_samp_rate+d_timeout_rx);
            dout << "[USRP] Receiving Done: " << num_rx_samps << std::endl;

			// Save timestamp
			// d_time_val = pmt::make_tuple(pmt::from_uint64(d_metadata_rx.time_spec.get_full_secs()),pmt::from_double(d_metadata_rx.time_spec.get_frac_secs()));
            d_time_val = pmt::from_double( d_metadata_rx.time_spec.get_full_secs() + d_metadata_rx.time_spec.get_frac_secs() );

			// Handle the error code
			if (d_metadata_rx.error_code != uhd::rx_metadata_t::ERROR_CODE_NONE)
			{
                std::cerr << "[USRP] Receiver Error: " << d_metadata_rx.strerror() << std::endl;
				// throw std::runtime_error(str(boost::format("Receiver error %s") % d_metadata_rx.strerror()));
			}

			if (num_rx_samps < total_num_samps)
            {
                std::cerr << "[USRP] Receive timeout before all samples received..." << std::endl;
                std::cerr << "[USRP] fragment_offset: " << d_metadata_rx.fragment_offset<< std::endl;
            }

            // std::cout << "[USRP] num_rx_samps : " << num_rx_samps << std::endl;
            // std::cout << "[USRP] error_code: " << d_metadata_rx.error_code << std::endl;
            // std::cout << "[USRP] fragment_offset: " << d_metadata_rx.fragment_offset<< std::endl;
            // std::cout << "[USRP] more_fragments: " << d_metadata_rx.more_fragments<< std::endl;
            // std::cout << "[USRP] out_of_sequence: " << d_metadata_rx.out_of_sequence<< std::endl;
		}
		else
		{
			std::cerr << "[USRP] current_time: " << current_time << std::endl;
			std::cerr <<  "[USRP] Planned RX Time: " << planned_rx_time << std::endl;
			std::cerr << "[USRP] Timing Requirements cannot Met!! -> Time Diff: " << current_time - planned_rx_time << std::endl;
			d_noutput_items_rx = 0;
		}
	}

	int usrp_mimo_trx_impl::calculate_output_stream_length(const gr_vector_int &ninput_items)
	{
		int noutput_items = ninput_items[0];
		return noutput_items ;
	}

	void usrp_mimo_trx_impl::set_num_delay_samps(int num_samps){
		d_num_delay_samps = num_samps;
	}

	void usrp_mimo_trx_impl::set_rx_gain(float gain){
        for (int i_rx = 0; i_rx < d_N_rx; i_rx++ )
        {
            d_usrp_mimo->set_rx_gain(gain, i_rx);
        }
	}

	void usrp_mimo_trx_impl::set_tx_gain(float gain){
        for (int i_tx = 0; i_tx < d_N_tx; i_tx++ )
        {
            d_usrp_mimo->set_tx_gain(gain, i_tx);
        }
	}

  } /* namespace mimo_ofdm_jrc */
} /* namespace gr */

