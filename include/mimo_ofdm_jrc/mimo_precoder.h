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

#ifndef INCLUDED_MIMO_OFDM_JRC_MIMO_PRECODER_H
#define INCLUDED_MIMO_OFDM_JRC_MIMO_PRECODER_H

#include <mimo_ofdm_jrc/api.h>
#include <gnuradio/tagged_stream_block.h>

namespace gr {
  namespace mimo_ofdm_jrc {

    /*!
     * \brief <+description of block+>
     * \ingroup mimo_ofdm_jrc
     *
     */
    class MIMO_OFDM_JRC_API mimo_precoder : virtual public gr::tagged_stream_block
    {
     public:
      typedef boost::shared_ptr<mimo_precoder> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of mimo_ofdm_jrc::mimo_precoder.
       *
       * To avoid accidental use of raw pointers, mimo_ofdm_jrc::mimo_precoder's
       * constructor is in a private implementation
       * class. mimo_ofdm_jrc::mimo_precoder::make is the public interface for
       * creating new instances.
       */
      static sptr make(int fft_len,
                        int N_tx,
                        int N_ss,
                        const std::vector<int>& data_carriers,
                        const std::vector<int>& pilot_carriers,
                        const std::vector<std::vector<gr_complex>>& pilot_symbols,
                        const std::vector<std::vector<gr_complex>>& sync_words,
                        const std::vector<std::vector<gr_complex>>& mapped_ltf_symbols,
                        const std::string& chan_est_file,
                        bool chan_est_smoothing,
                        const std::string& radar_log_file,
                        bool radar_aided,
                        bool phased_steering,
                        bool use_radar_streams,
                        const std::string& len_tag_key = "packet_len",
                        bool debug = false);
    
        virtual void set_chan_est_smoothing(bool chan_est_smoothing) = 0;
        virtual void set_radar_aided(bool radar_aided) = 0;
        virtual void set_use_radar_streams(bool use_radar_streams) = 0;
        virtual void set_phased_steering(bool phased_steering) = 0;
    };

  } // namespace mimo_ofdm_jrc
} // namespace gr

#endif /* INCLUDED_MIMO_OFDM_JRC_MIMO_PRECODER_H */

