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

#ifndef INCLUDED_MIMO_OFDM_JRC_OFDM_FRAME_GENERATOR_H
#define INCLUDED_MIMO_OFDM_JRC_OFDM_FRAME_GENERATOR_H

#include <mimo_ofdm_jrc/api.h>
#include <gnuradio/tagged_stream_block.h>

namespace gr {
  namespace mimo_ofdm_jrc {

    class MIMO_OFDM_JRC_API ofdm_frame_generator : virtual public gr::tagged_stream_block
    {
     public:
      typedef boost::shared_ptr<ofdm_frame_generator> sptr;
      virtual std::string len_tag_key() = 0;
      virtual const int fft_len() = 0;
      virtual std::vector<std::vector<int>> occupied_carriers() = 0;

      /*!
       * \brief Return a shared_ptr to a new instance of mimo_ofdm_jrc::ofdm_frame_generator.
       *
       * To avoid accidental use of raw pointers, mimo_ofdm_jrc::ofdm_frame_generator's
       * constructor is in a private implementation
       * class. mimo_ofdm_jrc::ofdm_frame_generator::make is the public interface for
       * creating new instances.
       */
      static sptr make(int fft_len,
                     const std::vector<std::vector<int>>& occupied_carriers,
                     const std::vector<std::vector<int>>& pilot_carriers,
                     const std::vector<std::vector<gr_complex>>& pilot_symbols,
                     const std::vector<std::vector<gr_complex>>& sync_words,
                     int ltf_len,
                     const std::string& len_tag_key = "packet_len",
                     const bool output_is_shifted = true);
    };

  } // namespace mimo_ofdm_jrc
} // namespace gr

#endif /* INCLUDED_MIMO_OFDM_JRC_OFDM_FRAME_GENERATOR_H */

