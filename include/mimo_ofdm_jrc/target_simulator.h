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

#ifndef INCLUDED_MIMO_OFDM_JRC_TARGET_SIMULATOR_H
#define INCLUDED_MIMO_OFDM_JRC_TARGET_SIMULATOR_H

#include <mimo_ofdm_jrc/api.h>
#include <gnuradio/tagged_stream_block.h>

namespace gr {
  namespace mimo_ofdm_jrc {

    /*!
     * \brief <+description of block+>
     * \ingroup mimo_ofdm_jrc
     *
     */
    class MIMO_OFDM_JRC_API target_simulator : virtual public gr::tagged_stream_block
    {
     public:
      typedef boost::shared_ptr<target_simulator> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of mimo_ofdm_jrc::target_simulator.
       *
       * To avoid accidental use of raw pointers, mimo_ofdm_jrc::target_simulator's
       * constructor is in a private implementation
       * class. mimo_ofdm_jrc::target_simulator::make is the public interface for
       * creating new instances.
       */
      static sptr make(std::vector<float> range,
                     std::vector<float> velocity,
                     std::vector<float> rcs,
                     std::vector<float> azimuth,
                     std::vector<float> position_rx,
                     int samp_rate,
                     float center_freq,
                     float self_coupling_db,
                     bool rndm_phaseshift = false,
                     bool self_coupling = false,
                     const std::string& len_key = "packet_len",
                     bool debug = false);
                     
       
      virtual void setup_targets(std::vector<float> range,
                                std::vector<float> velocity,
                                std::vector<float> rcs,
                                std::vector<float> azimuth,
                                std::vector<float> position_rx,
                                int samp_rate,
                                float center_freq,
                                float self_coupling_db,
                                bool rndm_phaseshift,
                                bool self_coupling) = 0;
    };

  } // namespace mimo_ofdm_jrc
} // namespace gr

#endif /* INCLUDED_MIMO_OFDM_JRC_TARGET_SIMULATOR_H */

