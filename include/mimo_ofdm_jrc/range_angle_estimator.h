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

#ifndef INCLUDED_MIMO_OFDM_JRC_RANGE_ANGLE_ESTIMATOR_H
#define INCLUDED_MIMO_OFDM_JRC_RANGE_ANGLE_ESTIMATOR_H

#include <mimo_ofdm_jrc/api.h>
#include <gnuradio/tagged_stream_block.h>

namespace gr {
  namespace mimo_ofdm_jrc {

    /*!
     * \brief <+description of block+>
     * \ingroup mimo_ofdm_jrc
     *
     */
    class MIMO_OFDM_JRC_API range_angle_estimator : virtual public gr::tagged_stream_block
    {
     public:
      typedef boost::shared_ptr<range_angle_estimator> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of mimo_ofdm_jrc::range_angle_estimator.
       *
       * To avoid accidental use of raw pointers, mimo_ofdm_jrc::range_angle_estimator's
       * constructor is in a private implementation
       * class. mimo_ofdm_jrc::range_angle_estimator::make is the public interface for
       * creating new instances.
       */
      static sptr make(int vlen, 
                        std::vector<float> range_bins,
                        std::vector<float> angle_bins,
                        float noise_discard_range_m,
                        float noise_discard_angle_deg,
                        float snr_threshold, 
                        float power_threshold, 
                        const std::string& stats_path,
                        bool stats_record, 
                        const std::string& len_key = "packet_len",
                        bool debug = false);
     
     virtual void set_snr_threshold(float snr_threshold) = 0;
     virtual void set_power_threshold(float power_threshold) = 0;
     virtual void set_stats_record(bool stats_record) = 0;
     
    };

  } // namespace mimo_ofdm_jrc
} // namespace gr

#endif /* INCLUDED_MIMO_OFDM_JRC_RANGE_ANGLE_ESTIMATOR_H */

