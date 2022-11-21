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

#ifndef INCLUDED_MIMO_OFDM_JRC_RANGE_ANGLE_ESTIMATOR_IMPL_H
#define INCLUDED_MIMO_OFDM_JRC_RANGE_ANGLE_ESTIMATOR_IMPL_H

#include <mimo_ofdm_jrc/range_angle_estimator.h>
#include "utils.h"

namespace gr {
  namespace mimo_ofdm_jrc {

    class range_angle_estimator_impl : public range_angle_estimator
    {
     private:
      	int d_vlen;
        std::vector<float> d_range_bins;
        std::vector<float> d_angle_bins;
        float d_noise_discard_range_m;
        float d_noise_discard_angle_deg;
        float d_snr_threshold;
        float d_power_threshold; 
        bool d_stats_record; 
        bool d_debug;
        const std::string d_stats_path;
        bool d_new_stat_started;

        std::ofstream file_stream;

     protected:
      int calculate_output_stream_length(const gr_vector_int &ninput_items);

     public:
      range_angle_estimator_impl(int vlen, 
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
      ~range_angle_estimator_impl();
      
      void set_snr_threshold(float snr_threshold);
      void set_power_threshold(float power_threshold);
      void set_stats_record(bool stats_record);

      // Where all the action really happens
      int work(
              int noutput_items,
              gr_vector_int &ninput_items,
              gr_vector_const_void_star &input_items,
              gr_vector_void_star &output_items
      );
    };

  } // namespace mimo_ofdm_jrc
} // namespace gr

#endif /* INCLUDED_MIMO_OFDM_JRC_RANGE_ANGLE_ESTIMATOR_IMPL_H */

