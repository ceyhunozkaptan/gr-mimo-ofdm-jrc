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

#ifndef INCLUDED_MIMO_OFDM_JRC_MIMO_OFDM_EQUALIZER_H
#define INCLUDED_MIMO_OFDM_JRC_MIMO_OFDM_EQUALIZER_H

#include <mimo_ofdm_jrc/api.h>
#include <gnuradio/block.h>

enum ChannelEstimator {
	LS   = 0,
	STA  = 1,
};

enum Modulation {
	BPSK  = 0,
	QPSK = 1,
	QAM16 = 2,
};

namespace gr {
  namespace mimo_ofdm_jrc {

    /*!
     * \brief <+description of block+>
     * \ingroup mimo_ofdm_jrc
     *
     */
    class MIMO_OFDM_JRC_API mimo_ofdm_equalizer : virtual public gr::block
    {
     public:
      typedef boost::shared_ptr<mimo_ofdm_equalizer> sptr;
      
      virtual void set_estimator(ChannelEstimator algo) = 0;
      virtual void set_bandwidth(double bw) = 0;
      virtual void set_frequency(double freq) = 0;
      virtual void set_stats_record(bool stats_record) = 0;

      /*!
       * \brief Return a shared_ptr to a new instance of mimo_ofdm_jrc::mimo_ofdm_equalizer.
       *
       * To avoid accidental use of raw pointers, mimo_ofdm_jrc::mimo_ofdm_equalizer's
       * constructor is in a private implementation
       * class. mimo_ofdm_jrc::mimo_ofdm_equalizer::make is the public interface for
       * creating new instances.
       */
      static sptr make(ChannelEstimator estimator_algo, 
                            double freq, 
                            double bw, 
                            int fft_len, 
                            int cp_len, 
                            std::vector<int> data_carriers, 
                            std::vector<int> pilot_carriers, 
                            const std::vector<std::vector<gr_complex>>& pilot_symbols,
                            std::vector<gr_complex> long_seq, 
                            const std::vector<std::vector<gr_complex>>& mapped_ltf_symbols,
                            int n_mimo_ltf,
                            const std::string& chan_est_file,
                            const std::string& comm_log_file,
                            bool stats_record,                            
                            bool debug);
    };

  } // namespace mimo_ofdm_jrc
} // namespace gr

#endif /* INCLUDED_MIMO_OFDM_JRC_MIMO_OFDM_EQUALIZER_H */

