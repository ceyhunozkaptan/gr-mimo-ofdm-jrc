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

#ifndef INCLUDED_MIMO_OFDM_JRC_MOVING_AVG_H
#define INCLUDED_MIMO_OFDM_JRC_MOVING_AVG_H

#include <mimo_ofdm_jrc/api.h>
#include <gnuradio/sync_block.h>

namespace gr {
  namespace mimo_ofdm_jrc {

    /*!
     * \brief <+description of block+>
     * \ingroup mimo_ofdm_jrc
     *
     */
    class MIMO_OFDM_JRC_API moving_avg : virtual public gr::sync_block
    {
     public:
      typedef boost::shared_ptr<moving_avg> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of mimo_ofdm_jrc::moving_avg.
       *
       * To avoid accidental use of raw pointers, mimo_ofdm_jrc::moving_avg's
       * constructor is in a private implementation
       * class. mimo_ofdm_jrc::moving_avg::make is the public interface for
       * creating new instances.
       */
      static sptr make(int length, float scale, int max_iter, bool debug);
      
       virtual int length() const = 0;
       virtual float scale() const = 0;
       virtual void set_length_and_scale(int length, float scale) = 0;
       virtual void set_length(int length) = 0;
       virtual void set_scale(float scale) = 0;
    };

  } // namespace mimo_ofdm_jrc
} // namespace gr

#endif /* INCLUDED_MIMO_OFDM_JRC_MOVING_AVG_H */

