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

#ifndef INCLUDED_MIMO_OFDM_JRC_GUI_TIME_PLOT_IMPL_H
#define INCLUDED_MIMO_OFDM_JRC_GUI_TIME_PLOT_IMPL_H

#include <mimo_ofdm_jrc/gui_time_plot.h>
#include "time_plot.h"

namespace gr {
  namespace mimo_ofdm_jrc {

    class gui_time_plot_impl : public gui_time_plot
    {
     private:
      // Nothing to declare in this block.

     public:
      gui_time_plot_impl(int interval,
                          std::string y_tag,
                          std::string y_label,
                          std::vector<float> y_axis_limits,
                          float time_window,
                          std::string title_label);
      ~gui_time_plot_impl();
      void handle_msg(pmt::pmt_t msg);
      void run_gui();

      int d_interval;
      std::string d_y_label, d_y_tag;
      std::vector<float> d_y_axis_limits;
      float d_time_window;
      pmt::pmt_t d_port_id_in;
      std::vector<float> d_y;
      bool d_y_read;
      std::string d_title_label;

      int d_argc;
      char* d_argv;
      QApplication* d_qApplication;

      time_plot* d_main_gui;

    };

  } // namespace mimo_ofdm_jrc
} // namespace gr

#endif /* INCLUDED_MIMO_OFDM_JRC_GUI_TIME_PLOT_IMPL_H */

