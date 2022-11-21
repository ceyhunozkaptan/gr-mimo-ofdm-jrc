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

#ifndef INCLUDED_MIMO_OFDM_JRC_GUI_HEATMAP_PLOT_IMPL_H
#define INCLUDED_MIMO_OFDM_JRC_GUI_HEATMAP_PLOT_IMPL_H

#include <mimo_ofdm_jrc/gui_heatmap_plot.h>
#include "heatmap_plot.h"

namespace gr {
  namespace mimo_ofdm_jrc {

    class gui_heatmap_plot_impl : public gui_heatmap_plot
    {
     private:
      // Nothing to declare in this block.

     protected:
      int calculate_output_stream_length(const gr_vector_int &ninput_items);

     public:
      gui_heatmap_plot_impl(int vlen,
                                int interval,
                                std::string xlabel,
                                std::string ylabel,
                                std::string label,
                                std::vector<float> axis_x,
                                std::vector<float> axis_y,
                                float dynamic_range_db,
                                std::vector<float> x_axis_ticks,
                                std::vector<float> y_axis_ticks,
                                bool autoscale_z,
                                bool db_scale,
                                std::string len_key);
      ~gui_heatmap_plot_impl();
      
        int d_argc;
        char* d_argv;
        QApplication* d_qApplication;

        int d_vlen, d_interval;
        std::string d_xlabel, d_ylabel, d_label;
        std::vector<float> d_axis_x, d_axis_y;
        float d_dynamic_range_db;
        std::vector<float> d_x_axis_ticks, d_y_axis_ticks;
        
        std::vector<float> d_buffer;
        bool d_autoscale_z, d_db_scale;
        heatmap_plot* d_main_gui;

        void run_gui();

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

#endif /* INCLUDED_MIMO_OFDM_JRC_GUI_HEATMAP_PLOT_IMPL_H */

