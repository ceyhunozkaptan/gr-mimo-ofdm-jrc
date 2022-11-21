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
#include "gui_heatmap_plot_impl.h"

namespace gr {
  namespace mimo_ofdm_jrc {

    gui_heatmap_plot::sptr 
    gui_heatmap_plot::make(int vlen,
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
                                                    std::string len_key)
    {
        return gnuradio::get_initial_sptr(new gui_heatmap_plot_impl(vlen,
                                                                    interval,
                                                                    xlabel,
                                                                    ylabel,
                                                                    label,
                                                                    axis_x,
                                                                    axis_y,
                                                                    dynamic_range_db,
                                                                    x_axis_ticks,
                                                                    y_axis_ticks,
                                                                    autoscale_z,
                                                                    db_scale,
                                                                    len_key));
    }

    /*
    * The private constructor
    */
    gui_heatmap_plot_impl::gui_heatmap_plot_impl(int vlen,
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
                                                    std::string len_key)
        : gr::tagged_stream_block("gui_heatmap_plot",
                                gr::io_signature::make(1, 1, sizeof(float) * vlen),
                                gr::io_signature::make(0, 0, 0),
                                len_key)
    {
        d_vlen = vlen;
        d_interval = interval;
        d_xlabel = xlabel;
        d_ylabel = ylabel;
        d_label = label;
        d_axis_x = axis_x;
        d_axis_y = axis_y;
        // d_axis_z = axis_z;
        d_dynamic_range_db = dynamic_range_db;
        d_x_axis_ticks = x_axis_ticks;
        d_y_axis_ticks = y_axis_ticks;
        d_autoscale_z = autoscale_z;
        d_db_scale = db_scale;
        d_buffer.resize(0);

        // Setup GUI
        run_gui();
    }

    /*
    * Our virtual destructor.
    */
    gui_heatmap_plot_impl::~gui_heatmap_plot_impl() {}

    void gui_heatmap_plot_impl::run_gui()
    {
        // Set QT window
        if (qApp != NULL) {
            d_qApplication = qApp;
        } else {
            d_argc = 1;
            d_argv = new char;
            d_argv[0] = '\0';
            d_qApplication = new QApplication(d_argc, &d_argv);
        }

        // Set QWT plot widget
        d_main_gui = new heatmap_plot(d_interval,
                                        d_vlen,
                                        &d_buffer,
                                        d_xlabel,
                                        d_ylabel,
                                        d_label,
                                        d_axis_x,
                                        d_axis_y,
                                        d_dynamic_range_db,
                                        d_x_axis_ticks,
                                        d_y_axis_ticks,
                                        d_autoscale_z,
                                        d_db_scale);
        d_main_gui->show();
    }

    int gui_heatmap_plot_impl::calculate_output_stream_length(
        const gr_vector_int& ninput_items)
    {
        int noutput_items = 0;
        return noutput_items;
    }

    int gui_heatmap_plot_impl::work(int noutput_items,
                                        gr_vector_int& ninput_items,
                                        gr_vector_const_void_star& input_items,
                                        gr_vector_void_star& output_items)
    {
        const float* in = (const float*)input_items[0];

        // Copy data to shared buffer with GUI
        if (d_buffer.size() != ninput_items[0] * d_vlen) { // resize buffer if needed
            d_buffer.resize(ninput_items[0] * d_vlen);
        }
        memcpy(&d_buffer[0], in, sizeof(float) * ninput_items[0] * d_vlen);

        // Tell runtime system how many output items we produced.
        return 0;
    }

  } /* namespace mimo_ofdm_jrc */
} /* namespace gr */

