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

#include <QTimer>
#include <QWidget>

#include <qwt_plot.h>
#include <qwt_plot_spectrogram.h>
#include <qwt_matrix_raster_data.h>
#include <qwt_color_map.h>
#include "qcolormap.h"

#include <QApplication>
#include <qwt_scale_widget.h>
#include <qwt_scale_engine.h>
#include <complex>

#include <qwt_point_3d.h>
#include "range_angle_raster_data.h"

namespace gr {
  namespace mimo_ofdm_jrc {

    class heatmap_plot : public QWidget
{
    Q_OBJECT

    public:
        heatmap_plot(int interval,
                        int vlen,
                        std::vector<float>* buffer,
                        std::string label_x,
                        std::string label_y,
                        std::string label,
                        std::vector<float> axis_x,
                        std::vector<float> axis_y,
                        float dynamic_range_db,
                        std::vector<float> x_axis_ticks,
                        std::vector<float> y_axis_ticks,
                        bool autoscale_z,
                        bool db_scale,
                        QWidget* parent = 0);
        ~heatmap_plot();

    private:
        int d_interval, d_vlen;

        int d_x_axis_div, d_y_axis_div;
        QVector<double> d_axis_x, d_axis_y;
        QVector<double> d_x_axis_ticks, d_y_axis_ticks;
        float d_dynamic_range_db;

        std::vector<float>* d_buffer;
        bool d_autoscale_z;
        bool d_db_scale;
        QTimer* d_timer;

        QwtPlot* d_plot;
        QwtPlotSpectrogram* d_spectrogram;

        RangeAngleRasterData* d_data;

        QwtLinearColorMap* d_colormap;
        QwtScaleWidget* d_scale;

        QVector<double> d_plot_data;
        QVector<QwtPoint3D> d_plot_data_3D;

        QFont axis_font;
        QFont axis_title_font;

    protected:
        void resizeEvent(QResizeEvent* event);

    public slots:
        void refresh();
    };

  } // namespace mimo_ofdm_jrc
} // namespace gr

