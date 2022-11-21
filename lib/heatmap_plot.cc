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

#include "heatmap_plot.h"
#include <iostream>
#include <stdexcept>

namespace gr {
namespace mimo_ofdm_jrc {


    heatmap_plot::heatmap_plot(int interval,
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
                                    QWidget* parent)
        : QWidget(parent)
    {
        d_interval = interval;
        d_vlen = vlen;
        d_buffer = buffer;
        d_autoscale_z = autoscale_z;
        d_db_scale = db_scale;

        d_axis_x = QVector<double>::fromStdVector(std::vector<double>(axis_x.begin(), axis_x.end()));
        d_axis_y = QVector<double>::fromStdVector(std::vector<double>(axis_y.begin(), axis_y.end()));
        d_x_axis_ticks = QVector<double>::fromStdVector(std::vector<double>(x_axis_ticks.begin(), x_axis_ticks.end()));
        d_y_axis_ticks = QVector<double>::fromStdVector(std::vector<double>(y_axis_ticks.begin(), y_axis_ticks.end()));

        d_dynamic_range_db = dynamic_range_db;
        // Setup GUI
        resize(QSize(960, 960));
        setWindowTitle("Heatmap Plot");

        d_plot = new QwtPlot(this);               // make main plot
        d_spectrogram = new QwtPlotSpectrogram(); // make spectrogram
        d_spectrogram->attach(d_plot);            // attach spectrogram to plot

        d_data = new RangeAngleRasterData(); // make data structure

        d_colormap = new QwtLinearColorMap( QColor(0,0,189), QColor(132,0,0), QwtColorMap::RGB );
        double pos;
        pos = 1.0/13.0*1.0; d_colormap->addColorStop(pos, QColor(0,0,255));
        pos = 1.0/13.0*2.0; d_colormap->addColorStop(pos, QColor(0,66,255));
        pos = 1.0/13.0*3.0; d_colormap->addColorStop(pos, QColor(0,132,255));
        pos = 1.0/13.0*4.0; d_colormap->addColorStop(pos, QColor(0,189,255));
        pos = 1.0/13.0*5.0; d_colormap->addColorStop(pos, QColor(0,255,255));
        pos = 1.0/13.0*6.0; d_colormap->addColorStop(pos, QColor(66,255,189));
        pos = 1.0/13.0*7.0; d_colormap->addColorStop(pos, QColor(132,255,132));
        pos = 1.0/13.0*8.0; d_colormap->addColorStop(pos, QColor(189,255,66));
        pos = 1.0/13.0*9.0; d_colormap->addColorStop(pos, QColor(255,255,0));
        pos = 1.0/13.0*10.0; d_colormap->addColorStop(pos, QColor(255,189,0));
        pos = 1.0/13.0*12.0; d_colormap->addColorStop(pos, QColor(255,66,0));
        pos = 1.0/13.0*13.0; d_colormap->addColorStop(pos, QColor(189,0,0));

        d_spectrogram->setColorMap(d_colormap);

        QwtText title_text;
        QFont title_font;
        title_font.setPointSize(20); 
        title_font.setBold(true);
        title_text.setFont(title_font);
        title_text.setText(label.c_str());
        d_plot->setTitle(title_text);

        QwtText axis_title_text;
        axis_title_font.setPointSize(20); 
        axis_title_font.setBold(true);
        axis_title_text.setFont(axis_title_font);
        axis_title_text.setText(label_x.c_str());

        d_plot->setAxisTitle(QwtPlot::xBottom, axis_title_text);
        axis_title_text.setText(label_y.c_str());
        d_plot->setAxisTitle(QwtPlot::yLeft, axis_title_text);

        axis_font.setPointSize(15); 
        d_plot->setAxisFont(QwtPlot::xBottom, axis_font);
        d_plot->setAxisFont(QwtPlot::yLeft, axis_font);

        if (d_x_axis_ticks.size() != 3 || d_y_axis_ticks.size() != 3)
        {
            throw std::invalid_argument("[HEATMAP PLOT] x-axis and y-axis ticks should be as [start, stop, step size]!");
        }

        d_plot->setAxisScale(QwtPlot::xBottom, d_x_axis_ticks[0], d_x_axis_ticks[1], d_x_axis_ticks[2]);
        d_plot->setAxisScale(QwtPlot::yLeft, d_y_axis_ticks[0], d_y_axis_ticks[1], d_y_axis_ticks[2]);

        // Do replot
        d_plot->replot();

        // Setup timer and connect refreshing plot
        d_timer = new QTimer(this);
        connect(d_timer, SIGNAL(timeout()), this, SLOT(refresh()));
        d_timer->start(d_interval);
    }

    heatmap_plot::~heatmap_plot() {}

    void heatmap_plot::resizeEvent(QResizeEvent* event)
    {
        d_plot->setGeometry(0, 0, this->width(), this->height());
    }

    void heatmap_plot::refresh()
    {
        // Fetch new data and push to matrix
        d_plot_data.clear();

        float maximum, minimum; // get maximum and minimum
        float min_display;
        if (d_buffer->size() != 0) {
            d_plot_data.resize(d_buffer->size());

            // d_plot_data = QVector<double>(d_buffer->begin(), d_buffer->end());
            d_plot_data = QVector<double>::fromStdVector(std::vector<double>(d_buffer->begin(), d_buffer->end()));

            minimum = *std::min_element( std::begin(d_plot_data), std::end(d_plot_data) );
            maximum = *std::max_element( std::begin(d_plot_data), std::end(d_plot_data) );
            
            if (std::isnan(minimum) || std::isnan(maximum)){
                throw std::runtime_error("minimum or maximum for z axis is NaN");
            }
                
            // Get rows and columns
            int columns, rows;
            columns = d_vlen;
            rows = d_buffer->size() / d_vlen;

            // Fill data in spectrogram
            d_data->setValueMatrix(d_plot_data, d_axis_x, d_axis_y);
            d_data->setResampleMode(RangeAngleRasterData::BilinearInterpolation); 

            QwtInterval intervall = QwtInterval(d_axis_x.front(), d_axis_x.back());
            
            d_data->setInterval(Qt::XAxis, QwtInterval(d_axis_x.front(), d_axis_x.back()));
            d_data->setInterval(Qt::YAxis, QwtInterval(d_axis_y.front(), d_axis_y.back()));

            if (d_db_scale)
            {
                min_display = maximum-d_dynamic_range_db;
            }
            else
            {
                min_display = maximum/std::pow(10,d_dynamic_range_db/10);
            }

            if (d_autoscale_z) 
            {
                d_data->setInterval(Qt::ZAxis, QwtInterval(minimum, maximum));
            } 
            else 
            {
                d_data->setInterval(Qt::ZAxis, QwtInterval(min_display, maximum));
            }
            d_spectrogram->setData(d_data);

            QFont colorbar_font;
            colorbar_font.setPointSize(10); 

            // Set colorbar
            d_scale = d_plot->axisWidget(QwtPlot::yRight);
            d_scale->setColorBarEnabled(true);
            d_scale->setColorBarWidth(20);
            d_scale->setFont(colorbar_font);

            if (d_autoscale_z) 
            {
                d_scale->setColorMap(QwtInterval(minimum, maximum), d_colormap);
                d_plot->setAxisScale(QwtPlot::yRight, minimum, maximum);
            } 
            else 
            {
                d_scale->setColorMap(QwtInterval(min_display, maximum), d_colormap);
                d_plot->setAxisScale(QwtPlot::yRight, min_display, maximum);
            }
            d_plot->enableAxis(QwtPlot::yRight);
        
            // Do replot
            d_plot->replot();
        }
    }
} // namespace mimo_ofdm_jrc
} // namespace gr

