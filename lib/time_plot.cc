#include "time_plot.h"
#include <iostream>

namespace gr {
namespace mimo_ofdm_jrc {

time_plot::time_plot(int interval,
                        std::string y_tag,
                        std::string y_label,
                        std::vector<float> axis_y,
                        float time_window,
                        std::vector<float>* y,
                        bool* y_read,
                        std::string title_label,
                        QWidget* parent)
    : QWidget(parent)
{
    d_interval = interval;
    d_axis_y = axis_y;
    d_y = y;
    d_y_read = y_read;
    d_y_label = y_label;
    d_y_tag = y_tag;
    d_time_window = time_window;
    d_refresh_counter = 0;
    d_marker.resize((int)(time_window / interval * 1000));

    // Setup GUI
    resize(QSize(960, 470));
    
    d_plot = new QwtPlot(this);
    d_symbol = new QwtSymbol(QwtSymbol::Ellipse, Qt::blue, Qt::NoPen, QSize(12, 12));
    d_grid = new QwtPlotGrid;

    // setWindowTitle(label_title.c_str());
    setWindowTitle("Time Plot");

    std::string label_title = " vs. Time";

    if (title_label != "") 
    {
        label_title = title_label;
        // label_title.append(" (");
        // label_title.append(title_label);
        // label_title.append(")");
    }
    else
    {
        label_title.insert(0, d_y_tag);
    }
        
    QwtText title_text;
    QFont title_font;
    title_font.setPointSize(16); 
    title_font.setBold(true);
    title_text.setFont(title_font);
    title_text.setText(label_title.c_str());
    d_plot->setTitle(title_text);

    QFont axis_tick_font;
    axis_tick_font.setPointSize(15); 

    QwtText axis_label_text;
    QFont axis_label_font;
    axis_label_font.setPointSize(16); 
    axis_label_font.setBold(true);
    axis_label_text.setFont(axis_label_font);
    axis_label_text.setText("Time [s]");

    d_plot->setAxisFont(QwtPlot::xBottom, axis_tick_font);
    d_plot->setAxisScale(QwtPlot::xBottom, 0, d_time_window);
    d_plot->setAxisTitle(QwtPlot::xBottom, axis_label_text);
    
    axis_label_text.setText(d_y_label.c_str());
    d_plot->setAxisFont(QwtPlot::yLeft, axis_tick_font);
    d_plot->setAxisScale(QwtPlot::yLeft, d_axis_y[0], d_axis_y[1]);
    d_plot->setAxisTitle(QwtPlot::yLeft, axis_label_text);

    // Add grid
    d_grid->setPen(QPen(QColor(119, 136, 153), 0.5, Qt::DashLine));
    d_grid->attach(d_plot);

    d_plot->setCanvasBackground(Qt::white); // nice blue

    // Plot grid and axis
    d_plot->replot();

    // Setup timer and connect refreshing plot
    d_timer = new QTimer(this);
    connect(d_timer, SIGNAL(timeout()), this, SLOT(refresh()));
    d_timer->start(d_interval);
}

time_plot::~time_plot() {}

void time_plot::resizeEvent(QResizeEvent* event)
{
    d_plot->setGeometry(0, 0, this->width(), this->height());
}

void time_plot::refresh()
{
    // Reset axis
    if (d_time_window <= d_refresh_counter * float(d_interval) / 1000.0) {
        d_plot->setAxisScale(QwtPlot::xBottom,
                             d_refresh_counter * float(d_interval) / 1000.0 -
                                 d_time_window,
                             d_refresh_counter * float(d_interval) / 1000.0);
    }

    // Detach old markers
    int marker_index = d_refresh_counter % d_marker.size();
    for (int k = 0; k < d_marker[marker_index].size(); k++) 
    {
        d_marker[marker_index][k]->detach();
    }

    if (!(*d_y_read)) 
    {
        // Set new marker
        for (int k = 0; k < d_y->size(); k++) 
        {
            if (k < d_marker[marker_index].size()) 
            {
                d_marker[marker_index][k]->setValue( QPointF(d_refresh_counter * float(d_interval) / 1000.0, (*d_y)[k]) );
                d_marker[marker_index][k]->attach(d_plot);
            } 
            else 
            {
                d_marker[marker_index].push_back(new QwtPlotMarker);
                d_marker[marker_index][k]->setSymbol(d_symbol);
                d_marker[marker_index][k]->setValue( QPointF(d_refresh_counter * float(d_interval) / 1000.0, (*d_y)[k]) );
                d_marker[marker_index][k]->attach(d_plot);
            }
        }
        *d_y_read = true; // set points as read
    }
    // Replot and increment counter
    d_plot->replot();
    d_refresh_counter++;
}

} // namespace radar
} // namespace gr
