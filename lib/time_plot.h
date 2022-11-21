#include <QApplication>
#include <QTimer>

#include <qwt_plot.h>
#include <qwt_plot_curve.h>
#include <qwt_plot_grid.h>
#include <qwt_plot_marker.h>
#include <qwt_symbol.h>

namespace gr {
  namespace mimo_ofdm_jrc {

    class time_plot : public QWidget
{
    Q_OBJECT

public:
    time_plot(int interval,
                std::string y_tag,
                std::string y_label,
                std::vector<float> axis_y,
                float time_window,
                std::vector<float>* y,
                bool* y_read,
                std::string title_label,
                QWidget* parent = 0);
    ~time_plot();

private:
    int d_interval;
    std::string d_y_label, d_y_tag;
    std::vector<float> d_axis_y;
    std::vector<float>* d_y;
    bool* d_y_read;
    float d_time_window;

    QwtPlot* d_plot;
    QwtSymbol* d_symbol;
    QwtPlotGrid* d_grid;
    std::vector<std::vector<QwtPlotMarker*>> d_marker;

    QTimer* d_timer;

    int d_refresh_counter;

protected:
    void resizeEvent(QResizeEvent* event);

public slots:
    void refresh();



};

}
}
