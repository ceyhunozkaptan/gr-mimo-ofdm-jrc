#include "qwt_global.h"
#include "qwt_raster_data.h"
#include <qvector.h>

class RangeAngleRasterData : public QwtRasterData
{
  public:
    /*!
       \brief Resampling algorithm
       The default setting is NearestNeighbour;
     */
    enum ResampleMode
    {
        /*!
           Return the value from the matrix, that is nearest to the
           the requested position.
         */
        NearestNeighbour,

        /*!
           Interpolate the value from the distances and values of the
           4 surrounding values in the matrix,
         */
        BilinearInterpolation,
    };

    RangeAngleRasterData();
    virtual ~RangeAngleRasterData();

    void setResampleMode(ResampleMode mode);
    ResampleMode resampleMode() const;

    void setInterval( Qt::Axis, const QwtInterval & );
    // virtual QwtInterval interval( Qt::Axis axis) const;

    void setValueMatrix( const QVector< double > &values, const QVector< double > &x_axis, const QVector< double > &y_axis);
    const QVector< double > valueMatrix() const;

    void setValue( int row, int col, double value );

    int numColumns() const;
    int numRows() const;

    virtual QRectF pixelHint( const QRectF & ) const;

    virtual double value( double x, double y ) const;

    long indexOfClosest(const QVector<double>& sorted_array, double x) const;

  private:
    // void update();
    
    class PrivateData;
    PrivateData *m_data;
};