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

#include "range_angle_raster_data.h"
// #include "qwt_interval.h"
#include <qnumeric.h>
#include <qmath.h>
#include <iostream>

class RangeAngleRasterData::PrivateData
{
  public:
    PrivateData()
        : resampleMode( RangeAngleRasterData::NearestNeighbour )
        , numColumns(0)
    {
    }

    inline double value(int row, int col) const
    {
        return values.data()[ row * numColumns + col ];
    }

    RangeAngleRasterData::ResampleMode resampleMode;

    QVector< double > values;

    QVector< double > x_axis; // columns
    QVector< double > y_axis; // row

    int numColumns;
    int numRows;
};

//! Constructor
RangeAngleRasterData::RangeAngleRasterData()
{
    m_data = new PrivateData();
    // update();
}

//! Destructor
RangeAngleRasterData::~RangeAngleRasterData()
{
    delete m_data;
}

/*!
   \brief Set the resampling algorithm
   \param mode Resampling mode
   \sa resampleMode(), value()
 */
void RangeAngleRasterData::setResampleMode( ResampleMode mode )
{
    m_data->resampleMode = mode;
}

/*!
   \return resampling algorithm
   \sa setResampleMode(), value()
 */
RangeAngleRasterData::ResampleMode RangeAngleRasterData::resampleMode() const
{
    return m_data->resampleMode;
}

/*!
   \brief Assign the bounding interval for an axis
   Setting the bounding intervals for the X/Y axis is mandatory
   to define the positions for the values of the value matrix.
   The interval in Z direction defines the possible range for
   the values in the matrix, what is f.e used by QwtPlotSpectrogram
   to map values to colors. The Z-interval might be the bounding
   interval of the values in the matrix, but usually it isn't.
   ( f.e a interval of 0.0-100.0 for values in percentage )
   \param axis X, Y or Z axis
   \param interval Interval
   \sa QwtRasterData::interval(), setValueMatrix()
 */
void RangeAngleRasterData::setInterval( Qt::Axis axis, const QwtInterval& interval )
{
    if ( axis >= 0 && axis <= 2 )
    {
        // m_data->intervals[axis] = interval;
        QwtRasterData::setInterval( axis, interval );
        // update();
    }
}

// /*!
//    \return Bounding interval for an axis
//    \sa setInterval
//  */
// QwtInterval RangeAngleRasterData::interval( Qt::Axis axis ) const
// {
//     if ( axis >= 0 && axis <= 2 )
//         return m_data->intervals[ axis ];

//     return QwtInterval();
// }

/*!
   \brief Assign a value matrix
   The positions of the values are calculated by dividing
   the bounding rectangle of the X/Y intervals into equidistant
   rectangles ( pixels ). Each value corresponds to the center of
   a pixel.
   \param values Vector of values
   \param numColumns Number of columns
   \sa valueMatrix(), numColumns(), numRows(), setInterval()()
 */
void RangeAngleRasterData::setValueMatrix(const QVector< double >& values, const QVector< double >& x_axis, const QVector< double >& y_axis)
{
    m_data->values = values;
    m_data->x_axis = x_axis;
    m_data->y_axis = y_axis;

    m_data->numColumns = m_data->x_axis.size();
    m_data->numRows = m_data->y_axis.size();

    // update();
}

/*!
   \return Value matrix
   \sa setValueMatrix(), numColumns(), numRows(), setInterval()
 */
const QVector< double > RangeAngleRasterData::valueMatrix() const
{
    return m_data->values;
}

/*!
   \brief Change a single value in the matrix
   \param row Row index
   \param col Column index
   \param value New value
   \sa value(), setValueMatrix()
 */
void RangeAngleRasterData::setValue( int row, int col, double value )
{
    if ( row >= 0 && row < m_data->numRows &&
        col >= 0 && col < m_data->numColumns )
    {
        const int index = row * m_data->numColumns + col;
        m_data->values.data()[ index ] = value;
    }
}

/*!
   \return Number of columns of the value matrix
   \sa valueMatrix(), numRows(), setValueMatrix()
 */
int RangeAngleRasterData::numColumns() const
{
    return m_data->numColumns;
}

/*!
   \return Number of rows of the value matrix
   \sa valueMatrix(), numColumns(), setValueMatrix()
 */
int RangeAngleRasterData::numRows() const
{
    return m_data->numRows;
}

/*!
   \brief Calculate the pixel hint
   pixelHint() returns the geometry of a pixel, that can be used
   to calculate the resolution and alignment of the plot item, that is
   representing the data.
   - NearestNeighbour\n
     pixelHint() returns the surrounding pixel of the top left value
     in the matrix.
   - BilinearInterpolation\n
     Returns an empty rectangle recommending
     to render in target device ( f.e. screen ) resolution.
   \param area Requested area, ignored
   \return Calculated hint
   \sa ResampleMode, setMatrix(), setInterval()
 */
QRectF RangeAngleRasterData::pixelHint( const QRectF& area ) const
{
    Q_UNUSED( area )

    QRectF rect;
    // if ( m_data->resampleMode == NearestNeighbour )
    // {
    //     const QwtInterval intervalX = interval( Qt::XAxis );
    //     const QwtInterval intervalY = interval( Qt::YAxis );
    //     if ( intervalX.isValid() && intervalY.isValid() )
    //     {
    //         rect = QRectF( intervalX.minValue(), intervalY.minValue(),
    //             m_data->dx, m_data->dy );
    //     }
    // }

    return rect;
}

long RangeAngleRasterData::indexOfClosest(const QVector<double>& sorted_array, double x) const
{

    auto iter_geq = std::lower_bound(
        sorted_array.begin(), 
        sorted_array.end(), 
        x
    );

    if (iter_geq == sorted_array.begin()) {
        return 0;
    }

    double a = *(iter_geq - 1);
    double b = *(iter_geq);

    if (std::abs(x - a) < std::abs(x - b)) {
        return iter_geq - sorted_array.begin() - 1;
    }

    return iter_geq - sorted_array.begin();

}

/*!
   \return the value at a raster position
   \param x X value in plot coordinates
   \param y Y value in plot coordinates
   \sa ResampleMode
 */
double RangeAngleRasterData::value( double x, double y ) const
{
    // std::cout << "[RANGE-ANGLE] x: " << x << ", y: " << y << std::endl;

    const QwtInterval xInterval = interval( Qt::XAxis );
    const QwtInterval yInterval = interval( Qt::YAxis );

    if ( !( xInterval.contains(x) && yInterval.contains(y) ) )
        return qQNaN();

    double value;

    auto i_col = indexOfClosest(m_data->x_axis, x);
    auto i_row = indexOfClosest(m_data->y_axis, y);

    switch( m_data->resampleMode )
    {
        case BilinearInterpolation:
        {
            // int col1 = qRound( ( x - xInterval.minValue() ) / m_data->dx ) - 1;
            // int row1 = qRound( ( y - yInterval.minValue() ) / m_data->dy ) - 1;

            int col1 = i_col;
            int row1 = i_row;

            int col2 = col1 + 1;
            int row2 = row1 + 1;

            if ( col1 < 0 )
                col1 = col2;
            else if ( col2 >= m_data->numColumns )
                col2 = col1;

            if ( row1 < 0 )
                row1 = row2;
            else if ( row2 >= m_data->numRows )
                row2 = row1;

            const double f11 = m_data->value( row1, col1 );
            const double f21 = m_data->value( row1, col2 );
            const double f12 = m_data->value( row2, col1 );
            const double f22 = m_data->value( row2, col2 );

            
            // const double x2 = xInterval.minValue() + ( col2 + 0.5 ) * m_data->dx;
            // const double y2 = yInterval.minValue() + ( row2 + 0.5 ) * m_data->dy;

            // const double rx2 = ( x2 - x ) / m_data->dx;
            // const double ry2 = ( y2 - y ) / m_data->dy;


            const double x1 = m_data->x_axis[col1];
            const double y1 = m_data->y_axis[row1];

            const double x2 = m_data->x_axis[col2];
            const double y2 = m_data->y_axis[row2];

            const double rx1 = ( x - x1 ) / (x2 - x1);
            const double rx2 = ( x2 - x ) / (x2 - x1);

            const double ry1 = ( y - y1 ) / (y2 - y1);
            const double ry2 = ( y2 - y ) / (y2 - y1);

            // const double vr1 = rx2 * f11 + ( 1.0 - rx2 ) * f21;
            // const double vr2 = rx2 * f12 + ( 1.0 - rx2 ) * f22;
            // value = ry2 * vr1 + ( 1.0 - ry2 ) * vr2;

            value = ry2 * (rx2 * f11 + rx1 * f21) + ry1 * ( rx2 * f12 + rx1 * f22 );

            break;
        }
        case NearestNeighbour:
        default:
        {
            // int row = int( ( y - yInterval.minValue() ) / m_data->dy );
            // int col = int( ( x - xInterval.minValue() ) / m_data->dx );

            int row = i_row;
            int col = i_col;

            // In case of intervals, where the maximum is included
            // we get out of bound for row/col, when the value for the
            // maximum is requested. Instead we return the value
            // from the last row/col

            if ( row >= m_data->numRows )
                row = m_data->numRows - 1;

            if ( col >= m_data->numColumns )
                col = m_data->numColumns - 1;

            value = m_data->value( row, col );
        }
    }

    return value;
}

