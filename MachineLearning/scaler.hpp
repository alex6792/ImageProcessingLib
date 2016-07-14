/*!
 * \author Alexandre Krebs
 * \file scaler.hpp
 * \brief Tools to scale data
 */


#pragma once

#ifndef SCALER_HPP
#define SCALER_HPP


#include "../matrix.hpp"


/*!
 * \class MinMaxScaler
 * \brief scale data according the minimum and maximum
 */
class MinMaxScaler
{
    private :
        Matrix<float> data_min;// min matrix
        Matrix<float> data_max;// max matrix
        Matrix<float> data_range, data_range_inv;
        Matrix<float> min;
        float range_min, range_max;
        int axis;

    public :
        MinMaxScaler(float range_min = 0.0f, float range_max = 1.0f, int axis = 0);
        void fit(const Matrix<float>&);
        Matrix<float> fit_transform(const Matrix<float>&);
        Matrix<float> inverse_transform(const Matrix<float>&);
        Matrix<float> transform(const Matrix<float>&);
};


/*!
 * \class StandardScaler
 * \brief scale data by removing the mean and divide by the std
 */
class StandardScaler
{
    private :
        Matrix<float> data_mean;// mean matrix
        Matrix<float> data_std;// std matrix
        Matrix<float> data_std_inv;
        int axis;
        bool with_mean, with_std;

    public :
        StandardScaler(bool with_mean = true, bool with_std = true, int axis = 0);
        void fit(const Matrix<float>&);
        Matrix<float> fit_transform(const Matrix<float>&);
        Matrix<float> inverse_transform(const Matrix<float>&);
        Matrix<float> transform(const Matrix<float>&);
};


#endif // SCALER_HPP
