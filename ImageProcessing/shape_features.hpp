/*!
 * \author Alexandre Krebs
 * \file shape_features.hpp
 * \brief shape features
 */


#pragma once

#ifndef SHAPE_FEATURES_HPP
#define SHAPE_FEATURES_HPP


#include "../matrix.hpp"
#include <complex>

float circularity(const Matrix<bool>&);
float major_axis(const Matrix<bool> &);

class Zernike
{
    private :
        int m, n;
        float (*rho_fct)(float);

    public :
        Zernike(int, int);
        std::complex<float> polynom(float, float);
};


#endif // SHAPE_FEATURES_HPP
