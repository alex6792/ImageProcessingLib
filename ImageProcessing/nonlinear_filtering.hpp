/*!
 * \author Alexandre Krebs
 * \file nonlinear_filtering.hpp
 * \brief non linear filters
 */


#pragma once

#ifndef NONLINEAR_FILTERING_HPP
#define NONLINEAR_FILTERING_HPP


#include "../matrix.hpp"


Matrix<float> variance_filter(const Matrix<float>&, std::size_t = 3);
Matrix<float> geometric_mean_filter(const Matrix<float>&, std::size_t = 3);
Matrix<float> harmonic_mean_filter(const Matrix<float>&, std::size_t = 3);
Matrix<float> quadratic_mean_filter(const Matrix<float>& , std::size_t = 3);

Matrix<float> imguided_filtering(const Matrix<float>&, const Matrix<float>& , std::size_t = 3);

Matrix<unsigned char> bilateral(const Matrix<unsigned char>&, std::size_t = 3, float = 1.0f, float = 1.0f);
Matrix<unsigned char> despeckle(const Matrix<unsigned char>&, std::size_t = 3);
Matrix<unsigned char> nagao(const Matrix<unsigned char>&, std::size_t = 3);


#endif // NONLINEAR_FILTERING_HPP
