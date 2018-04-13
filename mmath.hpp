/*!
 * \author Alexandre Krebs
 * \file mmath.hpp
 * \brief mathematical functions for matrices
 */


#pragma once

#ifndef MMATH_HPP
#define MMATH_HPP


#include "matrix.hpp"


#define PI 4*std::atan(1.0)


// mathematical functions
template <class T> Matrix<T> abs(const Matrix<T>&);
template <class T> Matrix<T> max(const Matrix<T>&, const Matrix<T>&);
template <class T> Matrix<T> min(const Matrix<T>&, const Matrix<T>&);

// trigonometric functions
template <class T> Matrix<T> cos(const Matrix<T>&);
template <class T> Matrix<T> sin(const Matrix<T>&);
template <class T> Matrix<T> tan(const Matrix<T>&);
template <class T> Matrix<T> acos(const Matrix<T>&);
template <class T> Matrix<T> asin(const Matrix<T>&);
template <class T> Matrix<T> atan(const Matrix<T>&);
template <class T> Matrix<T> atan2(const Matrix<T>&, const Matrix<T>&);

// hyperbolic functions
template <class T> Matrix<T> cosh(const Matrix<T>&);
template <class T> Matrix<T> sinh(const Matrix<T>&);
template <class T> Matrix<T> tanh(const Matrix<T>&);
template <class T> Matrix<T> acosh(const Matrix<T>&);
template <class T> Matrix<T> asinh(const Matrix<T>&);
template <class T> Matrix<T> atanh(const Matrix<T>&);

// exponential and logarithmic functions
template <class T> Matrix<T> exp(const Matrix<T>&);
template <class T> Matrix<T> exp2(const Matrix<T>&);
template <class T> Matrix<T> log(const Matrix<T>&);
template <class T> Matrix<T> log2(const Matrix<T>&);
template <class T> Matrix<T> log10(const Matrix<T>&);

// power functions
template <class T> Matrix<T> pow(const Matrix<T>&, const T&);
template <class T> Matrix<T> sqrt(const Matrix<T>&);
template <class T> Matrix<T> cbrt(const Matrix<T>&);
template <class T> Matrix<T> hypot(const Matrix<T>&, const Matrix<T>&);

// gamma functions
template <class T> Matrix<T> tgamma(const Matrix<T>&);
template <class T> Matrix<T> lgamma(const Matrix<T>&);

// rounding functions
template <class T> Matrix<T> ceil(const Matrix<T>&);
template <class T> Matrix<T> floor(const Matrix<T>&);
template <class T> Matrix<T> trunc(const Matrix<T>&);
template <class T> Matrix<T> round(const Matrix<T>&);

// classification functions
template <class T> Matrix<bool> isfinite(const Matrix<T>&);
template <class T> Matrix<bool> isinf(const Matrix<T>&);
template <class T> Matrix<bool> isnan(const Matrix<T>&);


#include "mmath.tcc"


#endif // MMATH_HPP
