/*!
 * \author Alexandre Krebs
 * \file fourier.hpp
 * \brief FFT implementation
 */


#pragma once

#ifndef FOURIER_HPP
#define FOURIER_HPP


#include <complex>
#include "../matrix.hpp"


Matrix<std::complex<float> > FFT(const Matrix<float>&);
Matrix<std::complex<float> > iFFT(const Matrix<std::complex<float> >&);


#endif // FOURIER_HPP
