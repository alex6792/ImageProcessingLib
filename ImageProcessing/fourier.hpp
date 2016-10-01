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
Matrix<std::complex<float> > FFTshift(const Matrix<std::complex<float> >&);
Matrix<float> FFTfreq(std::size_t, float = 1.0f);


#endif // FOURIER_HPP
