/*!
 * \author Alexandre Krebs
 * \file linear_filtering.hpp
 * \brief linear filters
 */


#pragma once

#ifndef LINEAR_FILTERING_HPP
#define LINEAR_FILTERING_HPP


#include "../matrix.hpp"


Matrix<float> average(std::size_t = 3);
Matrix<float> binomial(std::size_t = 3);
Matrix<float> disk(std::size_t = 2);
Matrix<float> gaussian(std::size_t = 7, float = 1.0f);
Matrix<float> isotropicx();
Matrix<float> isotropicy();
Matrix<float> kirch();
Matrix<float> laplacian(float);
Matrix<float> log(std::size_t = 5, float = 1.0f);
Matrix<float> mdifx(std::size_t = 7, float = 1.0f);
Matrix<float> mdify(std::size_t = 7, float = 1.0f);
Matrix<float> mdifxx(std::size_t = 7, float = 1.0f);
Matrix<float> mdifxy(std::size_t = 7, float = 1.0f);
Matrix<float> mdifyy(std::size_t = 7, float = 1.0f);
Matrix<float> prewittx();
Matrix<float> prewitty();
Matrix<float> pyramidal(std::size_t = 5);
Matrix<float> robertsx();
Matrix<float> robertsy();
Matrix<float> robinson();
Matrix<float> savgol(std::size_t = 5, std::size_t = 2, std::size_t = 0, std::size_t = 0);
Matrix<float> scharrx();
Matrix<float> scharry();
Matrix<float> sobelx();
Matrix<float> sobely();

Matrix<unsigned char> gradient(const Matrix<unsigned char>&, std::string = "sobel");
Matrix<unsigned char> filter(const Matrix<unsigned char>&, const Matrix<float>&, int = 0);
Matrix<float> conv(const Matrix<float>&, const Matrix<float>&, std::string = "same");
Matrix<float> xcorr(const Matrix<float>&, const Matrix<float>&, std::string = "same");
Matrix<float> interp1(const Matrix<float>&, const Matrix<float>&, const Matrix<float>&, std::string = "linear");

#endif // LINEAR_FILTERING_HPP
