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
Matrix<float> disk(std::size_t = 2);
Matrix<float> gaussian(std::size_t = 5, float = 1.0f);
Matrix<float> isotropicx();
Matrix<float> isotropicy();
Matrix<float> kirch();
Matrix<float> laplacian(float);
Matrix<float> log(int = 5, float = 1.0f);
Matrix<float> prewittx();
Matrix<float> prewitty();
Matrix<float> robertsx();
Matrix<float> robertsy();
Matrix<float> sobelx();
Matrix<float> sobely();

Matrix<unsigned char> gradient(Matrix<unsigned char>, std::string);
Matrix<unsigned char> filter(Matrix<unsigned char>, Matrix<float>, int = 0);
Matrix<float> conv(Matrix<float>, Matrix<float>, int = 0);


#endif // LINEAR_FILTERING_HPP
