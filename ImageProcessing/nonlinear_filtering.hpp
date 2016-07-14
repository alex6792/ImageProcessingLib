/*!
 * \author Alexandre Krebs
 * \file nonlinear_filtering.hpp
 * \brief non linear filters
 */


#pragma once

#ifndef NONLINEAR_FILTERING_HPP
#define NONLINEAR_FILTERING_HPP


#include "../matrix.hpp"


Matrix<unsigned char> bilateral(Matrix<unsigned char>, int = 3, float = 1.0f, float = 1.0f);
Matrix<unsigned char> despeckle(Matrix<unsigned char>, int = 3);
Matrix<unsigned char> nagao(Matrix<unsigned char>, int = 3);


#endif // NONLINEAR_FILTERING_HPP
