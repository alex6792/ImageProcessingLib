/*!
 * \author Alexandre Krebs
 * \file regression.hpp
 * \brief polynomial regression
 */


#pragma once

#ifndef REGRESSION_HPP
#define REGRESSION_HPP


#include "matrix.hpp"


template <class T> Matrix<T> poly_regression(const Matrix<T>&, const Matrix<T>&, int);


#include "regression.tcc"


#endif // REGRESSION_HPP
