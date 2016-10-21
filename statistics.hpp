/*!
 * \author Alexandre Krebs
 * \file statistics.hpp
 * \brief statistics tools for matrices
 */


#pragma once

#ifndef STATISTICS_HPP
#define STATISTICS_HPP


#include "matrix.hpp"


template <class T> Matrix<std::size_t> argmax(const Matrix<T>&);
template <class T> Matrix<std::size_t> argmin(const Matrix<T>&);
template <class T> Matrix<std::size_t> argsort(const Matrix<T>&);
template <class T> Matrix<T> axismax(const Matrix<T>&, int);
template <class T> Matrix<T> axismean(const Matrix<T>&, int);
template <class T> Matrix<T> axismin(const Matrix<T>&, int);
template <class T> Matrix<T> axisprod(const Matrix<T>&, int);
template <class T> Matrix<T> axisstdev(const Matrix<T>&, int);
template <class T> Matrix<T> axissum(const Matrix<T>&, int);
template <class T> Matrix<T> axisvar(const Matrix<T>&, int);
template <class T> T central_moment(const Matrix<T>&, int);
template <class T> T corrected_var(const Matrix<T>&);
template <class T> Matrix<T> cumsum(const Matrix<T>&);
template <class T> T geometric_mean(const Matrix<T>&);
template <class T> T max(const Matrix<T>&);
template <class T> T mean(const Matrix<T>&);
template <class T> T median(const Matrix<T>&);
template <class T> T min(const Matrix<T>&);
template <class T> T moment(const Matrix<T>&, int);
template <class T> T moment(const Matrix<T>&, int, int);
template <class T> T prod(const Matrix<T>&);
template <class T> T stdev(const Matrix<T>&);
template <class T> int quickselect(const Matrix<T>&, int);
template <class T> Matrix<T> quicksort(const Matrix<T>&);
template <class T> Matrix<T> sort(const Matrix<T>&);
template <class T> T sum(const Matrix<T>&);
template <class T> T var(const Matrix<T>&);
template <class T> static int partition(Matrix<T>&, int, int, int);
template <class T> static int quickselect(Matrix<T>&, int, int, int);
template <class T> static void quicksort(Matrix<T>&, int, int);


#include "statistics.tcc"


#endif // STATISTICS_HPP
