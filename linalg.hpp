/*!
 * \author Alexandre Krebs
 * \file linalg.hpp
 * \brief linear algebra tools
 */


#pragma once

#ifndef LINALG_HPP
#define LINALG_HPP


#include "matrix.hpp"


template <class T> Matrix<T> dot(const Matrix<T>&, const Matrix<T>&);
template <class T> Matrix<T> inv(const Matrix<T>&);
template <class T> T det(const Matrix<T>&);
template <class T> T norm(const Matrix<T>&);
template <class T> Matrix<T> bwdsub(const Matrix<T>&, const Matrix<T>&);
template <class T> Matrix<T> fwdsub(const Matrix<T>&, const Matrix<T>&);
template <class T> Matrix<T> cholesky(const Matrix<T>&);
template <class T> std::pair<Matrix<T>, Matrix<T> > crout(const Matrix<T>&);
template <class T> std::tuple<Matrix<T>, Matrix<T>, Matrix<T> > lu(const Matrix<T>&);
template <class T> std::pair<Matrix<T>, Matrix<T> > qr(const Matrix<T>&);
template <class T> std::pair<Matrix<T>, Matrix<T> > svd(const Matrix<T>&);
template <class T> T trace(const Matrix<T>&);
template <class T> std::pair<Matrix<T>, Matrix<T> > Jacobi(const Matrix<T>&);


#include "linalg.tcc"


#endif // LINALG_HPP
