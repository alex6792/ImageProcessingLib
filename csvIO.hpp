/*!
 * \author Alexandre Krebs
 * \file csvIO.hpp
 * \brief CSV files reader/writer
 */


#pragma once

#ifndef CSVIO_HPP
#define CSVIO_HPP


#include "matrix.hpp"


template <class T> Matrix<T> read_csv(std::string, char = ',');
template <class T> void save_csv(std::string, const Matrix<T>&, char = ',');


#include "csvIO.tcc"


#endif // CSVIO_HPP
