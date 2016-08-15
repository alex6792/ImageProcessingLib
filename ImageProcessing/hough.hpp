/*!
 * \author Alexandre Krebs
 * \file hough.hpp
 * \brief hough transform
 */


#pragma once

#ifndef HOUGH_HPP
#define HOUGH_HPP


#include "../matrix.hpp"


Matrix<std::size_t> hough_transform(const Matrix<bool>&);


#endif // HOUGH_HPP
