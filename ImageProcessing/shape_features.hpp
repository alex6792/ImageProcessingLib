/*!
 * \author Alexandre Krebs
 * \file shape_features.hpp
 * \brief shape features
 */


#pragma once

#ifndef SHAPE_FEATURES_HPP
#define SHAPE_FEATURES_HPP


#include "../matrix.hpp"

float circularity(const Matrix<bool>&);
float major_axis(const Matrix<bool> &);


#endif // SHAPE_FEATURES_HPP
