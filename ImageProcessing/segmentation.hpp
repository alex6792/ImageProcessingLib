/*!
 * \author Alexandre Krebs
 * \file segmentation.hpp
 * \brief Segmentation Algorithms
 */


#pragma once

#ifndef SEGMENTATION_HPP
#define SEGMENTATION_HPP


#include "../matrix.hpp"


Matrix<std::size_t> watershed(const Matrix<unsigned char>&, const Matrix<std::size_t>&);


#endif // SEGMENTATION_HPP
