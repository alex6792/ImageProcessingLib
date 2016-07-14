/*!
 * \author Alexandre Krebs
 * \file morphology.hpp
 * \brief morphology tools
 */


#pragma once

#ifndef MORPHOLOGY_HPP
#define MORPHOLOGY_HPP


#include "../matrix.hpp"


typedef Matrix<bool> Mask;

Mask circle(int);
Mask cross(int);
Mask ellipse(int, int);
Mask diamond(int);
Mask rect(int, int);
Mask square(int);


Matrix<bool> BTH(const Matrix<bool>&, const Mask& = square(3));
Matrix<bool> close(const Matrix<bool>&, const Mask& = square(3));
Matrix<bool> closing_by_reconstruction(const Matrix<bool>&, const Mask& = square(3));
Matrix<bool> conservative_smoothing(const Matrix<bool>&);
Matrix<bool> erode(const Matrix<bool>&, const Mask& = square(3));
Matrix<bool> dilate(const Matrix<bool>&, const Mask& = square(3));
Matrix<bool> gradient(const Matrix<bool>&, const Mask& = square(3));
Matrix<bool> HitOrMiss(const Matrix<bool>&, const Mask&, const Mask&);
Matrix<bool> inner_gradient(const Matrix<bool>&, const Mask& = square(3));
Matrix<bool> median_filter(const Matrix<bool>&, const Mask& = square(3));
Matrix<bool> open(const Matrix<bool>&, const Mask& = square(3));
Matrix<bool> opening_by_reconstruction(const Matrix<bool>&, const Mask& = square(3));
Matrix<bool> reconstruct(const Matrix<bool>&, const Matrix<bool>& , const Mask& = square(3));
Matrix<bool> outer_gradient(const Matrix<bool>&, const Mask& = square(3));
Matrix<bool> skeleton(const Matrix<bool>&, const Mask& = square(3));
Matrix<bool> WTH(const Matrix<bool>&, const Mask& = square(3));

Matrix<int> label(const Matrix<bool>&, const Mask& = square(3));


Matrix<unsigned char> BTH(const Matrix<unsigned char>&, const Mask& = square(3));
Matrix<unsigned char> close(const Matrix<unsigned char>&, const Mask& = square(3));
Matrix<unsigned char> closing_by_reconstruction(const Matrix<unsigned char>&, const Mask& = square(3));
Matrix<unsigned char> conservative_smoothing(const Matrix<unsigned char>&);
Matrix<unsigned char> erode(const Matrix<unsigned char>&, const Mask& = square(3));
Matrix<unsigned char> dilate(const Matrix<unsigned char>&, const Mask& = square(3));
Matrix<unsigned char> gradient(const Matrix<unsigned char>&, const Mask& = square(3));
Matrix<unsigned char> HitOrMiss(const Matrix<unsigned char>&, const Mask&, const Mask&);
Matrix<unsigned char> inner_gradient(const Matrix<unsigned char>&, const Mask& = square(3));
Matrix<unsigned char> median_filter(const Matrix<unsigned char>&, const Mask& = square(3));
Matrix<unsigned char> open(const Matrix<unsigned char>&, const Mask& = square(3));
Matrix<unsigned char> opening_by_reconstruction(const Matrix<unsigned char>&, const Mask& = square(3));
Matrix<unsigned char> outer_gradient(const Matrix<unsigned char>&, const Mask& = square(3));
Matrix<unsigned char> reconstruct(const Matrix<unsigned char>&, const Matrix<unsigned char>& , const Mask& = square(3));
Matrix<unsigned char> skeleton(const Matrix<unsigned char>&, const Mask& = square(3));
Matrix<unsigned char> WTH(const Matrix<unsigned char>&, const Mask& = square(3));

Matrix<int> label(const Matrix<unsigned char>&, const Mask& = square(3));


#endif // MORPHOLOGY_HPP
