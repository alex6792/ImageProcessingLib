/*!
 * \author Alexandre Krebs
 * \file imgconverter.hpp
 * \brief binary, grayscale, colored image conversion
 */


#pragma once

#ifndef IMGCONVERTER_HPP
#define IMGCONVERTER_HPP


#include "color.hpp"
#include "../matrix.hpp"


Matrix<bool> gray2bwimage(const Matrix<unsigned char>&, unsigned char = 128);
Matrix<bool> color2bwimage(const Matrix<Color>&, unsigned char = 128);
Matrix<unsigned char> bw2grayimage(const Matrix<bool>&);
Matrix<unsigned char> color2grayimage(const Matrix<Color>&);
Matrix<Color> bw2colorimage(const Matrix<bool>&);
Matrix<Color> gray2colorimage(const Matrix<unsigned char>&);

Matrix<Color> array2colorimage(const Matrix<std::size_t>&);

Matrix<std::size_t> histogram(const Matrix<bool>&);
Matrix<std::size_t> histogram(const Matrix<unsigned char>&);
Matrix<std::size_t> histogram(const Matrix<std::size_t>&);


#endif // IMGCONVERTER_HPP
