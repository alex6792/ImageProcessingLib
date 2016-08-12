/*!
 * \author Alexandre Krebs
 * \file imageIO.hpp
 * \brief images reader/writer
 */


#pragma once

#ifndef IMAGEIO_HPP
#define IMAGEIO_HPP


#include "../../matrix.hpp"
#include "../color.hpp"


Matrix<Color> read_img(std::string);
void save_img(std::string, const Matrix<Color>&);


#endif // IMAGEIO_HPP
