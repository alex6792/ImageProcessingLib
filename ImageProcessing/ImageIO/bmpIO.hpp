/*!
 * \author Alexandre Krebs
 * \file bmpIO.hpp
 * \brief BMP images reader/writer
 */


#pragma once

#ifndef BMPIO_HPP
#define BMPIO_HPP


#include "../../matrix.hpp"
#include "../color.hpp"


Matrix<Color> read_bmp(std::string);
void save_bmp(std::string, const Matrix<Color>&);


#endif // BMPIO_HPP
