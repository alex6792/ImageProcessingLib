/*!
 * \author Alexandre Krebs
 * \file tgaIO.hpp
 * \brief TGA images reader/writer
 */


#pragma once

#ifndef TGAIO_HPP
#define TGAIO_HPP


#include "../../matrix.hpp"
#include "../color.hpp"


Matrix<Color> read_tga(std::string);
void save_tga(std::string, const Matrix<Color>&);


#endif // TGAIO_HPP
