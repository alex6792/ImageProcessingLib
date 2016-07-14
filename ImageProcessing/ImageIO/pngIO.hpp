/*!
 * \author Alexandre Krebs
 * \file pngIO.hpp
 * \brief PNG images reader/writer
 */


#pragma once

#ifndef PNGIO_HPP
#define PNGIO_HPP


#include "../../matrix.hpp"
#include "../color.hpp"


Matrix<Color> read_png(std::string);
void save_png(std::string, const Matrix<Color>&);


#endif // PNGIO_HPP
