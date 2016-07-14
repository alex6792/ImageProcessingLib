/*!
 * \author Alexandre Krebs
 * \file gifIO.hpp
 * \brief GIF images reader/writer
 */


#pragma once

#ifndef GIFIO_HPP
#define GIFIO_HPP


#include "../../matrix.hpp"
#include "../color.hpp"


Matrix<Color> read_gif(std::string);
void save_gif(std::string, const Matrix<Color>&);


#endif // GIFIO_HPP
