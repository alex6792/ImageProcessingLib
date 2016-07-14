/*!
 * \author Alexandre Krebs
 * \file tiffIO.hpp
 * \brief TIFF images reader/writer
 */


#pragma once

#ifndef TIFFIO_HPP
#define TIFFIO_HPP


#include "../../matrix.hpp"
#include "../color.hpp"


Matrix<Color> read_tiff(std::string);
void save_tiff(std::string, const Matrix<Color>&);


#endif // TIFFIO_HPP
