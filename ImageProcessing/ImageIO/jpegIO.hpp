/*!
 * \author Alexandre Krebs
 * \file jpegIO.hpp
 * \brief JPEG images reader/writer
 */


#pragma once

#ifndef JPEGIO_HPP
#define JPEGIO_HPP


#include "../../matrix.hpp"
#include "../color.hpp"


Matrix<Color> read_jpeg(std::string);
void save_jpeg(std::string, const Matrix<Color>&);


#endif // JPEGIO_HPP
