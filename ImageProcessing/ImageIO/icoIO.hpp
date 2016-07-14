/*!
 * \author Alexandre Krebs
 * \file icoIO.hpp
 * \brief ICO images reader/writer
 */


#pragma once

#ifndef ICOIO_HPP
#define ICOIO_HPP


#include "../../matrix.hpp"
#include "../color.hpp"


Matrix<Color> read_ico(std::string);
void save_ico(std::string, const Matrix<Color>&);


#endif // ICOIO_HPP
