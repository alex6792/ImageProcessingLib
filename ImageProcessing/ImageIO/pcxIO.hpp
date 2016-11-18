/*!
 * \author Alexandre Krebs
 * \file pcxIO.hpp
 * \brief PCX images reader/writer
 */


#pragma once

#ifndef PCXIO_HPP
#define PCXIO_HPP


#include "../../matrix.hpp"
#include "../color.hpp"


Matrix<Color> read_pcx(std::string);
void save_pcx(std::string, const Matrix<Color>&);


#endif // PCXIO_HPP
