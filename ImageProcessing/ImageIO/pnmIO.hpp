/*!
 * \author Alexandre Krebs
 * \file pnmIO.hpp
 * \brief Portable Pixmap images reader/writer
 */


#pragma once

#ifndef PNMIO_HPP
#define PNMIO_HPP


#include "../../matrix.hpp"
#include "../color.hpp"


Matrix<Color> read_pnm(std::string);
Matrix<bool> read_pbm(std::string);
Matrix<unsigned char> read_pgm(std::string);
Matrix<Color> read_ppm(std::string);

void save_pnm(std::string, const Matrix<Color>&);
void save_pbm(std::string, const Matrix<bool>&);
void save_pgm(std::string, const Matrix<unsigned char>&);
void save_ppm(std::string, const Matrix<Color>&);


#endif // PNMIO_HPP
