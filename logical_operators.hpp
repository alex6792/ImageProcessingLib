/*!
 * \author Alexandre Krebs
 * \file logical_operators.hpp
 * \brief logical operators for matrices
 */


#pragma once

#ifndef lOGICAL_OPERATORS_HPP
#define lOGICAL_OPERATORS_HPP


#include "matrix.hpp"


// all, any and none functions
bool all(const Matrix<bool>&);
bool any(const Matrix<bool>&);
bool none(const Matrix<bool>&);

// logical operators
bool NOT(bool);
bool AND(bool, bool);
bool CDT(bool, bool);
bool EQ(bool, bool);
bool NAND(bool, bool);
bool NOR(bool, bool);
bool OR(bool, bool);
bool XOR(bool, bool);
bool XNOR(bool, bool);

// logical operators on matrix
Matrix<bool> NOT(const Matrix<bool>&);
Matrix<bool> AND(const Matrix<bool>&, const Matrix<bool>&);
Matrix<bool> CDT(const Matrix<bool>&, const Matrix<bool>&);
Matrix<bool> EQ(const Matrix<bool>&, const Matrix<bool>&);
Matrix<bool> NAND(const Matrix<bool>&, const Matrix<bool>&);
Matrix<bool> NOR(const Matrix<bool>&, const Matrix<bool>&);
Matrix<bool> OR(const Matrix<bool>&, const Matrix<bool>&);
Matrix<bool> XOR(const Matrix<bool>&, const Matrix<bool>&);
Matrix<bool> XNOR(const Matrix<bool>&, const Matrix<bool>&);


#endif // lOGICAL_OPERATORS_HPP
