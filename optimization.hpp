/*!
 * \author Alexandre Krebs
 * \file optimization.hpp
 * \brief optimization tools
 */


#pragma once

#ifndef OPTIMIZATION_HPP
#define OPTIMIZATION_HPP


#include "matrix.hpp"

Matrix<float> steepest_descent(std::pair<float, Matrix<float> > (*)(const Matrix<float>&), const Matrix<float>&);
Matrix<float> bfgs(std::pair<float, Matrix<float> > (*)(const Matrix<float>&), const Matrix<float>&);
Matrix<float> lbfgs(std::pair<float, Matrix<float> > (*)(const Matrix<float>&), const Matrix<float>&);

Matrix<float> linprog_can(Matrix<float> f, Matrix<float> A, Matrix<float> b);
Matrix<float> linprog_std(Matrix<float> f, Matrix<float> A, Matrix<float> b);
Matrix<float> linprog(Matrix<float> f, Matrix<float> A, Matrix<float> b, Matrix<float> Aeq, Matrix<float> beq);
Matrix<float> quadprog(Matrix<float> H, Matrix<float> f, Matrix<float> A, Matrix<float> b, Matrix<float> Aeq, Matrix<float> beq, Matrix<float> X0);

// predefined function
std::pair<float, Matrix<float> > rosenbrock(const Matrix<float>&);
std::pair<float, Matrix<float> > himmelblau(const Matrix<float>&);
std::pair<float, Matrix<float> > rastrigin(const Matrix<float>&);

float ArmijoLineSearch(std::pair<float, Matrix<float> > (*ptr)(const Matrix<float>&), const Matrix<float>& X0, float err, const Matrix<float>& g0, const Matrix<float>& p);
float WolfeLineSearch(std::pair<float, Matrix<float> > (*ptr)(const Matrix<float>&), const Matrix<float>& X0, float err, const Matrix<float>& g0, const Matrix<float>& p);
Matrix<float> conjugate_gradient(const Matrix<float>& A, const Matrix<float>& b, const Matrix<float>& X0);

float zoom(std::pair<float, Matrix<float> > (*ptr)(const Matrix<float>&), const Matrix<float>& X0, const Matrix<float>& p, float alphal, float alphah);
float StrongWolfeLineSearch(std::pair<float, Matrix<float> > (*ptr)(const Matrix<float>&), const Matrix<float>& X0, float err, const Matrix<float>& g0, const Matrix<float>& p);


#endif // OPTIMIZATION_HPP
