/*!
 * \author Alexandre Krebs
 * \file optimization.hpp
 * \brief optimization tools
 */


#pragma once

#ifndef OPTIMIZATION_HPP
#define OPTIMIZATION_HPP


#include <functional>
#include "matrix.hpp"

typedef std::function<std::pair<float, Matrix<float> >(const Matrix<float>&) > opti_func;

//standard algorithms
Matrix<float> steepest_descent(opti_func, const Matrix<float>&);
Matrix<float> bfgs(opti_func, const Matrix<float>&);
Matrix<float> lbfgs(opti_func, const Matrix<float>&);

//linear programming and quadratic programming
Matrix<float> linprog_can(Matrix<float> f, Matrix<float> A, Matrix<float> b);
Matrix<float> linprog_std(Matrix<float> f, Matrix<float> A, Matrix<float> b);
Matrix<float> linprog(Matrix<float> f, Matrix<float> A, Matrix<float> b, Matrix<float> Aeq, Matrix<float> beq);
Matrix<float> quadprog(Matrix<float> H, Matrix<float> f, Matrix<float> A, Matrix<float> b, Matrix<float> Aeq, Matrix<float> beq, Matrix<float> X0);
Matrix<float> conjugate_gradient(const Matrix<float>& A, const Matrix<float>& b, const Matrix<float>& X0);
Matrix<float> solve(const Matrix<float>& A, const Matrix<float>& b, const Matrix<float>& X0);
Matrix<float> solve(const Matrix<float>& A, const Matrix<float>& b, const Matrix<float>& X0, const Matrix<float>& lb, const Matrix<float>& ub);

//linear search
float ArmijoLineSearch(opti_func, const Matrix<float>& X0, float err, const Matrix<float>& g0, const Matrix<float>& p);
float WolfeLineSearch(opti_func, const Matrix<float>& X0, float err, const Matrix<float>& g0, const Matrix<float>& p);
float zoom(opti_func, const Matrix<float>& X0, const Matrix<float>& p, float alphal, float alphah);
float StrongWolfeLineSearch(opti_func, const Matrix<float>& X0, float err, const Matrix<float>& g0, const Matrix<float>& p);

// predefined function
std::pair<float, Matrix<float> > rosenbrock(const Matrix<float>&);
std::pair<float, Matrix<float> > himmelblau(const Matrix<float>&);
std::pair<float, Matrix<float> > rastrigin(const Matrix<float>&);


#endif // OPTIMIZATION_HPP
