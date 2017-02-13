#include "matrix.hpp"

Matrix<float> linprog_can(Matrix<float> f, Matrix<float> A, Matrix<float> b);
Matrix<float> linprog_std(Matrix<float> f, Matrix<float> A, Matrix<float> b);
Matrix<float> linprog(Matrix<float> f, Matrix<float> A, Matrix<float> b, Matrix<float> Aeq, Matrix<float> beq);
Matrix<float> quadprog(Matrix<float> H, Matrix<float> f, Matrix<float> A, Matrix<float> b, Matrix<float> Aeq, Matrix<float> beq, Matrix<float> X0);
