#pragma once

#ifndef POLYNOMIAL_EQUATION_SOLVER_HPP
#define POLYNOMIAL_EQUATION_SOLVER_HPP


#include <complex>
#include "matrix.hpp"


template <class T> Matrix<std::complex<T> > solve(const T& a, const T& b, const T& c)
{
    if(a!=0)
    {
        T delta = b*b-4*a*c;
        if(delta==0)
        {
            Matrix<std::complex<T> > M(1);
            M(0, 0) = std::complex<T>(-b/(2*a), T(0));
            return M;
        }
        else if(delta>0)
        {
            Matrix<std::complex<T> > M(2, 1);
            if(b<0)
            {
                M(0, 0) = std::complex<T>((-b+std::sqrt(delta))/(2*a), T(0));
                M(1, 0) = std::complex<T>(2*c/(-b+std::sqrt(delta)), T(0));
            }
            else
            {
                M(0, 0) = std::complex<T>((-b-std::sqrt(delta))/(2*a), T(0));
                M(1, 0) = std::complex<T>(2*c/(-b-std::sqrt(delta)), T(0));
            }
            return M;
        }
        else
        {
            Matrix<std::complex<T> > M(2, 1);
            M(0, 0) = std::complex<T>((-b)/(2*a), -std::sqrt(-delta)/(2*a));
            M(1, 0) = std::complex<T>((-b)/(2*a), std::sqrt(-delta)/(2*a));
            return M;
        }
    }
    else if(b!=0)
    {
        Matrix<std::complex<T> > M(1);
        M(0, 0) = std::complex<T>(-c/b, T(0));
        return M;
    }
    else
    {
        std::cout<<"No solution found"<<std::endl;
        return zeros<std::complex<T> >(1);
    }
}


#endif // POLYNOMIAL_EQUATION_SOLVER_HPP
