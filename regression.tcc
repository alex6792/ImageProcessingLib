#include "linalg.hpp"
#include "optimization.hpp"
#include "statistics.hpp"


template <class T> Matrix<T> poly_regression(const Matrix<T>& Y, const Matrix<T>& X, std::size_t degree)
{
    if(X.colNb()==1 && X.rowNb()==Y.rowNb() && X.colNb()==Y.colNb())
    {
        if(degree>=0)
        {
            std::size_t n = X.rowNb();
            if(n>degree)
            {
                Matrix<T> A(n, degree+1);
                A.setCol(0, ones<T>(n, 1));
                for(std::size_t i=1;i<=degree;++i)
                    A.setCol(i, pow(X, T(i)));
                Matrix<T> X0 = zeros<T>(n,1);
                return solve(A, Y, X0);
            }
            else
                std::cout<<"not enough elements"<<std::endl;
        }
        else
            std::cout<<"invalid degree"<<std::endl;
    }
    else
        std::cout<<"invalid matrix size"<<std::endl;
    return Y;
}
