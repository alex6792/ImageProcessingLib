#include "matrix.hpp"
#include "mmath.hpp"


template <class T> Matrix<T> abs(const Matrix<T>& M)
{
    return apply<T,T>(M, std::abs);
}

template <class T> Matrix<T> max(const Matrix<T>& A, const Matrix<T>& B)
{
    return where(A>B, A, B);
}

template <class T> Matrix<T> min(const Matrix<T>& A, const Matrix<T>& B)
{
    return where(A<B, A, B);
}


// trigonometric functions
template <class T> Matrix<T> cos(const Matrix<T>& M)
{
    return apply<T,T>(M, std::cos);
}

template <class T> Matrix<T> sin(const Matrix<T>& M)
{
    return apply<T,T>(M, std::sin);
}

template <class T> Matrix<T> tan(const Matrix<T>& M)
{
    return apply<T,T>(M, std::tan);
}

template <class T> Matrix<T> acos(const Matrix<T>& M)
{
    return apply<T,T>(M, std::acos);
}

template <class T> Matrix<T> asin(const Matrix<T>& M)
{
    return apply<T,T>(M, std::asin);
}

template <class T> Matrix<T> atan(const Matrix<T>& M)
{
    return apply<T,T>(M, std::atan);
}

template <class T> Matrix<T> atan2(const Matrix<T>& My, const Matrix<T>& Mx)
{
    Matrix<T> new_mat(My.rowNb(), My.colNb());
    if(My.rowNb()==Mx.rowNb() && My.rowNb()==Mx.colNb())
        std::transform(My.begin(), My.end(), Mx.begin(), new_mat.begin(), std::atan2);
    else
        std::cout<<"dimension mismatch"<<std::endl;
    return new_mat;
}


// hyperbolic functions
template <class T> Matrix<T> cosh(const Matrix<T>& M)
{
    return apply<T,T>(M, std::cosh);
}

template <class T> Matrix<T> sinh(const Matrix<T>& M)
{
    return apply<T,T>(M, std::sinh);
}

template <class T> Matrix<T> tanh(const Matrix<T>& M)
{
    return apply<T,T>(M, std::tanh);
}

template <class T> Matrix<T> acosh(const Matrix<T>& M)
{
    return apply<T,T>(M, std::acosh);
}

template <class T> Matrix<T> asinh(const Matrix<T>& M)
{
    return apply<T,T>(M, std::sinh);
}

template <class T> Matrix<T> atanh(const Matrix<T>& M)
{
    return apply<T,T>(M, std::atanh);
}


// exponential and logarithmic functions
template <class T> Matrix<T> exp(const Matrix<T>& M)
{
    return apply<T,T>(M, std::exp);
}

template <class T> Matrix<T> exp2(const Matrix<T>& M)
{
    return apply<T,T>(M, std::exp2);
}

template <class T> Matrix<T> log(const Matrix<T>& M)
{
    return apply<T,T>(M, std::log);
}

template <class T> Matrix<T> log2(const Matrix<T>& M)
{
    return apply<T,T>(M, std::log2);
}

template <class T> Matrix<T> log10(const Matrix<T>& M)
{
    return apply<T,T>(M, std::log10);
}


// power functions
template <class T> Matrix<T> pow(const Matrix<T>& M, const T& exponent)
{
    Matrix<T> new_mat(M.rowNb(), M.colNb());
    std::transform(M.cbegin(), M.cend(), new_mat.begin(), [exponent](T value){return std::pow(value, exponent);});
    return new_mat;
}

template <class T> Matrix<T> sqrt(const Matrix<T>& M)
{
    return apply<T,T>(M, std::sqrt);
}
