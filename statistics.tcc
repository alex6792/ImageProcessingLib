#include <algorithm>
#include <cmath>
#include "mmath.hpp"


template <class T> Matrix<std::size_t> argmax(const Matrix<T>& M)
{
    Matrix<std::size_t> argmaxM = zeros<std::size_t>(1, 2);
    auto it = std::max_element(M.cbegin(), M.cend());
    std::size_t argmax_idx = it-M.cbegin();
    std::size_t w = M.colNb();
    argmaxM(0, 0) = argmax_idx/w;
    argmaxM(0, 1) = argmax_idx%w;
    return argmaxM;
}

template <class T> Matrix<std::size_t> argmin(const Matrix<T>& M)
{
    Matrix<std::size_t> argminM = zeros<std::size_t>(1, 2);
    auto it = std::min_element(M.cbegin(), M.cend());
    std::size_t argmin_idx = it-M.cbegin();
    std::size_t w = M.colNb();
    argminM(0, 0) = argmin_idx/w;
    argminM(0, 1) = argmin_idx%w;
    return argminM;
}

template <class T> Matrix<std::size_t> argsort(const Matrix<T>& M)
{
    Matrix<std::size_t> temp = arange<std::size_t>(0, M.size());
    std::size_t w = M.colNb();
    auto it = M.cbegin();
    std::sort(temp.begin(), temp.end(), [it](int i1, int i2) {return *(it+i1) < *(it+i2);});
    Matrix<std::size_t> args = Matrix<std::size_t>(M.size(), 2);
    args.setCol(0, temp/w);
    args.setCol(1, temp%w);
    return args;
}

template <class T> Matrix<T> axismax(const Matrix<T>& M, int axis)
{
    if(axis==0)
        return max(M)*ones<T>(1);
    else if(axis==1)
    {
        Matrix<T> maximum = M.getCol(0);
        for(std::size_t j=1, J=M.colNb();j<J;++j)
            maximum = max(maximum, M.getCol(j));
        return maximum;
    }
    else if(axis==2)
    {
        Matrix<T> maximum = M.getRow(0);
        for(std::size_t i=1, I=M.rowNb();i<I;++i)
            maximum = max(maximum, M.getRow(i));
        return maximum;
    }
    else
    {
        std::cout<<"invalid axis"<<std::endl;
        return zeros<T>(1);
    }
}

template <class T> Matrix<T> axismean(const Matrix<T>& M, int axis)
{
    if(axis==0)
        return mean(M)*ones<T>(1);
    else if(axis==1)
        return axissum(M, 1)/(T)M.colNb();
    else if(axis==2)
        return axissum(M, 2)/(T)M.rowNb();
    else
    {
        std::cout<<"invalid axis"<<std::endl;
        return zeros<T>(1);
    }
}

template <class T> Matrix<T> axismin(const Matrix<T>& M, int axis)
{
    if(axis==0)
        return min(M)*ones<T>(1);
    else if(axis==1)
    {
        Matrix<T> minimum = M.getCol(0);
        for(std::size_t j=1, J=M.colNb();j<J;++j)
            minimum = min(minimum, M.getCol(j));
        return minimum;
    }
    else if(axis==2)
    {
        Matrix<T> minimum = M.getRow(0);
        for(std::size_t i=1, I=M.rowNb();i<I;++i)
            minimum = min(minimum, M.getRow(i));
        return minimum;
    }
    else
    {
        std::cout<<"invalid axis"<<std::endl;
        return zeros<T>(1);
    }
}

template <class T> Matrix<T> axisprod(const Matrix<T>& M, int axis)
{
    if(axis==0)
        return prod(M)*ones<T>(1);
    else if(axis==1)
    {
        Matrix<T> prod = M.getCol(0);
        for(std::size_t j=1;j<M.colNb();++j)
            prod*=M.getCol(j);
        return prod;
    }
    else if(axis==2)
    {
        Matrix<T> prod = M.getRow(0);
        for(std::size_t i=1;i<M.rowNb();++i)
            prod*=M.getRow(i);
        return prod;
    }
    else
    {
        std::cout<<"invalid axis"<<std::endl;
        return zeros<T>(1);
    }
}

template <class T> Matrix<T> axisstdev(const Matrix<T>& M, int axis)
{
    return sqrt(axisvar(M, axis));
}

template <class T> Matrix<T> axissum(const Matrix<T>& M, int axis)
{
    if(axis==0)
        return sum(M)*ones<T>(1);
    else if(axis==1)
    {
        Matrix<T> sum = M.getCol(0);
        for(std::size_t j=1;j<M.colNb();++j)
            sum+=M.getCol(j);
        return sum;
    }
    else if(axis==2)
    {
        Matrix<T> sum = M.getRow(0);
        for(std::size_t i=1;i<M.rowNb();++i)
            sum+=M.getRow(i);
        return sum;
    }
    else
    {
        std::cout<<"invalid axis"<<std::endl;
        return zeros<T>(1);
    }
}

template <class T> Matrix<T> axisvar(const Matrix<T>& M, int axis)
{
    if(axis==0)
        return var(M)*ones<T>(1);
    else if(axis==1)
        return axismean(M*M, 1)-pow(axismean(M, 1), T(2));
    else if(axis==2)
        return axismean(M*M, 2)-pow(axismean(M, 2), T(2));
    else
    {
        std::cout<<"invalid axis"<<std::endl;
        return zeros<T>(1);
    }
}

template <class T> T central_moment(const Matrix<T>& M, int order)
{
    return moment(M-mean(M), order);
}

template <class T> T central_moment(const Matrix<T>& M, int orderx, int ordery)
{
    auto XY = meshgrid<T>(M.rowNb(), M.colNb());
    Matrix<T>& X = XY.first;
    Matrix<T>& Y = XY.second;
    T M00 = sum(M);
    if(std::abs(M00)!=T(0))
    {
        T M10 = sum(M*X);
        T M01 = sum(M*Y);
        X-=M10/M00;
        Y-=M01/M00;
        return sum(M*pow(X, T(orderx))*pow(Y, T(ordery)));
    }
    else
        return T(0);
}

template <class T> T corrected_var(const Matrix<T>& M)
{
    if(M.size()==1)
        return T(0);
    else
        return var(M)*M.size()/(T)(M.size()-1);
}

template <class T> Matrix<T> cumsum(const Matrix<T>& M)
{
    Matrix<T> new_mat(M.rowNb(), M.colNb());
    std::partial_sum(M.cbegin(), M.cend(), new_mat.begin());
    return new_mat;
}


template <class T> T geometric_mean(const Matrix<T>& M)
{
    return std::pow(prod(M),T(1)/T(M.size()));
}

template <class T> T max(const Matrix<T>& M)
{
    return *std::max_element(M.cbegin(), M.cend());
}

template <class T> T mean(const Matrix<T>& M)
{
    return sum(M)/(T)M.size();
}

template <class T> T median(const Matrix<T>& M)
{
    std::deque<T> M_copy(M.cbegin(), M.cend());
    if(M_copy.size()%2==1)
    {
        std::nth_element(M_copy.begin(), M_copy.begin()+(M_copy.size()-1)/2, M_copy.end());
        return M_copy[(M_copy.size()-1)/2];
    }
    else
    {
        std::nth_element(M_copy.begin(), M_copy.begin()+(M_copy.size())/2-1, M_copy.end());
        std::nth_element(M_copy.begin()+(M_copy.size())/2, M_copy.begin()+(M_copy.size())/2, M_copy.end());
        return (M_copy[M_copy.size()/2]+M_copy[M_copy.size()/2-1])/2.0;
    }
}

template <class T> T min(const Matrix<T>& M)
{
    return *std::min_element(M.cbegin(), M.cend());
}

template <class T> T moment(const Matrix<T>& M, int order)
{
    return mean(pow(M, order));
}

template <class T> T moment(const Matrix<T>& M, int orderx, int ordery)
{
    auto XY = meshgrid(M.rowNb(), M.colNb());
    Matrix<T>& X = XY.first;
    Matrix<T>& Y = XY.second;
    return sum(M*pow(X, orderx)*pow(Y, ordery));
}

template <class T> T prod(const Matrix<T>& M)
{
    return std::accumulate(M.cbegin(), M.cend(), T(1), std::multiplies<T>());
}

template <class T> T stdev(const Matrix<T>& M)
{
    return std::sqrt(var(M));
}

template <class T> Matrix<T> sort(const Matrix<T>& M)
{
    Matrix<T> Mcopy = M;
    std::sort(Mcopy.begin(), Mcopy.end());
    return Mcopy;
}

template <class T> T sum(const Matrix<T>& M)
{
    return std::accumulate(M.cbegin(), M.cend(), T(0));
}

template <class T> T var(const Matrix<T>& M)
{
    T m = mean(M);
    return mean(M*M)-m*m;
}
