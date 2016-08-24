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
        Matrix<T> maximum = M.getCol(0);
        for(std::size_t j=1, J=M.colNb();j<J;++j)
            maximum = min(maximum, M.getCol(j));
        return maximum;
    }
    else if(axis==2)
    {
        Matrix<T> maximum = M.getRow(0);
        for(std::size_t i=1, I=M.rowNb();i<I;++i)
            maximum = min(maximum, M.getRow(i));
        return maximum;
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
        for(std::size_t i=0, I=M.rowNb();i<I;++i)
        {
            for(std::size_t j=1;j<M.colNb();++j)
                prod(i, 0)*=M(i, j);
        }
        return prod;
    }
    else if(axis==2)
    {
        Matrix<T> prod = M.getRow(0);
        for(std::size_t i=1;i<M.rowNb();++i)
        {
            for(std::size_t j=0, J=M.colNb();j<J;++j)
                prod(0, j)*=M(i, j);
        }
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
        for(std::size_t i=0, I=M.rowNb();i<I;++i)
        {
            for(std::size_t j=1;j<M.colNb();++j)
                sum(i, 0)+=M(i, j);
        }
        return sum;
    }
    else if(axis==2)
    {
        Matrix<T> sum = M.getRow(0);
        for(std::size_t i=1;i<M.rowNb();++i)
        {
            for(std::size_t j=0, J=M.colNb();j<J;++j)
                sum(0, j)+=M(i, j);
        }
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
        return full<T>(1, var(M));
    else if(axis==1)
        return axissum(M*M, 1)/(T)M.colNb()-pow(axismean(M, 1), T(2));
    else if(axis==2)
        return axissum(M*M, 2)/(T)M.rowNb()-pow(axismean(M, 2), T(2));
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
    auto XY = meshgrid(M.rowNb(), M.colNb());
    Matrix<T>& X = XY.first;
    Matrix<T>& Y = XY.second;
    T M00 = sum(M);
    if(std::abs(M00)!=T(0))
    {
        T M10 = sum(M*X);
        T M01 = sum(M*Y);
        X-=M10/M00;
        Y-=M01/M00;
        return sum(M*pow(X, orderx)*pow(Y, ordery));
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
    std::deque<T> M_copy(M.begin(), M.end());
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

template <class T> int quickselect(const Matrix<T>& M, int order)
{
    Matrix<T> Mcopy = M;
    return quickselect(Mcopy, 0, M.size()-1, order);
}

template <class T> Matrix<T> quicksort(const Matrix<T>& M)
{
    Matrix<T> Mcopy = M;
    quicksort(Mcopy, 0, M.size()-1);
    return Mcopy;
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


template <class T> static int partition(Matrix<T>& M, int first, int last, int pivot)
{
    std::swap(M(pivot/M.colNb(), pivot%M.colNb()), M(last/M.colNb(), last%M.colNb()));
    int j = first;
    for(int i=first;i<last;++i)
    {
        if(M(i/M.colNb(), i%M.colNb())<M(last/M.colNb(), last%M.colNb()))
        {
            std::swap(M(i/M.colNb(), i%M.colNb()), M(j/M.colNb(), j%M.colNb()));
            ++j;
        }
    }
    std::swap(M(last/M.colNb(), last%M.colNb()), M(j/M.colNb(), j%M.colNb()));
    return j;
}

template <class T> static int quickselect(Matrix<T>& M, int first, int last, int order)
{
    if(first==last)
        return M(first/M.colNb(), first%M.colNb());

    int pivot = partition(M, first, last, first);
    if(order==pivot)
        return M(order/M.colNb(), order%M.colNb());
    else if(order<pivot)
        return quickselect(M, first, pivot-1, order);
    else
        return quickselect(M, pivot+1, last, order);
}

template <class T> static void quicksort(Matrix<T>& M, int first, int last)
{
    if(first<last)
    {
        int pivot = partition(M, first, last, first);
        quicksort(M, first, pivot-1);
        quicksort(M, pivot+1, last);
    }
}
