#include <cmath>
#include <tuple>
#include <utility>
#include "mmath.hpp"
#include "logical_operators.hpp"
#include "statistics.hpp"


template <class T> Matrix<T> dot(const Matrix<T>& A, const Matrix<T>& B)
{
    if(A.colNb()==B.rowNb())
    {
        Matrix<T> my_matrix = zeros<T>(A.rowNb(), B.colNb());
        Matrix<T> Bt = transpose(B);

        auto A_begin = A.cbegin();
        auto B_begin = Bt.cbegin();
        std::size_t n = A.colNb();

        for(std::size_t i=0;i<my_matrix.rowNb();++i)
        {
            for(std::size_t j=0;j<my_matrix.colNb();++j)
                my_matrix(i, j) = std::inner_product(A_begin+i*n, A_begin+(i+1)*n, B_begin+j*n, T(0));
        }

        return my_matrix;
    }
    std::cout<<"dimension mismatch"<<std::endl;
    return A;
}

template <class T> Matrix<T> inv(const Matrix<T>& M)
{
    if(M.rowNb()!=M.colNb())
    {
        std::cout<<"the matrix is not a square matrix"<<std::endl;
        return M;
    }
    T Mdet = det(M);
    if(Mdet==0)
    {
        std::cout<<"the matrix is not invertible"<<std::endl;
        return M;
    }
    else
    {
        Matrix<T> new_mat(M.rowNb());
        if(M.rowNb()==1)
            new_mat(0, 0) = T(1);
        else if(M.rowNb()==2)
        {
            new_mat(0, 0) = M(1, 1);
            new_mat(0, 1) = -M(0, 1);
            new_mat(1, 0) = -M(1, 0);
            new_mat(1, 1) = M(0, 0);
        }
        else
        {
            for(std::size_t i=0;i<M.rowNb();++i)
            {
                for(std::size_t j=0;j<M.colNb();++j)
                {
                    Matrix<T> temp = M;
                    temp.delRow(i);
                    temp.delCol(j);
                    new_mat(j, i) = ((i+j)%2==0)?det(temp):-det(temp);
                }
            }
        }
        return new_mat/Mdet;
    }
}

template <class T> T det(const Matrix<T>& M)
{
    if(M.rowNb()!=M.colNb())
    {
        std::cout<<"the matrix is not a square matrix"<<std::endl;
        return T(0);
    }
    else if(M.rowNb()==1)
        return M(0, 0);
    else if(M.rowNb()==2)
        return M(0, 0)*M(1, 1)-M(0, 1)*M(1, 0);
    else
    {
        T determinant = T(1);
        Matrix<T> Mcopy = M;
        for(std::size_t i=0;i<M.rowNb()-2;++i)
        {
            Matrix<T> cur_col = Mcopy.getSubmat(i, Mcopy.rowNb(), i, Mcopy.rowNb());
            cur_col = abs(cur_col);
            Matrix<std::size_t> maxidx = argmax(cur_col);
            if(cur_col(maxidx(0, 0), maxidx(0, 1)) == T(0))
                return T(0);
            else
            {
                if(maxidx(0, 0)!=0)
                {
                    Mcopy.swaprow(i, maxidx(0, 0)+i);
                    determinant = -determinant;
                }
                if(maxidx(0, 1)!=0)
                {
                    Mcopy.swapcol(i, maxidx(0, 1)+i);
                    determinant = -determinant;
                }
                determinant*=Mcopy(i, i);
                for(std::size_t j=i+1;j<Mcopy.rowNb();++j)
                    Mcopy.setRow(j, (Mcopy.getRow(j)*Mcopy(i, i)-Mcopy.getRow(i)*Mcopy(j, i))/Mcopy(i, i));
            }
        }
        return determinant*det(Mcopy.getSubmat(Mcopy.rowNb()-2, Mcopy.rowNb(), Mcopy.rowNb()-2, Mcopy.rowNb()));
    }
}

template <class T> T norm(const Matrix<T>& M)
{
    return std::sqrt(sum(M*M));
}

template <class T> Matrix<T> bwdsub(const Matrix<T>& U, const Matrix<T>& B)
{
    if(U.rowNb()==U.colNb() && B.rowNb()==U.rowNb())
    {
        Matrix<T> X(B.rowNb(), B.colNb());
        std::size_t n = B.rowNb()-1;
        X.setRow(n, B.getRow(n)/U(n, n));
        for(std::size_t i=n;i>0;--i)
        {
            Matrix<T> new_row = B.getRow(i-1);
            for(std::size_t j=n;j>i-1;--j)
                new_row -= X.getRow(j)*U(i-1, j);
            X.setRow(i-1, new_row/U(i-1, i-1));
        }
        return X;
    }
    std::cout<<"dimension mismatch"<<std::endl;
    return zeros<T>(1);
}

template <class T> Matrix<T> fwdsub(const Matrix<T>& L, const Matrix<T>& B)
{
    if(L.rowNb()==L.colNb() && B.rowNb()==L.rowNb())
    {
        Matrix<T> X(B.rowNb(), B.colNb());
        X.setRow(0, B.getRow(0)/L(0,0));
        for(std::size_t i=1;i<X.rowNb();++i)
        {
            X.setRow(i, B.getRow(i));
            for(std::size_t j=0;j<i;++j)
                X.setRow(i, X.getRow(i)-X.getRow(j)*L(i, j));
            X.setRow(i, X.getRow(i)/L(i, i));
        }
        return X;
    }
    std::cout<<"dimension mismatch"<<std::endl;
    return zeros<T>(1);
}

template <class T> Matrix<T> cholesky(const Matrix<T>& M)
{
    if(M.rowNb()==M.colNb())
    {
        Matrix<T> L(M.rowNb(), M.colNb());
        for(std::size_t j=0;j<M.colNb();++j)
        {
            Matrix<T> V = M.getSubmat(j, M.rowNb(), j, j+1);
            for(std::size_t i=0;i<j;++i)
                V-=L(j, i)*L.getSubmat(j, M.rowNb(), i, i+1);
            if(V(0, 0)>0)
                L.setSubmat(j, j, V/sqrt(V(0, 0)));
            else
            {
                std::cout<<"the matrix is not positive definite"<<std::endl;
                return zeros<T>(1);
            }
        }
        return L;
    }
    std::cout<<"the matrix is not a square matrix"<<std::endl;
    return zeros<T>(1);
}

template <class T> std::pair<Matrix<T>, Matrix<T> > crout(const Matrix<T>& M)
{
    if(M.rowNb()==M.colNb())
    {
        Matrix<T> D(M.rowNb(), M.colNb());
        Matrix<T> L = id<T>(M.rowNb(), M.colNb());
        D(0, 0) = M(0, 0);
        L.setSubmat(1, 0, (M.getSubmat(1, M.rowNb(), 0, 1))/M(0, 0));
        for(std::size_t j=1;j<M.colNb();++j)
        {
            Matrix<T> V(j+1, 1);
            for(std::size_t i=0;i<j;++i)
                V(i, 0) = L(j, i)*D(i, i);
            V(j, 0) = M(j, j)-dot(L.getSubmat(j, j+1, 0, j+1), V.getRows(0, j+1))(0, 0);
            D(j, j) = V(j, 0);
            if(j!=M.colNb()-1)
                L.setSubmat(j+1, j, (M.getSubmat(j+1, M.rowNb(), j, j+1)-dot(L.getSubmat(j+1, M.rowNb(), 0, j), V.getSubmat(0, j, 0, 1)))/V(j, 0));
        }
        return std::make_pair(L, D);
    }
    std::cout<<"the matrix is not square"<<std::endl;
    return std::make_pair(zeros<T>(1), zeros<T>(1));
}

template <class T> std::tuple<Matrix<T>, Matrix<T>, Matrix<T> > lu(const Matrix<T>& M)
{
    Matrix<T> L(M.rowNb(), M.colNb());
    Matrix<T> U = M;
    Matrix<T> P = id<T>(M.rowNb(), M.colNb());
    for(std::size_t i=0;i<M.rowNb()-1;++i)
    {
        Matrix<T> cur_col = abs(U.getSubmat(i, M.rowNb(), i, i+1));
        Matrix<std::size_t> maxidx = argmax(cur_col);
        P.swaprow(i, maxidx(0, 0)+i);
        U.swaprow(i, maxidx(0, 0)+i);
        L.swaprow(i, maxidx(0, 0)+i);
        L.setSubmat(i+1, i, U.getSubmat(i+1, M.rowNb(), i, i+1)/U(i, i));
        L(i, i) = T(1);
        for(std::size_t j=i+1;j<M.rowNb();++j)
            U.setRow(j, (U.getRow(j)*U(i, i)-U.getRow(i)*U(j, i))/U(i, i));
    }
    L(M.rowNb()-1, M.colNb()-1) = T(1);
    return std::make_tuple(L, U, P);
}

template <class T> std::pair<Matrix<T>, Matrix<T> > qr(const Matrix<T>& M)
{
    Matrix<T> Q(M.rowNb(), M.colNb());
    Matrix<T> R(M.colNb(), M.colNb());
    for(std::size_t i=0;i<M.colNb();++i)
    {
        Matrix<T> U = M.getCol(i);
        for(std::size_t j=0;j<i;++j)
        {
            Matrix<T> cur_col = Q.getCol(j);
            R(j, i) = sum(U*cur_col);
            U-=R(j, i)*Q.getCol(j);
        }
        U/=norm(U);
        R(i, i) = sum(M.getCol(i)*U);
        Q.setCol(i, U);
    }
    return std::make_pair(Q, R);
}

template <class T> std::pair<Matrix<T>, Matrix<T> > svd(const Matrix<T>& M)
{
    Matrix<float> A = M;
    return std::make_pair(zeros<T>(1), zeros<T>(1));
}

template <class T> T trace(const Matrix<T>& M)
{
    std::size_t minimum = M.rowNb()>M.colNb()?M.colNb():M.rowNb();
    T tr = T(0);
    for(std::size_t i=0;i<minimum;++i)
        tr+=M(i, i);
    return tr;
}

template <class T> std::pair<Matrix<T>, Matrix<T> > Jacobi(const Matrix<T>& M)
{
    Matrix<T> A = M;
    Matrix<T> P = id<T>(M.rowNb(), M.colNb());
    while(true)
    {
        Matrix<std::size_t> Sortedvalues = argsort(abs(A));
        std::size_t i = Sortedvalues.rowNb()-1;
        Matrix<std::size_t> cur_value = Sortedvalues.getRow(i);
        while(cur_value(0, 0) >= cur_value(0,1))
        {
            --i;
            cur_value = Sortedvalues.getRow(i);
        }
        std::size_t p = cur_value(0, 0);
        std::size_t q = cur_value(0, 1);
        if(std::abs(A(p,q))<10e-10)
        {
            A(p, q) = 0;
            A(q, p) = 0;
            return std::make_pair(A, P);
        }

        T e = (A(q, q)-A(p, p))/(2*A(p, q));
        T t;
        if(e>0)
            t = -e+std::sqrt(1+e*e);
        else if(e<0)
            t = -e-std::sqrt(1+e*e);
        else
            t = 1;
        T cos = 1.0/std::sqrt(1+t*t);
        T sin = t/std::sqrt(1+t*t);

        Matrix<T> rowp = A.getRow(p);
        Matrix<T> rowq = A.getRow(q);
        Matrix<T> colp = A.getCol(p);
        Matrix<T> colq = A.getCol(q);
        T aqq = A(q, q);
        T app = A(p, p);
        T apq = A(p, q);
        A.setRow(p, cos*rowp-sin*rowq);
        A.setCol(q, sin*colp+cos*colq);
        rowp.transpose();
        rowq.transpose();
        colp.transpose();
        colq.transpose();
        A.setCol(p, cos*rowp-sin*rowq);
        A.setRow(q, sin*colp+cos*colq);
        A(p, q) = 0;
        A(q, p) = 0;
        A(q, q) = aqq+t*apq;
        A(p, p) = app-t*apq;

        Matrix<T> Pp = P.getCol(p);
        Matrix<T> Pq = P.getCol(q);
        P.setCol(p, cos*Pp-sin*Pq);
        P.setCol(q, cos*Pq+sin*Pp);
    }
}
