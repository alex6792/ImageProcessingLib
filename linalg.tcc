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
        std::size_t H = A.rowNb(), W = B.colNb();
        Matrix<T> my_matrix = zeros<T>(H, W);
        Matrix<T> Bt = transpose(B);

        auto A_begin = A.cbegin();
        auto B_begin = Bt.cbegin();
        std::size_t n = A.colNb();

        for(std::size_t i=0;i<H;++i)
        {
            for(std::size_t j=0;j<W;++j)
                my_matrix(i, j) = std::inner_product(A_begin+i*n, A_begin+(i+1)*n, B_begin+j*n, T(0));
        }

        return my_matrix;
    }
    std::cout<<"dimension mismatch"<<std::endl;
    return A;
}

template <class T> Matrix<T> inv(const Matrix<T>& M)
{
    std::size_t H = M.rowNb(), W = M.colNb();
    if(H!=W)
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
        Matrix<T> new_mat(H);
        if(H==1)
            new_mat(0, 0) = T(1);
        else if(H==2)
        {
            new_mat(0, 0) = M(1, 1);
            new_mat(0, 1) = -M(0, 1);
            new_mat(1, 0) = -M(1, 0);
            new_mat(1, 1) = M(0, 0);
        }
        else
        {
            for(std::size_t i=0;i<H;++i)
            {
                for(std::size_t j=0;j<W;++j)
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
    std::size_t H = M.rowNb(), W = M.colNb();
    if(H!=W)
    {
        std::cout<<"the matrix is not a square matrix"<<std::endl;
        return T(0);
    }
    else if(H==1)
        return M(0, 0);
    else if(H==2)
        return M(0, 0)*M(1, 1)-M(0, 1)*M(1, 0);
    else
    {
        T determinant = T(1);
        Matrix<T> Mcopy = M;
        for(std::size_t i=0;i<H-2;++i)
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
        for(std::size_t i=1, I=X.rowNb();i<I;++i)
        {
            Matrix<T> new_row = B.getRow(i);
            for(std::size_t j=0;j<i;++j)
                new_row-=X.getRow(j)*L(i, j);
            X.setRow(i, new_row/L(i, i));
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
        for(std::size_t j=0, J=M.rowNb();j<J;++j)
        {
            Matrix<T> V = M.getSubmat(j, J, j, j+1);
            for(std::size_t i=0;i<j;++i)
                V-=L(j, i)*L.getSubmat(j, J, i, i+1);
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

template <class T> std::pair<Matrix<T>, Matrix<T> > crout(const Matrix<T>& M)//decomposition LDU -> see doolittle
{
    std::size_t H = M.rowNb(), W = M.colNb();
    if(H==W)
    {
        Matrix<T> D(H, W);
        Matrix<T> L = id<T>(H, W);
        D(0, 0) = M(0, 0);
        L.setSubmat(1, 0, (M.getSubmat(1, H, 0, 1))/M(0, 0));
        for(std::size_t j=1;j<W;++j)
        {
            Matrix<T> V(j+1, 1);
            for(std::size_t i=0;i<j;++i)
                V(i, 0) = L(j, i)*D(i, i);
            V(j, 0) = M(j, j)-dot(L.getSubmat(j, j+1, 0, j+1), V.getRows(0, j+1))(0, 0);
            D(j, j) = V(j, 0);
            if(j!=W-1)
                L.setSubmat(j+1, j, (M.getSubmat(j+1, H, j, j+1)-dot(L.getSubmat(j+1, H, 0, j), V.getSubmat(0, j, 0, 1)))/V(j, 0));
        }
        return std::make_pair(L, D);
    }
    std::cout<<"the matrix is not square"<<std::endl;
    return std::make_pair(zeros<T>(1), zeros<T>(1));
}

template <class T> std::tuple<Matrix<T>, Matrix<T>, Matrix<T> > lu(const Matrix<T>& M)
{
    std::size_t H = M.rowNb(), W = M.colNb();
    Matrix<T> L(H, W);
    Matrix<T> U = M;
    Matrix<T> P = id<T>(H, W);
    for(std::size_t i=0;i<H-1;++i)
    {
        Matrix<T> cur_col = abs(U.getSubmat(i, H, i, i+1));
        Matrix<std::size_t> maxidx = argmax(cur_col);
        P.swaprow(i, maxidx(0, 0)+i);
        U.swaprow(i, maxidx(0, 0)+i);
        L.swaprow(i, maxidx(0, 0)+i);
        L.setSubmat(i+1, i, U.getSubmat(i+1, H, i, i+1)/U(i, i));
        L(i, i) = T(1);
        Matrix<float> rowi = U.getRow(i);
        for(std::size_t j=i+1;j<H;++j)
            U.setRow(j, (U.getRow(j)*U(i, i)-rowi*U(j, i))/U(i, i));
    }
    L(H-1, W-1) = T(1);
    return std::make_tuple(L, U, P);
}


template <class T> Matrix<T> pinv(const Matrix<T>& M)
{
    std::size_t H = M.rowNb(), W = M.colNb();
    bool tr = false;
    Matrix<T> A;
    if(H<W)
    {
        tr = true;
        A = dot(M, transpose(M));
        W = H;
    }
    else
        A = dot(transpose(M), M);

    Matrix<T> da = A.getDiag();
    T tol = 10e-9;
    Matrix<T> L = zeros<T>(W, W);
    std::size_t r = 0;
    for(std::size_t k=0;k<W;++k)
    {
        Matrix<T> temp = A.getSubmat(k, W, k, k+1);
        if(r>0)
            temp-=dot(L.getSubmat(k, W, 0, r), transpose(L.getSubmat(k, k+1, 0, r)));
        L.setSubmat(k, r, temp);
        if(L(k,r)>tol)
        {
            L(k,r) = std::sqrt(L(k,r));
            if(k+1<W)
                L.setSubmat(k+1, r, L.getSubmat(k+1, W, r, r+1)/L(k, r));
            ++r;
        }
    }

    L = L.getCols(0,r);
    Matrix<T> L_t = transpose(L);
    Matrix<T> N = inv(dot(L_t, L));
    Matrix<T> X = dot(L, N);
    Matrix<T> Y = dot(N, L_t);
    Matrix<T> Z = dot(X, Y);
    if(tr)
        return dot(transpose(M), Z);
    else
        return dot(Z, transpose(M));
}

template <class T> std::pair<Matrix<T>, Matrix<T> > qr(const Matrix<T>& M)
{
    std::size_t H = M.rowNb(), W = M.colNb();
    Matrix<T> Q, R;
    Q = zeros<T>(H, W);
    R = zeros<T>(W, W);
    std::vector<Matrix<float> > V(W);
    for(std::size_t i=0;i<W;++i)
    {
        V[i] = M.getCol(i);
    }
    for(std::size_t i=0;i<W;++i)
    {
        R(i, i) = norm(V[i]);
        Q.setCol(i, V[i]/R(i, i));
        Matrix<T> cur_col = Q.getCol(i);
        for(std::size_t j=i+1;j<W;++j)
        {
            R(i, j) = sum(V[j]*cur_col);
            V[j]-=R(i, j)*cur_col;
        }
    }
    return std::make_pair(Q, R);
}

template <class T> std::pair<Matrix<T>, Matrix<T> > rq(const Matrix<T>& M)
{
    std::size_t H = M.rowNb(), W = M.colNb();
    Matrix<T> Q, R;
    Q = zeros<T>(H, W);
    R = zeros<T>(H, H);
    std::vector<Matrix<float> > V(H);
    for(std::size_t i=0;i<H;++i)
    {
        V[i] = M.getRow(i);
    }
    for(std::size_t i=H;i>0;--i)
    {
        R(i-1, i-1) = norm(V[i-1]);
        Q.setRow(i-1, V[i-1]/R(i-1, i-1));
        Matrix<T> cur_row = Q.getRow(i-1);
        for(std::size_t j=H;j>i-1;--j)
        {
            R(i-1, j-1) = sum(V[j-1]*cur_row);
            V[j-1]-=R(i-1, j-1)*cur_row;
        }
    }
    return std::make_pair(R, Q);
}

template <class T> std::pair<Matrix<T>, Matrix<T> > ql(const Matrix<T>& M)
{
    std::size_t H = M.rowNb(), W = M.colNb();
    Matrix<T> Q, L;
    Q = zeros<T>(H, W);
    L = zeros<T>(W, W);
    std::vector<Matrix<float> > V(W);
    for(std::size_t i=0;i<W;++i)
    {
        V[i] = M.getCol(i);
    }
    for(std::size_t i=W;i>0;--i)
    {
        L(i-1, i-1) = norm(V[i-1]);
        Q.setCol(i-1, V[i-1]/L(i-1, i-1));
        Matrix<T> cur_col = Q.getCol(i-1);
        for(std::size_t j=W;j>i-1;--j)
        {
            L(j-1, i-1) = sum(V[j-1]*cur_col);
            V[j-1]-=L(j-1, i-1)*cur_col;
        }
    }
    return std::make_pair(Q, L);
}

template <class T> std::pair<Matrix<T>, Matrix<T> > lq(const Matrix<T>& M)
{
    std::size_t H = M.rowNb(), W = M.colNb();
    Matrix<T> Q, L;
    Q = zeros<T>(H, W);
    L = zeros<T>(H, H);
    std::vector<Matrix<float> > V(H);
    for(std::size_t i=0;i<H;++i)
    {
        V[i] = M.getRow(i);
    }
    for(std::size_t i=0;i<H;++i)
    {
        L(i, i) = norm(V[i]);
        Q.setRow(i, V[i]/L(i, i));
        Matrix<T> cur_row = Q.getRow(i);
        for(std::size_t j=i+1;j<H;++j)
        {
            L(j, i) = sum(V[j]*cur_row);
            V[j]-=L(j, i)*cur_row;
        }
    }
    return std::make_pair(L, Q);
}

template <class T> std::tuple<Matrix<T>, Matrix<T>, Matrix<T> > svd(const Matrix<T>& M)
{
    std::size_t H = M.rowNb();
    std::size_t W = M.colNb();//H>=W

    Matrix<float> Mt = transpose(M);
    Matrix<float> temp1 = dot(Mt, M);
    Matrix<float> temp2 = dot(M, Mt);

    Matrix<float> U(H, H);
    Matrix<float> V(W, W);

    auto J1 = jacobi(temp1);
    auto J2 = jacobi(temp2);

    if(H>W)
    {
        Matrix<float> S(W, W);
        std::tie(S, V) = J1;
        U = J2.second;
        for(std::size_t i=0;i<H-W;++i)
            U.delCol(W);
        return std::make_tuple(U,sqrt(S),V);
    }
    else
    {
        Matrix<float> S(H, H);
        std::tie(S, U) = J2;
        V = J1.second;
        for(std::size_t i=0;i<W-H;++i)
            S.newCol();
        return std::make_tuple(U,sqrt(S),V);
    }
}

template <class T> T trace(const Matrix<T>& M)
{
    return sum(M.getDiag());
}

template <class T> std::pair<Matrix<T>, Matrix<T> > jacobi(const Matrix<T>& M)
{
    std::size_t H = M.rowNb();
    std::size_t W = M.colNb();
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
        if(std::abs(A(p,q))<10e-9)
        {
            A(p, q) = 0;
            A(q, p) = 0;
            Matrix<T> dA = A.getDiag();
            std::size_t n = dA.size();
            A*=id<T>(n);
            Matrix<std::size_t> sorted_idx = argsort(dA);
            sorted_idx.flipud();
            for(std::size_t i=0;i<n;++i)
            {
                std::size_t cur_idx = sorted_idx(i, 0);
                std::swap(A(i,i), A(cur_idx, cur_idx));
                P.swapcol(i, cur_idx);

                Matrix<T> dA = A.getDiag();
                sorted_idx = argsort(dA);
                sorted_idx.flipud();
            }

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


template <class T> Matrix<T> vander(const Matrix<T>& M)
{
    std::size_t H = M.rowNb(), W = M.colNb();
    if(H==1)
    {
        Matrix<T> V(W, W);
        V.setRow(0, ones<T>(1, W));
        for(std::size_t i=1;i<W;++i)
            V.setRow(i, pow(M, T(i)));
        return V;
    }
    else if(W==1)
    {
        Matrix<T> V(H, H);
        V.setCol(0, ones<T>(H, 1));
        for(std::size_t i=1;i<H;++i)
            V.setCol(i, pow(M, T(i)));
        return V;
    }
    else
    {
        std::cout<<"The argument must be a vector"<<std::endl;
        return Matrix<T>();
    }
}


template <class T> std::tuple<Matrix<T>, Matrix<T>, Matrix<T> > bidiag_reduction(const Matrix<T>& M)
{
    std::size_t m = M.rowNb(), n = M.colNb();

    Matrix<T> B = M;
    Matrix<T> U = id<T>(m);
    Matrix<T> V = id<T>(n);
    Matrix<T> H;

    for(std::size_t i=0;i<n;++i)
    {
        H = householder(B.getCol(i),i);
        B = dot(H, B);
        U = dot(U, H);
        if(i<n-2)
        {
            H = householder(transpose(B.getRow(i)),i+1);
            B = dot(B, H);
            V = dot(H, V);
        }
    }

    return std::make_tuple(U, B, V);
}


template <class T> Matrix<T> householder(const Matrix<T>& M, std::size_t idx)
{
    std::size_t m = M.rowNb(), n = M.colNb();
    Matrix<T> d;
    if(m>1)
        d = M.getRows(idx, m);
    else
        d = M.getCols(idx, n);

    T alpha;
    if(d(0,0)>=0)
        alpha = -norm(d);
    else
        alpha = norm(d);

    if(std::abs(alpha)<10e-9)
        return id<T>(m*n);

    std::size_t len_d = d.size();
    Matrix<T> v = zeros<T>(len_d, 1);
    v(0, 0) = std::sqrt(0.5*(1-d(0,0)/alpha));
    T p = -alpha*v(0,0);
    v.setRows(1, d.getRows(1,len_d)/(2*p));//warning
    Matrix<T> w = zeros<T>(m*n, 1);
    w.setRows(idx, v);
    return id<T>(m*n)-T(2)*dot(w, transpose(w));
}
