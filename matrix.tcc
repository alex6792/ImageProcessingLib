#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <numeric>
#include <random>
#include <vector>



// constructors & accessors
template <class T> Matrix<T>::Matrix() : Matrix<T>::Matrix(1)
{
}

template <class T> Matrix<T>::Matrix(int a) : Matrix<T>::Matrix(a, a)
{
}

template <class T> Matrix<T>::Matrix(int a, int b)
{
    if(a>0 && b>0)
    {
        sizex = a;
        sizey = b;
        mat = std::deque<T>(a*b);
    }
    else
        std::cout<<"invalid matrix size {"<<a<<", "<<b<<"}"<<std::endl;
}

template <class T> template <typename Type> Matrix<T>::Matrix(const Matrix<Type>& M)
{
    sizex = M.rowNb();
    sizey = M.colNb();
    mat = std::deque<T>(M.cbegin(), M.cend());
}

template <class T> Matrix<T>::Matrix(std::initializer_list<std::initializer_list<T> > tab)
{
    auto it = tab.begin();
    mat = std::deque<T>(*it);
    ++it;
    sizex = 1;
    sizey = mat.size();
    while(it!=tab.end())
    {
        mat.insert(mat.end(), *it);
        ++it;
        ++sizex;
        if(mat.size()!= std::size_t(sizex*sizey) || mat.size()==0)
        {
            std::cout<<"invalid array"<<std::endl;
            sizex = 1;
            sizey = 1;
            mat = std::deque<T>(1);
            break;
        }
    }
}

template <class T> typename std::deque<T>::const_iterator Matrix<T>::cbegin() const
{
    return mat.cbegin();
}

template <class T> typename std::deque<T>::const_iterator Matrix<T>::cend() const
{
    return mat.cend();
}

template <class T> typename std::deque<T>::iterator Matrix<T>::begin()
{
    return mat.begin();
}

template <class T> typename std::deque<T>::iterator Matrix<T>::end()
{
    return mat.end();
}


template <class T> T Matrix<T>::operator()(int x, int y) const
{
    if(x<sizex && y<sizey && x>=0 && y>=0)
        return mat[x*sizey+y];
    std::cout<<"invalid indices {"<<x<<", "<<y<<"}"<<std::endl;
    return T();
}

template <class T> T& Matrix<T>::operator()(int x, int y)
{
    if(x<sizex && y<sizey && x>=0 && y>=0)
        return mat[x*sizey+y];
    std::cout<<"invalid indices {"<<x<<", "<<y<<"}"<<std::endl;
    return mat[0];
}


// submatrix
template <class T> Matrix<T> Matrix<T>::getCol(int a) const
{
    return getSubmat(0, sizex, a, a+1);
}

template <class T> Matrix<T> Matrix<T>::getCols(int a, int b) const
{
    return getSubmat(0, sizex, a, b);
}

template <class T> Matrix<T> Matrix<T>::getRow(int a) const
{
    return getSubmat(a, a+1, 0, sizey);
}

template <class T> Matrix<T> Matrix<T>::getRows(int a, int b) const
{
    return getSubmat(a, b, 0, sizey);
}

template <class T> Matrix<T> Matrix<T>::getSubmat(int a, int b, int c, int d) const
{
    if(a>=0 && a<b && b<=sizex && c>=0 && c<d && d<=sizey)
    {
        Matrix<T> submat(b-a, d-c);
        for (int i=0;i<b-a;++i)
            std::copy(cbegin()+(i+a)*sizey+c,
                      cbegin()+(i+a)*sizey+d,
                      submat.begin()+i*(d-c));
        return submat;
    }
    std::cout<<"invalid indices {"<<a<<", "<<b<<"}"<<" {"<<c<<", "<<d<<"}"<<std::endl;
    return *this;
}

template <class T> void Matrix<T>::setCol(int a, const Matrix<T>& M)
{
    setSubmat(0, a, M);
}

template <class T> void Matrix<T>::setCols(int a, const Matrix<T>& M)
{
    setSubmat(0, a, M);
}

template <class T> void Matrix<T>::setRow(int a, const Matrix<T>& M)
{
    setSubmat(a, 0, M);
}

template <class T> void Matrix<T>::setRows(int a, const Matrix<T>& M)
{
    setSubmat(a, 0, M);
}

template <class T> void Matrix<T>::setSubmat(int a, int b, const Matrix<T>& M)
{
    if(a>=0 && b>=0 && M.rowNb()+a<=sizex && M.colNb()+b<=sizey)
    {
        for (int i=0;i<M.rowNb();++i)
        {
            std::copy(M.cbegin()+i*M.colNb(),
                      M.cbegin()+i*M.colNb()+M.colNb(),
                      begin()+(i+a)*sizey+b);
        }
    }
    else
        std::cout<<"invalid indices {"<<a<<", "<<b<<"}"<<std::endl;
}


// size
template <class T> int Matrix<T>::colNb() const
{
    return sizey;
}

template <class T> int Matrix<T>::rowNb() const
{
    return sizex;
}

template <class T> int Matrix<T>::size() const
{
    return mat.size();
}


// matrix transformation
template <class T> void Matrix<T>::fliplr()
{
    for(int i=0;i<sizex;++i)
        std::reverse(begin()+i*sizey, begin()+i*sizey+sizey);
}

template <class T> void Matrix<T>::flipud()
{
    for(int i=0;i<sizex/2;++i)
    {
        std::swap_ranges(begin()+i*sizey,
                         begin()+i*sizey+sizey,
                         begin()+(sizex-1-i)*sizey);
    }
}

template <class T> void Matrix<T>::reshape(int a, int b)
{
    if(a*b==sizex*sizey && a>0)
    {
        sizex = a;
        sizey = b;
    }
    else
        std::cout<<"invalid matrix size {"<<a<<", "<<b<<"}"<<std::endl;
}

template <class T> void Matrix<T>::rot90()
{
    Matrix<T> rot_mat(sizey, sizex);
    for(int i=0;i<sizex;++i)
    {
        for(int j=0;j<sizey;++j)
            rot_mat(sizey-1-j, i) = mat[i*sizey+j];
    }
    *this = rot_mat;
}

template <class T> void Matrix<T>::rot180()
{
    std::reverse(begin(), end());
}

template <class T> void Matrix<T>::rot270()
{
    Matrix<T> rot_mat(sizey, sizex);
    for(int i=0;i<sizex;++i)
    {
        for(int j=0;j<sizey;++j)
            rot_mat(j, sizex-1-i) = mat[i*sizey+j];
    }
    *this = rot_mat;
}

template <class T> void Matrix<T>::swapcol(int a, int b)
{
    if(a>=0 && b>=0 && a<sizey && b<sizey)
    {
        for(int i=0;i<sizex;++i)
            std::swap(mat[i*sizey+a], mat[i*sizey+b]);
    }
    else
        std::cout<<"invalid column number"<<std::endl;
}

template <class T> void Matrix<T>::swaprow(int a, int b)
{
    if(a>=0 && b>=0 && a<sizex && b<sizex)
    {
        std::swap_ranges(begin()+a*sizey,
                         begin()+a*sizey+sizey,
                         begin()+b*sizey);
    }
    else
        std::cout<<"invalid row number"<<std::endl;
}

template <class T> void Matrix<T>::transpose()
{
    Matrix<T> t_mat(sizey, sizex);
    for(int i=0;i<sizex;++i)
    {
        for(int j=0;j<sizey;++j)
            t_mat(j, i) = mat[i*sizey+j];
    }
    *this = t_mat;
}


// add/remove a row/column
template <class T> void Matrix<T>::delCol(int c)
{
    if(c<sizey && c>=0 && sizey>1)
    {
        for(int i=0;i<sizex;++i)
            mat.erase(mat.begin()+i*sizey+c-i);
        --sizey;
    }
    else if(sizey<=1)
        std::cout<<"the matrix has only one column"<<std::endl;
    else
        std::cout<<"invalid column number"<<std::endl;
}

template <class T> void Matrix<T>::delRow(int r)
{
    if(r<sizex && r>=0 && sizex>1)
    {
        mat.erase(begin()+r*sizey, begin()+r*sizey+sizey);
        --sizex;
    }
    else if(sizex<=1)
        std::cout<<"the matrix has only one row"<<std::endl;
    else
        std::cout<<"invalid row number"<<std::endl;
}

template <class T> void Matrix<T>::newCol()
{
    newCol(sizey);
}

template <class T> void Matrix<T>::newCol(int c)
{
    if(c<=sizey && c>=0)
    {
        for(int i=0;i<sizex;++i)
            mat.insert(mat.begin()+i*sizey+c+i, 0);
        ++sizey;
    }
    else
        std::cout<<"invalid column number"<<std::endl;
}

template <class T> void Matrix<T>::newRow()
{
    newRow(sizex);
}

template <class T> void Matrix<T>::newRow(int r)
{
    if(r<=sizex && r>=0)
    {
        mat.insert(mat.begin()+r*sizey, sizey,0);
        ++sizex;
    }
    else
        std::cout<<"invalid row number"<<std::endl;
}


// operators
template <class T> void Matrix<T>::operator+=(const Matrix<T>& M)
{
    if(sizex==M.rowNb() && sizey==M.colNb())
        std::transform(cbegin(), cend(), M.cbegin(), begin(), std::plus<T>());
    else
        std::cout<<"dimension mismatch"<<std::endl;
}

template <class T> void Matrix<T>::operator-=(const Matrix<T>& M)
{
    if(sizex==M.rowNb() && sizey==M.colNb())
        std::transform(cbegin(), cend(), M.cbegin(), begin(), std::minus<T>());
    else
        std::cout<<"dimension mismatch"<<std::endl;
}

template <class T> void Matrix<T>::operator*=(const Matrix<T>& M)
{
    if(sizex==M.rowNb() && sizey==M.colNb())
        std::transform(cbegin(), cend(), M.cbegin(), begin(), std::multiplies<T>());
    else
        std::cout<<"dimension mismatch"<<std::endl;
}

template <class T> void Matrix<T>::operator+=(const T& value)
{
    std::for_each(begin(), end(), [value](T& x){x+=value;});
}

template <class T> void Matrix<T>::operator-=(const T& value)
{
    std::for_each(begin(), end(), [value](T& x){x-=value;});
}

template <class T> void Matrix<T>::operator*=(const T& value)
{
    std::for_each(begin(), end(), [value](T& x){x*=value;});
}

template <class T> void Matrix<T>::operator/=(const T& value)
{
    if(std::abs(value)>10e-9)
        std::for_each(begin(), end(), [value](T& x){x/=value;});
    else
        std::cout<<"division by zero"<<std::endl;
}

template <class T> void Matrix<T>::operator%=(const T& value)
{
    if(std::abs(value)>10e-9)
        std::for_each(begin(), end(), [value](T& x){x%=value;});
    else
        std::cout<<"division by zero"<<std::endl;
}

template <class T> Matrix<T> Matrix<T>::operator+(const Matrix<T>& M) const
{
    Matrix<T> newM(sizex, sizey);
    if(sizex==M.rowNb() && sizey==M.colNb())
        std::transform(cbegin(), cend(), M.cbegin(), newM.begin(), std::plus<T>());
    else
        std::cout<<"dimension mismatch"<<std::endl;
    return newM;
}

template <class T> Matrix<T> Matrix<T>::operator-(const Matrix<T>& M) const
{
    Matrix<T> newM(sizex, sizey);
    if(sizex==M.rowNb() && sizey==M.colNb())
        std::transform(cbegin(), cend(), M.cbegin(), newM.begin(), std::minus<T>());
    else
        std::cout<<"dimension mismatch"<<std::endl;
    return newM;
}

template <class T> Matrix<T> Matrix<T>::operator*(const Matrix<T>& M) const
{
    Matrix<T> newM(sizex, sizey);
    if(sizex==M.rowNb() && sizey==M.colNb())
        std::transform(cbegin(), cend(), M.cbegin(), newM.begin(), std::multiplies<T>());
    else
        std::cout<<"dimension mismatch"<<std::endl;
    return newM;
}

template <class T> Matrix<T> Matrix<T>::operator+(const T& value) const
{
    Matrix<T> newM(sizex, sizey);
    std::transform(cbegin(), cend(), std::deque<T>(size(), value).begin(), newM.begin(), std::plus<T>());
    return newM;
}

template <class T> Matrix<T> Matrix<T>::operator-(const T& value) const
{
    Matrix<T> newM(sizex, sizey);
    std::transform(cbegin(), cend(), std::deque<T>(size(), value).begin(), newM.begin(), std::minus<T>());
    return newM;
}

template <class T> Matrix<T> Matrix<T>::operator*(const T& value) const
{
    Matrix<T> newM(sizex, sizey);
    std::transform(cbegin(), cend(), std::deque<T>(size(), value).begin(), newM.begin(), std::multiplies<T>());
    return newM;
}

template <class T> Matrix<T> Matrix<T>::operator/(const T& value) const
{
    Matrix<T> newM(sizex, sizey);
    if(std::abs(value)>10e-9)
        std::transform(cbegin(), cend(), std::deque<T>(size(), value).begin(), newM.begin(), std::divides<T>());
    else
        std::cout<<"division by zero"<<std::endl;
    return newM;
}

template <class T> Matrix<T> Matrix<T>::operator%(const T& value) const
{
    Matrix<T> newM(sizex, sizey);
    if(std::abs(value)>10e-9)
        std::transform(cbegin(), cend(), std::deque<T>(size(), value).begin(), newM.begin(), std::modulus<T>());
    else
        std::cout<<"division by zero"<<std::endl;
    return newM;
}

template <class T> Matrix<T> Matrix<T>::operator-() const
{
    Matrix<T> newM(sizex, sizey);
    std::transform(cbegin(), cend(), newM.begin(), std::negate<T>());
    return newM;
}

template <class T> Matrix<bool> Matrix<T>::operator==(const Matrix<T>& M) const
{
    Matrix<bool> newM(sizex, sizey);
    if(sizex==M.rowNb() && sizey==M.colNb())
        std::transform(cbegin(), cend(), M.cbegin(), newM.begin(), std::equal_to<T>());
    else
        std::cout<<"dimension mismatch"<<std::endl;
    return newM;
}

template <class T> Matrix<bool> Matrix<T>::operator!=(const Matrix<T>& M) const
{
    Matrix<bool> newM(sizex, sizey);
    if(sizex==M.rowNb() && sizey==M.colNb())
        std::transform(cbegin(), cend(), M.cbegin(), newM.begin(), std::not_equal_to<T>());
    else
        std::cout<<"dimension mismatch"<<std::endl;
    return newM;
}

template <class T> Matrix<bool> Matrix<T>::operator>=(const Matrix<T>& M) const
{
    Matrix<bool> newM(sizex, sizey);
    if(sizex==M.rowNb() && sizey==M.colNb())
        std::transform(cbegin(), cend(), M.cbegin(), newM.begin(), std::greater_equal<T>());
    else
        std::cout<<"dimension mismatch"<<std::endl;
    return newM;
}

template <class T> Matrix<bool> Matrix<T>::operator<=(const Matrix<T>& M) const
{
    Matrix<bool> newM(sizex, sizey);
    if(sizex==M.rowNb() && sizey==M.colNb())
        std::transform(cbegin(), cend(), M.cbegin(), newM.begin(), std::less_equal<T>());
    else
        std::cout<<"dimension mismatch"<<std::endl;
    return newM;
}

template <class T> Matrix<bool> Matrix<T>::operator>(const Matrix<T>& M) const
{
    Matrix<bool> newM(sizex, sizey);
    if(sizex==M.rowNb() && sizey==M.colNb())
        std::transform(cbegin(), cend(), M.cbegin(), newM.begin(), std::greater<T>());
    else
        std::cout<<"dimension mismatch"<<std::endl;
    return newM;
}

template <class T> Matrix<bool> Matrix<T>::operator<(const Matrix<T>& M) const
{
    Matrix<bool> newM(sizex, sizey);
    if(sizex==M.rowNb() && sizey==M.colNb())
        std::transform(cbegin(), cend(), M.cbegin(), newM.begin(), std::less<T>());
    else
        std::cout<<"dimension mismatch"<<std::endl;
    return newM;
}

template <class T> Matrix<bool> Matrix<T>::operator==(const T& value) const
{
    Matrix<bool> newM(sizex, sizey);
    std::transform(cbegin(), cend(), std::deque<T>(size(), value).begin(), newM.begin(), std::equal_to<T>());
    return newM;
}

template <class T> Matrix<bool> Matrix<T>::operator!=(const T& value) const
{
    Matrix<bool> newM(sizex, sizey);
    std::transform(cbegin(), cend(), std::deque<T>(size(), value).begin(), newM.begin(), std::not_equal_to<T>());
    return newM;
}

template <class T> Matrix<bool> Matrix<T>::operator>=(const T& value) const
{
    Matrix<bool> newM(sizex, sizey);
    std::transform(cbegin(), cend(), std::deque<T>(size(), value).cbegin(), newM.begin(), std::greater_equal<T>());
    return newM;
}

template <class T> Matrix<bool> Matrix<T>::operator<=(const T& value) const
{
    Matrix<bool> newM(sizex, sizey);
    std::transform(cbegin(), cend(), std::deque<T>(size(), value).cbegin(), newM.begin(), std::less_equal<T>());
    return newM;
}

template <class T> Matrix<bool> Matrix<T>::operator>(const T& value) const
{
    Matrix<bool> newM(sizex, sizey);
    std::transform(cbegin(), cend(), std::deque<T>(size(), value).cbegin(), newM.begin(), std::greater<T>());
    return newM;
}

template <class T> Matrix<bool> Matrix<T>::operator<(const T& value) const
{
    Matrix<bool> newM(sizex, sizey);
    std::transform(cbegin(), cend(), std::deque<T>(size(), value).cbegin(), newM.begin(), std::less<T>());
    return newM;
}



// specific constructors
template <class T> Matrix<T> arange(T a, T b, T step)
{
    if(step == T(0))
    {
        std::cout<<"invalid step argument"<<std::endl;
        return zeros<T>(1);
    }
    int row_nb = (int)((b-a)/step);
    if(row_nb<1)
    {
        std::cout<<"unable to create the matrix"<<std::endl;
        return zeros<T>(1);
    }
    Matrix<T> new_mat(row_nb, 1);
    std::iota(new_mat.begin(), new_mat.end(), T(0));
    return a+step*new_mat;
}

template <class T> Matrix<T> full(int a, T value)
{
    return full<T>(a, a, value);
}

template <class T> Matrix<T> full(int a, int b, T value)
{
    Matrix<T> new_mat(a, b);
    std::fill(new_mat.begin(), new_mat.end(), value);
    return new_mat;
}

template <class T> Matrix<T> id(int a)
{
    return id<T>(a, a);
}

template <class T> Matrix<T> id(int a, int b)
{
    Matrix<T> new_mat = zeros<T>(a, b);
    int minimum = a<b?a:b;
    for(int i=0;i<minimum;++i)
        new_mat(i, i) = T(1);
    return new_mat;
}


template <class T> std::pair<Matrix<T>, Matrix<T> > meshgrid(int a)
{
    return meshgrid<T>(a, a);
}

template <class T> std::pair<Matrix<T>, Matrix<T> > meshgrid(int a, int b)
{
    Matrix<T> X(a, b);
    Matrix<T> Y(a, b);
    for(int i=0;i<a;++i)
        std::fill(X.begin()+i*b, X.begin()+(i+1)*b, T(i));
    for(int i=0;i<a;++i)
        std::iota(Y.begin()+i*b, Y.begin()+(i+1)*b, T(0));
    return std::make_pair(X, Y);
}

template <class T> Matrix<T> rand(int a)
{
    return rand<T>(a, a);
}

template <class T> Matrix<T> rand(int a, int b)
{
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_int_distribution<T> distribution(0);
    Matrix<T> new_mat(a, b);
    auto gen = std::bind(distribution, generator);
    std::generate(new_mat.begin(), new_mat.end(), gen);
    return new_mat;
}


template <class T> Matrix<T> randn(int a, T mean, T stddev)
{
    return randn<T>(a, a, mean, stddev);
}

template <class T> Matrix<T> randn(int a, int b, T mean, T stddev)
{
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::normal_distribution<T> distribution(mean, stddev);
    Matrix<T> new_mat(a, b);
    auto gen = std::bind(distribution, generator);
    generate(new_mat.begin(), new_mat.end(), gen);
    return new_mat;
}

template <class T> Matrix<T> ones(int a)
{
    return ones<T>(a, a);
}

template <class T> Matrix<T> ones(int a, int b)
{
    Matrix<T> new_mat(a, b);
    std::fill(new_mat.begin(), new_mat.end(), T(1));
    return new_mat;
}

template <class T> Matrix<T> zeros(int a)
{
    return zeros<T>(a, a);
}

template <class T> Matrix<T> zeros(int a, int b)
{
    Matrix<T> new_mat(a, b);
    std::fill(new_mat.begin(), new_mat.end(), T(0));
    return new_mat;
}

// functions
template <class T, typename Type> Matrix<Type> apply(const Matrix<T>& M, Type (*ptr)(T))
{
    Matrix<Type> new_mat(M.rowNb(), M.colNb());
    std::transform(M.cbegin(), M.cend(), new_mat.begin(), *ptr);
    return new_mat;
}

template <class T, typename Type> Matrix<Type> apply(const Matrix<T>& M, Type (T::*ptr)() const)
{
    Matrix<Type> new_mat(M.rowNb(), M.colNb());
    std::transform(M.cbegin(), M.cend(), new_mat.begin(), [ptr](T x){return (x.*ptr)();});
    return new_mat;
}

// construct a matrix depending on a condition
template <class T> Matrix<T> where(const Matrix<bool>& cdt, const T& a, const T& b)
{
    Matrix<T> my_matrix(cdt.rowNb(), cdt.colNb());
    std::transform(cdt.cbegin(), cdt.cend(), my_matrix.begin(), [a, b](bool cdt){return cdt?a:b;});
    return my_matrix;
}

template <class T> Matrix<T> where(const Matrix<bool>& cdt, const Matrix<T>& A, const Matrix<T>& B)
{
    if(A.rowNb()==B.rowNb() && A.colNb()==B.colNb() && cdt.rowNb()==B.rowNb() && cdt.colNb()==B.colNb())
    {
        Matrix<T> my_matrix(cdt.rowNb(), cdt.colNb());
        auto A_it = A.cbegin();
        auto B_it = B.cbegin();
        auto result_it = my_matrix.begin();
        auto cdt_it = cdt.cbegin();
        auto it_end = cdt.cend();
        while(cdt_it != it_end)
        {
            *result_it = (*cdt_it)?(*A_it):(*B_it);
            ++result_it; ++cdt_it; ++A_it; ++B_it;
        }
        return my_matrix;
    }
    std::cout<<"dimension mismatch"<<std::endl;
    return A;
}

template <class T> Matrix<T> compress(const Matrix<bool>& cdt, const Matrix<T>& M, int axis)
{
    if(axis==1)
    {
        if(cdt.rowNb()==M.rowNb() && cdt.colNb()==1)
        {
            Matrix<T> new_matrix = zeros<T>(1, M.colNb());
            for(int i=0;i<cdt.rowNb();++i)
            {
                if(cdt(i, 0))
                {
                    new_matrix.newRow();
                    new_matrix.setRow(new_matrix.rowNb()-1, M.getRow(i));
                }
            }
            new_matrix.delRow(0);
            return new_matrix;
        }
        std::cout<<"dimension mismatch"<<std::endl;
    }
    else if(axis==2)
    {
        if(cdt.colNb()==M.colNb() && cdt.rowNb()==1)
        {
            Matrix<T> new_matrix = zeros<T>(M.rowNb(), 1);
            for(int i=0;i<cdt.colNb();++i)
            {
                if(cdt(0, i))
                {
                    new_matrix.newCol();
                    new_matrix.setCol(new_matrix.colNb()-1, M.getCol(i));
                }
            }
            new_matrix.delCol(0);
            return new_matrix;
        }
        std::cout<<"dimension mismatch"<<std::endl;
    }
    else
        std::cout<<"invalid axis"<<std::endl;
    return zeros<T>(1);
}

template <class T> int count(const Matrix<T>& M, const T& value)
{
    return std::count(M.cbegin(), M.cend(), value);
}

template <class T> int count_nonzero(const Matrix<T>& M)
{
    return M.size()-std::count(M.cbegin(), M.cend(), T(0));
}

template <class T> void replace(Matrix<T>& M, const T& a, const T& b)
{
    std::replace(M.begin(), M.end(), a, b);
}

template <class T> void replace_if(Matrix<T>& M, const Matrix<bool>& cdt, const T& value)
{
    if(cdt.rowNb()==M.rowNb() && cdt.colNb()==M.colNb())
    {
        auto it = cdt.cbegin();
        std::replace_if(M.begin(), M.end(), [&it](T){return *it++;}, value);
    }
    else
        std::cout<<"dimension mismatch"<<std::endl;

}

template <class T> Matrix<int> argwhere(const Matrix<bool>& M)
{
    std::vector<int> indices(M.size());
    std::iota(indices.begin(), indices.end(), 0);
    auto it_begin = M.cbegin();
    auto it = std::remove_if(indices.begin(), indices.end(), [&it_begin](int){return !(*it_begin++);});
    indices.resize(std::distance(indices.begin(), it));
    int w = M.colNb();
    Matrix<int> temp(indices.size(), 1);
    std::move(indices.begin(), indices.end(), temp.begin());
    Matrix<int> new_mat(indices.size(), 2);
    new_mat.setCol(0, temp/w);
    new_mat.setCol(1, temp%w);
    return new_mat;
}

template <class T> Matrix<T> unique(const Matrix<T>& M)
{
    std::vector<T> vec(M.cbegin(), M.cend());
    std::sort(vec.begin(), vec.end());
    auto it = std::unique(vec.begin(), vec.end());
    int n = std::distance(vec.begin(), it);
    vec.resize(n);
    Matrix<T> new_mat(n, 1);
    std::move(vec.begin(), vec.end(), new_mat.begin());
    return new_mat;
}


template <class T> Matrix<T> diag(const Matrix<T>& M)
{
    int row_nb = std::min(M.rowNb(), M.colNb());
    Matrix<T> new_mat(row_nb, 1);
    for(int i=0;i<row_nb;++i)
        new_mat(i, 0) = M(i, i);
    return new_mat;
}

template <class T> Matrix<T> transpose(const Matrix<T>& M)
{
    Matrix<T> t_mat(M.colNb(), M.rowNb());
    for(int i=0;i<M.rowNb();++i)
    {
        for(int j=0;j<M.colNb();++j)
            t_mat(j, i) = M(i, j);
    }
    return t_mat;
}

// operators
template <class T> Matrix<T> operator+(const T& value, const Matrix<T>& M)
{
    return M+value;
}

template <class T> Matrix<T> operator-(const T& value, const Matrix<T>& M)
{
    return -M+value;
}

template <class T> Matrix<T> operator*(const T& value, const Matrix<T>& M)
{
    return M*value;
}

template <class T> Matrix<T> operator/(const T& value, const Matrix<T>& M)
{
    Matrix<T> new_mat = M;
    std::for_each(new_mat.begin(), new_mat.end(), [value](T& x){x=value/x;});
    return new_mat;
}

template <class T> std::ostream& operator<<(std::ostream& s, const Matrix<T>& my_matrix)
{
    for(int i=0;i<my_matrix.rowNb();++i)
    {
        for(int j=0;j<my_matrix.colNb();++j)
            s<<my_matrix(i, j)<<'\t';
        s<<std::endl;
    }
    return s;
}
