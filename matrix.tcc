#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <numeric>
#include <random>


// constructors & accessors
template <class T> Matrix<T>::Matrix() : Matrix<T>::Matrix(1)
{
}

template <class T> Matrix<T>::Matrix(std::size_t a) : Matrix<T>::Matrix(a, a)
{
}

template <class T> Matrix<T>::Matrix(std::size_t a, std::size_t b)
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


template <class T> T Matrix<T>::operator()(std::size_t x, std::size_t y) const
{
    if(x<sizex && y<sizey)
        return mat[x*sizey+y];
    std::cout<<"invalid indices {"<<x<<", "<<y<<"}"<<std::endl;
    return T();
}

template <class T> T& Matrix<T>::operator()(std::size_t x, std::size_t y)
{
    if(x<sizex && y<sizey)
        return mat[x*sizey+y];
    std::cout<<"invalid indices {"<<x<<", "<<y<<"}"<<std::endl;
    return mat[0];
}


// submatrix
template <class T> Matrix<T> Matrix<T>::getCol(std::size_t a) const
{
    return getSubmat(0, sizex, a, a+1);
}

template <class T> Matrix<T> Matrix<T>::getCols(std::size_t a, std::size_t b) const
{
    return getSubmat(0, sizex, a, b);
}

template <class T> Matrix<T> Matrix<T>::getDiag() const
{
    std::size_t row_nb = std::min(sizex, sizey);
    Matrix<T> new_mat(row_nb, 1);
    auto it_this = cbegin();
    auto it_new_mat =  new_mat.begin();
    for(std::size_t i=0;i<row_nb;++i)
    {
        *it_new_mat = *it_this;
        it_this+=sizey+1;
        ++it_new_mat;
    }
    return new_mat;
}

template <class T> Matrix<T> Matrix<T>::getRow(std::size_t a) const
{
    return getSubmat(a, a+1, 0, sizey);
}

template <class T> Matrix<T> Matrix<T>::getRows(std::size_t a, std::size_t b) const
{
    return getSubmat(a, b, 0, sizey);
}

template <class T> Matrix<T> Matrix<T>::getSubmat(std::size_t a, std::size_t b, std::size_t c, std::size_t d) const
{
    if(a<b && b<=sizex && c<d && d<=sizey)
    {
        Matrix<T> submat(b-a, d-c);
        for (std::size_t i=0;i<b-a;++i)
            std::copy(cbegin()+(i+a)*sizey+c,
                      cbegin()+(i+a)*sizey+d,
                      submat.begin()+i*(d-c));
        return submat;
    }
    std::cout<<"invalid indices {"<<a<<", "<<b<<"}"<<" {"<<c<<", "<<d<<"}"<<std::endl;
    return *this;
}

template <class T> void Matrix<T>::setCol(std::size_t a, const Matrix<T>& M)
{
    setSubmat(0, a, M);
}

template <class T> void Matrix<T>::setCols(std::size_t a, const Matrix<T>& M)
{
    setSubmat(0, a, M);
}

template <class T> void Matrix<T>::setDiag(const Matrix<T>& M)
{
    std::size_t row_nb = std::min(sizex, sizey);
    if(row_nb==M.size() && (M.rowNb()==1 || M.colNb()==1))
    {
        auto it_this = begin();
        auto it_m =  M.cbegin();
        for(std::size_t i=0;i<row_nb;++i)
        {
            *it_this = *it_m;
            it_this+=sizey+1;
            ++it_m;
        }
    }
    else
        std::cout<<"dimension mismatch"<<std::endl;
}

template <class T> void Matrix<T>::setRow(std::size_t a, const Matrix<T>& M)
{
    setSubmat(a, 0, M);
}

template <class T> void Matrix<T>::setRows(std::size_t a, const Matrix<T>& M)
{
    setSubmat(a, 0, M);
}

template <class T> void Matrix<T>::setSubmat(std::size_t a, std::size_t b, const Matrix<T>& M)
{
    if(M.rowNb()+a<=sizex && M.colNb()+b<=sizey)
    {
        for (std::size_t i=0;i<M.rowNb();++i)
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
template <class T> std::size_t Matrix<T>::colNb() const
{
    return sizey;
}

template <class T> std::size_t Matrix<T>::rowNb() const
{
    return sizex;
}

template <class T> std::size_t Matrix<T>::size() const
{
    return mat.size();
}


// matrix transformation
template <class T> void Matrix<T>::fliplr()
{
    for(std::size_t i=0;i<sizex;++i)
        std::reverse(begin()+i*sizey, begin()+i*sizey+sizey);
}

template <class T> void Matrix<T>::flipud()
{
    for(std::size_t i=0;i<sizex/2;++i)
    {
        std::swap_ranges(begin()+i*sizey,
                         begin()+i*sizey+sizey,
                         begin()+(sizex-1-i)*sizey);
    }
}

template <class T> void Matrix<T>::reshape(std::size_t a, std::size_t b)
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
    for(std::size_t i=0;i<sizex;++i)
    {
        for(std::size_t j=0;j<sizey;++j)
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
    for(std::size_t i=0;i<sizex;++i)
    {
        for(std::size_t j=0;j<sizey;++j)
            rot_mat(j, sizex-1-i) = mat[i*sizey+j];
    }
    *this = rot_mat;
}

template <class T> void Matrix<T>::swapcol(std::size_t a, std::size_t b)
{
    if(a<sizey && b<sizey)
    {
        for(std::size_t i=0;i<sizex;++i)
            std::swap(mat[i*sizey+a], mat[i*sizey+b]);
    }
    else
        std::cout<<"invalid column number"<<std::endl;
}

template <class T> void Matrix<T>::swaprow(std::size_t a, std::size_t b)
{
    if(a<sizex && b<sizex)
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
    if(sizex==1 || sizey==1)
    {
        std::swap(sizex, sizey);
    }
    else if(sizex==sizey)
    {
        for(std::size_t i=0;i<sizex-1;++i)
        {
            for(std::size_t j=i+1;j<sizey;++j)
                std::swap(mat[i*sizey+j], mat[j*sizey+i]);
        }
    }
    else
    {
        Matrix<T> t_mat(sizey, sizex);
        for(std::size_t i=0;i<sizex;++i)
        {
            for(std::size_t j=0;j<sizey;++j)
                t_mat(j, i) = mat[i*sizey+j];
        }
        *this = t_mat;
    }
}


// add/remove a row/column
template <class T> void Matrix<T>::delCol(std::size_t c)
{
    if(c<sizey && sizey>1)
    {
        for(std::size_t i=0;i<sizex;++i)
            mat.erase(mat.begin()+i*sizey+c-i);
        --sizey;
    }
    else if(sizey<=1)
        std::cout<<"the matrix has only one column"<<std::endl;
    else
        std::cout<<"invalid column number"<<std::endl;
}

template <class T> void Matrix<T>::delRow(std::size_t r)
{
    if(r<sizex && sizex>1)
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

template <class T> void Matrix<T>::newCol(std::size_t c)
{
    if(c<=sizey)
    {
        for(std::size_t i=0;i<sizex;++i)
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

template <class T> void Matrix<T>::newRow(std::size_t r)
{
    if(r<=sizex)
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
    std::transform(cbegin(), cend(), newM.begin(), [value](const T& v){return v+value;});
    return newM;
}

template <class T> Matrix<T> Matrix<T>::operator-(const T& value) const
{
    Matrix<T> newM(sizex, sizey);
    std::transform(cbegin(), cend(), newM.begin(), [value](const T& v){return v-value;});
    return newM;
}

template <class T> Matrix<T> Matrix<T>::operator*(const T& value) const
{
    Matrix<T> newM(sizex, sizey);
    std::transform(cbegin(), cend(), newM.begin(), [value](const T& v){return v*value;});
    return newM;
}

template <class T> Matrix<T> Matrix<T>::operator/(const T& value) const
{
    Matrix<T> newM(sizex, sizey);
    if(std::abs(value)>10e-9)
        std::transform(cbegin(), cend(), newM.begin(), [value](const T& v){return v/value;});
    else
        std::cout<<"division by zero"<<std::endl;
    return newM;
}

template <class T> Matrix<T> Matrix<T>::operator%(const T& value) const
{
    Matrix<T> newM(sizex, sizey);
    if(std::abs(value)>10e-9)
        std::transform(cbegin(), cend(), newM.begin(), [value](const T& v){return v%value;});
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
    std::transform(cbegin(), cend(), newM.begin(), [value](const T& v){return v==value;});
    return newM;
}

template <class T> Matrix<bool> Matrix<T>::operator!=(const T& value) const
{
    Matrix<bool> newM(sizex, sizey);
    std::transform(cbegin(), cend(), newM.begin(), [value](const T& v){return v!=value;});
    return newM;
}

template <class T> Matrix<bool> Matrix<T>::operator>=(const T& value) const
{
    Matrix<bool> newM(sizex, sizey);
    std::transform(cbegin(), cend(), newM.begin(), [value](const T& v){return v>=value;});
    return newM;
}

template <class T> Matrix<bool> Matrix<T>::operator<=(const T& value) const
{
    Matrix<bool> newM(sizex, sizey);
    std::transform(cbegin(), cend(), newM.begin(), [value](const T& v){return v<=value;});
    return newM;
}

template <class T> Matrix<bool> Matrix<T>::operator>(const T& value) const
{
    Matrix<bool> newM(sizex, sizey);
    std::transform(cbegin(), cend(), newM.begin(), [value](const T& v){return v>value;});
    return newM;
}

template <class T> Matrix<bool> Matrix<T>::operator<(const T& value) const
{
    Matrix<bool> newM(sizex, sizey);
    std::transform(cbegin(), cend(), newM.begin(), [value](const T& v){return v<value;});
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
    std::size_t row_nb = (std::size_t)((b-a)/step);
    if(row_nb<1)
    {
        std::cout<<"unable to create the matrix"<<std::endl;
        return zeros<T>(1);
    }
    Matrix<T> new_mat(row_nb, 1);
    std::iota(new_mat.begin(), new_mat.end(), T(0));
    return a+step*new_mat;
}

template <class T> Matrix<T> full(std::size_t a, T value)
{
    return full<T>(a, a, value);
}

template <class T> Matrix<T> full(std::size_t a, std::size_t b, T value)
{
    Matrix<T> new_mat(a, b);
    std::fill(new_mat.begin(), new_mat.end(), value);
    return new_mat;
}

template <class T> Matrix<T> id(std::size_t a)
{
    return id<T>(a, a);
}

template <class T> Matrix<T> id(std::size_t a, std::size_t b)
{
    Matrix<T> new_mat = zeros<T>(a, b);
    std::size_t minimum = a<b?a:b;
    for(std::size_t i=0;i<minimum;++i)
        new_mat(i, i) = T(1);
    return new_mat;
}


template <class T> std::pair<Matrix<T>, Matrix<T> > meshgrid(std::size_t a)
{
    return meshgrid<T>(a, a);
}

template <class T> std::pair<Matrix<T>, Matrix<T> > meshgrid(std::size_t a, std::size_t b)
{
    Matrix<T> X(a, b);
    Matrix<T> Y(a, b);
    for(std::size_t i=0;i<a;++i)
    {
        std::fill(X.begin()+i*b, X.begin()+(i+1)*b, T(i));
        std::iota(Y.begin()+i*b, Y.begin()+(i+1)*b, T(0));
    }
    return std::make_pair(X, Y);
}

template <class T> Matrix<T> rand(std::size_t a)
{
    return rand<T>(a, a);
}

template <class T> Matrix<T> rand(std::size_t a, std::size_t b)
{
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_int_distribution<T> distribution(0);
    Matrix<T> new_mat(a, b);
    auto gen = std::bind(distribution, generator);
    std::generate(new_mat.begin(), new_mat.end(), gen);
    return new_mat;
}


template <class T> Matrix<T> randn(std::size_t a, T mean, T stddev)
{
    return randn<T>(a, a, mean, stddev);
}

template <class T> Matrix<T> randn(std::size_t a, std::size_t b, T mean, T stddev)
{
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::normal_distribution<T> distribution(mean, stddev);
    Matrix<T> new_mat(a, b);
    auto gen = std::bind(distribution, generator);
    generate(new_mat.begin(), new_mat.end(), gen);
    return new_mat;
}

template <class T> Matrix<T> ones(std::size_t a)
{
    return ones<T>(a, a);
}

template <class T> Matrix<T> ones(std::size_t a, std::size_t b)
{
    Matrix<T> new_mat(a, b);
    std::fill(new_mat.begin(), new_mat.end(), T(1));
    return new_mat;
}

template <class T> Matrix<T> zeros(std::size_t a)
{
    return zeros<T>(a, a);
}

template <class T> Matrix<T> zeros(std::size_t a, std::size_t b)
{
    Matrix<T> new_mat(a, b);
    std::fill(new_mat.begin(), new_mat.end(), T(0));
    return new_mat;
}

// functions
template <typename Type, class T> Matrix<T> apply(const Matrix<Type>& M, T (*ptr)(Type))
{
    Matrix<T> new_mat(M.rowNb(), M.colNb());
    std::transform(M.cbegin(), M.cend(), new_mat.begin(), *ptr);
    return new_mat;
}


template <typename Type, class T> Matrix<T> apply(const Matrix<Type>& M, T (*ptr)(const Type&))
{
    Matrix<T> new_mat(M.rowNb(), M.colNb());
    std::transform(M.cbegin(), M.cend(), new_mat.begin(), *ptr);
    return new_mat;
}

template <typename Type, class T> Matrix<T> apply(const Matrix<Type>& M, T (Type::*ptr)() const)
{
    Matrix<T> new_mat(M.rowNb(), M.colNb());
    std::transform(M.cbegin(), M.cend(), new_mat.begin(), [ptr](Type x){return (x.*ptr)();});
    return new_mat;
}

template <typename T1, typename T2, class T> Matrix<T> apply(const Matrix<T1>& M1, const Matrix<T2>& M2, T (*ptr)(T1, T2))
{
    Matrix<T> new_mat(M1.rowNb(), M1.colNb());
    if(M1.rowNb()==M2.rowNb() && M1.colNb()==M2.colNb())
    {

        std::transform(M1.cbegin(), M1.cend(), M2.cbegin(), new_mat.begin(), *ptr);
        return new_mat;
    }
    std::cout<<"dimension mismatch"<<std::endl;
    return new_mat;
}

template <typename T1, typename T2, class T> Matrix<T> apply(const Matrix<T1>& M1, const Matrix<T2>& M2, T (*ptr)(const T1&, const T2&))
{
    Matrix<T> new_mat(M1.rowNb(), M1.colNb());
    if(M1.rowNb()==M2.rowNb() && M1.colNb()==M2.colNb())
    {
        std::transform(M1.cbegin(), M1.cend(), M2.cbegin(), new_mat.begin(), *ptr);
        return new_mat;
    }
    std::cout<<"dimension mismatch"<<std::endl;
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
            for(std::size_t i=0;i<cdt.rowNb();++i)
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
            for(std::size_t i=0;i<cdt.colNb();++i)
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

template <class T> std::size_t count(const Matrix<T>& M, const T& value)
{
    return std::count(M.cbegin(), M.cend(), value);
}

template <class T> std::size_t count_nonzero(const Matrix<T>& M)
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

template <class T> Matrix<std::size_t> argwhere(const Matrix<bool>& M)
{
    std::deque<std::size_t> indices(M.size());
    std::iota(indices.begin(), indices.end(), 0);
    auto it_begin = M.cbegin();
    auto it = std::remove_if(indices.begin(), indices.end(), [&it_begin](std::size_t){return !(*it_begin++);});
    indices.resize(std::distance(indices.begin(), it));
    std::size_t w = M.colNb();
    Matrix<std::size_t> temp(indices.size(), 1);
    std::move(indices.begin(), indices.end(), temp.begin());
    Matrix<std::size_t> new_mat(indices.size(), 2);
    new_mat.setCol(0, temp/w);
    new_mat.setCol(1, temp%w);
    return new_mat;
}

template <class T> Matrix<T> unique(const Matrix<T>& M)
{
    std::deque<T> vec(M.cbegin(), M.cend());
    std::sort(vec.begin(), vec.end());
    auto it = std::unique(vec.begin(), vec.end());
    std::size_t n = std::distance(vec.begin(), it);
    vec.resize(n);
    Matrix<T> new_mat(n, 1);
    std::move(vec.begin(), vec.end(), new_mat.begin());
    return new_mat;
}


template <class T> Matrix<T> transpose(const Matrix<T>& M)
{
    std::size_t H = M.rowNb(), W = M.colNb();
    Matrix<T> t_mat;
    if(H==1 || W==1)
    {
        t_mat=M;
        t_mat.reshape(W, H);
    }
    else
    {
        t_mat = Matrix<T>(W, H);
        for(std::size_t i=0;i<H;++i)
        {
            for(std::size_t j=0;j<W;++j)
                t_mat(j, i) = M(i, j);
        }
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
    for(std::size_t i=0, I=my_matrix.rowNb();i<I;++i)
    {
        for(std::size_t j=0, J=my_matrix.colNb();j<J;++j)
            s<<my_matrix(i, j)<<'\t';
        s<<std::endl;
    }
    return s;
}
