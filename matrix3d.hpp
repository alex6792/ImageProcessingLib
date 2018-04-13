/*!
 * \author Alexandre Krebs
 * \file matrix3d.hpp
 * \brief 3D matrix class
 */


#pragma once

#ifndef MATRIX3D_HPP
#define MATRIX3D_HPP


#include <deque>
#include <functional>
#include <iostream>


/*!
 * \class Matrix
 * \brief template class representing a matrix
 */
template <class T> class Matrix3D
{

    private :
        std::size_t sizex; /*!< number of rows*/
        std::size_t sizey; /*!< number of columns*/
        std::size_t sizez; /*!< depth*/
        std::deque<T> mat; /*!< data*/

    public :
        // constructors
        /*!
         * \brief Default constructor that construct a matrix 1X1X1
         */
        Matrix3D();

        /*!
         * \brief NxNxN matrix constructor
         * \param N the size of the matrix
         */
        explicit Matrix3D(std::size_t N);

        /*!
         * \brief MxNxP matrix constructor
         * \param M number of rows
         * \param N number of columns
         * \param P depth
         */
        Matrix3D(std::size_t M, std::size_t N, size_t P);

        /*!
         * \brief cast a matrix of type T1 to a matrix of type T2
         * \param M the matrix of type T1
         */
        template <typename Type> explicit Matrix3D(const Matrix<Type>& M);

        /*!
         * \brief construct a matrix from a C-style static array
         * \param A a static array
         */
        Matrix3D(std::initializer_list<std::initializer_list<std::initializer_list<T> > > A);

        // iterators
        /*!
         * \return a const iterator on the first element
         */
        typename std::deque<T>::const_iterator cbegin() const;

        /*!
         * \return a const iterator on the last element
         */
        typename std::deque<T>::const_iterator cend() const;

        /*!
         * \return an iterator on the first element
         */
        typename std::deque<T>::iterator begin();

        /*!
         * \return an iterator on the last element
         */
        typename std::deque<T>::iterator end();

        // accessors
        /*!
         * \brief getter on the element on the i-th row, the j-th column and the z-th plan
         * \param i row index
         * \param j column index
         * \param z depth index
         * \return the value of M(i,j,k)
         */
        T operator()(std::size_t i, std::size_t j, std::size_t k) const;

        /*!
         * \brief setter on the element on the i-th row, the j-th column and the z-th plan
         * \param i row index
         * \param j column index
         * \param z depth index
         * \return a reference on M(i,j,k)
         */
        T& operator()(std::size_t i, std::size_t j, std::size_t k);

        // accessors submatrix
        Matrix<T> getCol(std::size_t) const;
        Matrix3D<T> getCols(std::size_t, std::size_t) const;
        Matrix<T> getRow(std::size_t) const;
        Matrix3D<T> getRows(std::size_t, std::size_t) const;
        Matrix<T> getSubmat(std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t) const;
        void setCol(std::size_t, const Matrix<T>&);
        void setCols(std::size_t, const Matrix3D<T>&);
        void setRow(std::size_t, const Matrix<T>&);
        void setRows(std::size_t, const Matrix3D<T>&);
        void setSubmat(std::size_t, std::size_t, std::size_t, const Matrix3D<T>&);

        // size
        std::size_t colNb() const;
        std::size_t rowNb() const;
        std::size_t depth() const;
        std::size_t size() const;

        // matrix transformation
        void fliplr();
        void flipud();
        void flipfb();
        void reshape(std::size_t, std::size_t, std::size_t);
        void rot90();
        void rot180();
        void rot270();
        void swapcol(std::size_t, std::size_t);
        void swaprow(std::size_t, std::size_t);

        // add/remove a row/column
        void delCol(std::size_t);
        void delRow(std::size_t);
        void newCol();
        void newCol(std::size_t);
        void newRow();
        void newRow(std::size_t);

        // operators
        void operator+=(const Matrix3D<T>&);
        void operator-=(const Matrix3D<T>&);
        void operator*=(const Matrix3D<T>&);

        void operator+=(const T&);
        void operator-=(const T&);
        void operator*=(const T&);
        void operator/=(const T&);
        void operator%=(const T&);

        Matrix3D<T> operator+(const Matrix3D<T>&) const;
        Matrix3D<T> operator-(const Matrix3D<T>&) const;
        Matrix3D<T> operator*(const Matrix3D<T>&) const;

        Matrix3D<T> operator+(const T&) const;
        Matrix3D<T> operator-(const T&) const;
        Matrix3D<T> operator*(const T&) const;
        Matrix3D<T> operator/(const T&) const;
        Matrix3D<T> operator%(const T&) const;

        Matrix3D<T> operator-() const;

        Matrix3D<bool> operator==(const Matrix3D<T>&) const;
        Matrix3D<bool> operator!=(const Matrix3D<T>&) const;
        Matrix3D<bool> operator>=(const Matrix3D<T>&) const;
        Matrix3D<bool> operator<=(const Matrix3D<T>&) const;
        Matrix3D<bool> operator>(const Matrix3D<T>&) const;
        Matrix3D<bool> operator<(const Matrix3D<T>&) const;

        Matrix3D<bool> operator==(const T&) const;
        Matrix3D<bool> operator!=(const T&) const;
        Matrix3D<bool> operator>=(const T&) const;
        Matrix3D<bool> operator<=(const T&) const;
        Matrix3D<bool> operator>(const T&) const;
        Matrix3D<bool> operator<(const T&) const;

};


// specific constructors
template <class T> Matrix3D<T> arange(T, T, T = T(1));
template <class T> Matrix3D<T> full(std::size_t, T);
template <class T> Matrix3D<T> full(std::size_t, std::size_t, std::size_t, T);
template <class T> Matrix3D<T> id(std::size_t);
template <class T> Matrix3D<T> id(std::size_t, std::size_t, std::size_t);
template <class T> std::pair<Matrix3D<T>, Matrix3D<T>, Matrix3D<T> > meshgrid(std::size_t);
template <class T> std::pair<Matrix3D<T>, Matrix3D<T>, Matrix3D<T> > meshgrid(std::size_t, std::size_t, std::size_t);
template <class T> Matrix3D<T> rand(std::size_t);
template <class T> Matrix3D<T> rand(std::size_t, std::size_t, std::size_t);
template <class T> Matrix3D<T> randn(std::size_t, T = T(0), T = T(1));
template <class T> Matrix3D<T> randn(std::size_t, std::size_t, std::size_t, T = T(0), T = T(1));
template <class T> Matrix3D<T> ones(std::size_t);
template <class T> Matrix3D<T> ones(std::size_t, std::size_t, std::size_t);
template <class T> Matrix3D<T> zeros(std::size_t);
template <class T> Matrix3D<T> zeros(std::size_t, std::size_t);

// apply any function
template <typename Type, class T> Matrix<T> apply(const Matrix<Type>&, T (*)(Type));
template <typename Type, class T> Matrix<T> apply(const Matrix<Type>&, T (*)(const Type&));
template <typename Type, class T> Matrix<T> apply(const Matrix<Type>&, T (Type::*)() const);
//template <typename Type, class T> Matrix<T> apply(const Matrix<Type>&, std::function<T(Type)>);
//template <typename Type, class T> Matrix<T> apply(const Matrix<Type>&, std::function<T(const Type&)>);

template <typename T1, typename T2, class T> Matrix<T> apply(const Matrix<T1>&, const Matrix<T2>&, T (*)(T1, T2));
template <typename T1, typename T2, class T> Matrix<T> apply(const Matrix<T1>&, const Matrix<T2>&, T (*)(const T1&, const T2&));

// construct a matrix depending on a condition
template <class T> Matrix<T> where(const Matrix<bool>&, const T&, const T&);
template <class T> Matrix<T> where(const Matrix<bool>&, const Matrix<T>&, const Matrix<T>&);

//concatenation
template <class T> Matrix<T> horzcat(const Matrix<T>&, const Matrix<T>&);
template <class T> Matrix<T> vertcat(const Matrix<T>&, const Matrix<T>&);

// compression
template <class T> Matrix<T> compress(const Matrix<bool>&, const Matrix<T>&, int);

// count elements
template <class T> std::size_t count(const Matrix3D<T>&, const T&);
template <class T> std::size_t count_nonzero(const Matrix3D<T>&);

// replacement
template <class T> void replace(Matrix3D<T>&, const T&, const T&);
template <class T> void replace_if(Matrix3D<T>&, const Matrix3D<bool>&, const T&);

// argwhere
template <class T = std::size_t> Matrix3D<std::size_t> argwhere(const Matrix3D<bool>&);

// unique
template <class T> Matrix3D<T> unique(const Matrix3D<T>&);

// operators
template <class T> Matrix3D<T> operator+(const T&, const Matrix3D<T>&);
template <class T> Matrix3D<T> operator-(const T&, const Matrix3D<T>&);
template <class T> Matrix3D<T> operator*(const T&, const Matrix3D<T>&);
template <class T> Matrix3D<T> operator/(const T&, const Matrix3D<T>&);

//print
template <class T> std::ostream& operator<<(std::ostream&, const Matrix3D<T>&);


#include "matrix3d.tcc"


#endif // MATRIX_HPP
