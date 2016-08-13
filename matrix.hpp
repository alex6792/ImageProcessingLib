/*!
 * \author Alexandre Krebs
 * \file matrix.hpp
 * \brief template matrix class
 */


#pragma once

#ifndef MATRIX_HPP
#define MATRIX_HPP


#include <deque>
#include <iostream>


/*!
 * \class Matrix
 * \brief template class representing a matrix
 */
template <class T> class Matrix
{

    private :
        std::size_t sizex; /*!< number of rows*/
        std::size_t sizey; /*!< number of columns*/
        std::deque<T> mat; /*!< data*/

    public :
        // constructors
        /*!
         * \brief Default constructor that construct a matrix 1X1
         */
        Matrix();

        /*!
         * \brief Square NxN matrix constructor
         * \param N the size of the matrix
         */
        explicit Matrix(std::size_t N);

        /*!
         * \brief MxN matrix constructor
         * \param M number of rows
         * \param N number of columns
         */
        Matrix(std::size_t M, std::size_t N);

        /*!
         * \brief cast a matrix of type T1 to a matrix of type T2
         * \param M the matrix of type T1
         */
        template <typename Type> explicit Matrix(const Matrix<Type>& M);

        /*!
         * \brief construct a matrix from a C-style static array
         * \param A a static array
         */
        Matrix(std::initializer_list<std::initializer_list<T> > A);

        // iterators
        /*!
         * \return a const iterator on the first element (top left element)
         */
        typename std::deque<T>::const_iterator cbegin() const;

        /*!
         * \return a const iterator on the last element (bottom right element)
         */
        typename std::deque<T>::const_iterator cend() const;

        /*!
         * \return an iterator on the first element (top left element)
         */
        typename std::deque<T>::iterator begin();

        /*!
         * \return an iterator on the last element (bottom right element)
         */
        typename std::deque<T>::iterator end();

        // accessors
        /*!
         * \brief getter on the element on the i-th row and the j-th column
         * \param i row index
         * \param j column index
         * \return the value of M(i,j)
         */
        T operator()(std::size_t i, std::size_t j) const;

        /*!
         * \brief setter for the element on the i-th row and the j-th column
         * \param i row index
         * \param j column index
         * \return a reference on M(i,j)
         */
        T& operator()(std::size_t i, std::size_t j);

        // accessors submatrix
        Matrix<T> getCol(std::size_t) const;
        Matrix<T> getCols(std::size_t, std::size_t) const;
        Matrix<T> getRow(std::size_t) const;
        Matrix<T> getRows(std::size_t, std::size_t) const;
        Matrix<T> getSubmat(std::size_t, std::size_t, std::size_t, std::size_t) const;
        void setCol(std::size_t, const Matrix<T>&);
        void setCols(std::size_t, const Matrix<T>&);
        void setRow(std::size_t, const Matrix<T>&);
        void setRows(std::size_t, const Matrix<T>&);
        void setSubmat(std::size_t, std::size_t, const Matrix<T>&);

        // size
        std::size_t colNb() const;
        std::size_t rowNb() const;
        std::size_t size() const;

        // matrix transformation
        void fliplr();
        void flipud();
        void reshape(std::size_t, std::size_t);
        void rot90();
        void rot180();
        void rot270();
        void swapcol(std::size_t, std::size_t);
        void swaprow(std::size_t, std::size_t);
        void transpose();

        // add/remove a row/column
        void delCol(std::size_t);
        void delRow(std::size_t);
        void newCol();
        void newCol(std::size_t);
        void newRow();
        void newRow(std::size_t);

        // operators
        void operator+=(const Matrix<T>&);
        void operator-=(const Matrix<T>&);
        void operator*=(const Matrix<T>&);

        void operator+=(const T&);
        void operator-=(const T&);
        void operator*=(const T&);
        void operator/=(const T&);
        void operator%=(const T&);

        Matrix<T> operator+(const Matrix<T>&) const;
        Matrix<T> operator-(const Matrix<T>&) const;
        Matrix<T> operator*(const Matrix<T>&) const;

        Matrix<T> operator+(const T&) const;
        Matrix<T> operator-(const T&) const;
        Matrix<T> operator*(const T&) const;
        Matrix<T> operator/(const T&) const;
        Matrix<T> operator%(const T&) const;

        Matrix<T> operator-() const;

        Matrix<bool> operator==(const Matrix<T>&) const;
        Matrix<bool> operator!=(const Matrix<T>&) const;
        Matrix<bool> operator>=(const Matrix<T>&) const;
        Matrix<bool> operator<=(const Matrix<T>&) const;
        Matrix<bool> operator>(const Matrix<T>&) const;
        Matrix<bool> operator<(const Matrix<T>&) const;

        Matrix<bool> operator==(const T&) const;
        Matrix<bool> operator!=(const T&) const;
        Matrix<bool> operator>=(const T&) const;
        Matrix<bool> operator<=(const T&) const;
        Matrix<bool> operator>(const T&) const;
        Matrix<bool> operator<(const T&) const;

};

// specific constructors
template <class T> Matrix<T> arange(T, T, T = T(1));
template <class T> Matrix<T> full(std::size_t, T);
template <class T> Matrix<T> full(std::size_t, std::size_t, T);
template <class T> Matrix<T> id(std::size_t);
template <class T> Matrix<T> id(std::size_t, std::size_t);
template <class T> std::pair<Matrix<T>, Matrix<T> > meshgrid(std::size_t);
template <class T> std::pair<Matrix<T>, Matrix<T> > meshgrid(std::size_t, std::size_t);
template <class T> Matrix<T> rand(std::size_t);
template <class T> Matrix<T> rand(std::size_t, std::size_t);
template <class T> Matrix<T> randn(std::size_t, T = T(0), T = T(1));
template <class T> Matrix<T> randn(std::size_t, std::size_t, T = T(0), T = T(1));
template <class T> Matrix<T> ones(std::size_t);
template <class T> Matrix<T> ones(std::size_t, std::size_t);
template <class T> Matrix<T> zeros(std::size_t);
template <class T> Matrix<T> zeros(std::size_t, std::size_t);

// apply any function
template <class T, typename Type> Matrix<Type> apply(const Matrix<T>&, Type (*)(T));
template <class T, typename Type> Matrix<Type> apply(const Matrix<T>&, Type (T::*)() const);

// construct a matrix depending on a condition
template <class T> Matrix<T> where(const Matrix<bool>&, const T&, const T&);
template <class T> Matrix<T> where(const Matrix<bool>&, const Matrix<T>&, const Matrix<T>&);

// compression
template <class T> Matrix<T> compress(const Matrix<bool>&, const Matrix<T>&, int);

//count elements
template <class T> std::size_t count(const Matrix<T>&, const T&);
template <class T> std::size_t count_nonzero(const Matrix<T>&);

// replacement
template <class T> void replace(Matrix<T>&, const T&, const T&);
template <class T> void replace_if(Matrix<T>&, const Matrix<bool>&, const T&);

//argwhere
template <class T = std::size_t> Matrix<std::size_t> argwhere(const Matrix<bool>&);

// unique
template <class T> Matrix<T> unique(const Matrix<T>&);

// get diagonal
template <class T> Matrix<T> diag(const Matrix<T>&);

//non member function transpose
template <class T> Matrix<T> transpose(const Matrix<T>&);

// operators
template <class T> Matrix<T> operator+(const T&, const Matrix<T>&);
template <class T> Matrix<T> operator-(const T&, const Matrix<T>&);
template <class T> Matrix<T> operator*(const T&, const Matrix<T>&);
template <class T> Matrix<T> operator/(const T&, const Matrix<T>&);

//prstd::size_t
template <class T> std::ostream& operator<<(std::ostream&, const Matrix<T>&);


#include "matrix.tcc"


#endif // MATRIX_HPP
