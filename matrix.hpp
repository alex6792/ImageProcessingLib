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
        int sizex; /*!< number of rows*/
        int sizey; /*!< number of columns*/
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
        explicit Matrix(int N);

        /*!
         * \brief MxN matrix constructor
         * \param M number of rows
         * \param N number of columns
         */
        Matrix(int M, int N);

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
        T operator()(int i, int j) const;

        /*!
         * \brief setter for the element on the i-th row and the j-th column
         * \param i row index
         * \param j column index
         * \return a reference on M(i,j)
         */
        T& operator()(int i, int j);

        // accessors submatrix
        Matrix<T> getCol(int) const;
        Matrix<T> getCols(int, int) const;
        Matrix<T> getRow(int) const;
        Matrix<T> getRows(int, int) const;
        Matrix<T> getSubmat(int, int, int, int) const;
        void setCol(int, const Matrix<T>&);
        void setCols(int, const Matrix<T>&);
        void setRow(int, const Matrix<T>&);
        void setRows(int, const Matrix<T>&);
        void setSubmat(int, int, const Matrix<T>&);

        // size
        int colNb() const;
        int rowNb() const;
        int size() const;

        // matrix transformation
        void fliplr();
        void flipud();
        void reshape(int, int);
        void rot90();
        void rot180();
        void rot270();
        void swapcol(int, int);
        void swaprow(int, int);
        void transpose();

        // add/remove a row/column
        void delCol(int);
        void delRow(int);
        void newCol();
        void newCol(int);
        void newRow();
        void newRow(int);

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
template <class T> Matrix<T> full(int, T);
template <class T> Matrix<T> full(int, int, T);
template <class T> Matrix<T> id(int);
template <class T> Matrix<T> id(int, int);
template <class T> std::pair<Matrix<T>, Matrix<T> > meshgrid(int);
template <class T> std::pair<Matrix<T>, Matrix<T> > meshgrid(int, int);
template <class T> Matrix<T> rand(int);
template <class T> Matrix<T> rand(int, int);
template <class T> Matrix<T> randn(int, T = T(0), T = T(1));
template <class T> Matrix<T> randn(int, int, T = T(0), T = T(1));
template <class T> Matrix<T> ones(int);
template <class T> Matrix<T> ones(int, int);
template <class T> Matrix<T> zeros(int);
template <class T> Matrix<T> zeros(int, int);

// apply any function
template <class T, typename Type> Matrix<Type> apply(const Matrix<T>&, Type (*)(T));
template <class T, typename Type> Matrix<Type> apply(const Matrix<T>&, Type (T::*)() const);

// construct a matrix depending on a condition
template <class T> Matrix<T> where(const Matrix<bool>&, const T&, const T&);
template <class T> Matrix<T> where(const Matrix<bool>&, const Matrix<T>&, const Matrix<T>&);

// compression
template <class T> Matrix<T> compress(const Matrix<bool>&, const Matrix<T>&, int);

//count elements
template <class T> int count(const Matrix<T>&, const T&);
template <class T> int count_nonzero(const Matrix<T>&);

// replacement
template <class T> void replace(Matrix<T>&, const T&, const T&);
template <class T> void replace_if(Matrix<T>&, const Matrix<bool>&, const T&);

//argwhere
template <class T = int> Matrix<int> argwhere(const Matrix<bool>&);

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

//print
template <class T> std::ostream& operator<<(std::ostream&, const Matrix<T>&);


#include "matrix.tcc"


#endif // MATRIX_HPP
