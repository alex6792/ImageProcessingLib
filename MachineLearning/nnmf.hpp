/*!
 * \author Alexandre Krebs
 * \file nnmf.hpp
 * \brief Non-negative matrix factorization
 */


#pragma once

#ifndef NNMF_HPP
#define NNMF_HPP


#include "../matrix.hpp"


/*!
 * \class NNMF
 * \brief class for NNMF decomposition
 */
class NNMF
{
    private :
        std::size_t nb_components; /*!< number of components*/
        Matrix<float> W; /*!< W matrix*/
        Matrix<float> H; /*!< H matrix*/

    public :
        NNMF(std::size_t = 0);
        void fit(const Matrix<float>&);
        Matrix<float> fit_transform(const Matrix<float>&);
        Matrix<float> inverse_transform(const Matrix<float>&);
        Matrix<float> transform(const Matrix<float>&);
};


#endif // NNMF_HPP
