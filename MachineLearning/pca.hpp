/*!
 * \author Alexandre Krebs
 * \file pca.hpp
 * \brief Principal Component Analysis
 */


#pragma once

#ifndef PCA_HPP
#define PCA_HPP


#include "../matrix.hpp"


/*!
 * \class PCA
 * \brief class for Principal Component Analysis decomposition
 */
class PCA
{
    private :
        Matrix<float> mean; /*!< mean matrix*/
        Matrix<float> covar; /*!< covariance matrix*/
        Matrix<float> eigenvalues; /*!< eigenvalues*/
        Matrix<float> eigenvectors; /*!< eigenvectors*/

    public :
        PCA();
        void fit(const Matrix<float>&);
        Matrix<float> fit_transform(const Matrix<float>&);
        Matrix<float> inverse_transform(const Matrix<float>&);
        Matrix<float> transform(const Matrix<float>&);
};


#endif // PCA_HPP
