/*!
 * \author Alexandre Krebs
 * \file lda.hpp
 * \brief LDA classifier
 */


#pragma once

#ifndef LDA_HPP
#define LDA_HPP


#include "../../matrix.hpp"


/*!
 * \class LDA
 * \brief class representing a classifier based on Linear Discriminant Analysis
 */
class LDA
{
    private :
        std::size_t nb_classes;
        std::size_t nb_features;
        Matrix<float> means;
        Matrix<float> covars;
        Matrix<float> priors;
        Matrix<std::size_t> unique_labels;

    public :
        LDA();
        Matrix<float> getCenters();
        Matrix<float> getCovariances();
        Matrix<float> getPriors();
        void fit(const Matrix<float>&, const Matrix<std::size_t>&);
        Matrix<std::size_t> fit_predict(const Matrix<float>&, const Matrix<std::size_t>&);
        Matrix<std::size_t> predict(const Matrix<float>&);
        Matrix<float> predict_proba(const Matrix<float>&);
};


#endif // LDA_HPP
