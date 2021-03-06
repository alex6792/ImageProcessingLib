/*!
 * \author Alexandre Krebs
 * \file gaussian_naive_bayes.hpp
 * \brief Gaussian Naive Bayes classifier
 */


#pragma once

#ifndef GAUSSIAN_NAIVE_BAYES_HPP
#define GAUSSIAN_NAIVE_BAYES_HPP


#include "../../matrix.hpp"


/*!
 * \class GaussianNaiveBayes
 * \brief class representing a Gaussian Naive Bayes classifier
 */
class GaussianNaiveBayes
{
    private :
        std::size_t nb_class;
        Matrix<float> means;
        Matrix<float> stdevs;
        Matrix<float> priors;

    public :
        GaussianNaiveBayes();
        Matrix<float> getCenters();
        Matrix<float> getStddevs();
        Matrix<float> getPriors();
        void fit(const Matrix<float>&, const Matrix<std::size_t>&);
        Matrix<std::size_t> fit_predict(const Matrix<float>&, const Matrix<std::size_t>&);
        Matrix<std::size_t> predict(const Matrix<float>&);
        Matrix<float> predict_proba(const Matrix<float>&);
};


#endif // GAUSSIAN_NAIVE_BAYES_HPP
