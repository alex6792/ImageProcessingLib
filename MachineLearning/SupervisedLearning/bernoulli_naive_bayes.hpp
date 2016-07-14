/*!
 * \author Alexandre Krebs
 * \file bernoulli_naive_bayes.hpp
 * \brief Bernoulli Naive Bayes classifier
 */


#pragma once

#ifndef BERNOULLI_NAIVE_BAYES_HPP
#define BERNOULLI_NAIVE_BAYES_HPP


#include "../../matrix.hpp"


/*!
 * \class BernoulliNaiveBayes
 * \brief class representing a Bernoulli Naive Bayes classifier
 */
class BernoulliNaiveBayes
{
    private :
        int nb_class;
        Matrix<float> proba;
        Matrix<float> priors;

    public :
        BernoulliNaiveBayes();
        void fit(const Matrix<bool>&, const Matrix<int>&);
        Matrix<int> fit_predict(const Matrix<bool>&, const Matrix<int>&);
        Matrix<int> predict(const Matrix<bool>&);
        Matrix<float> predict_proba(const Matrix<bool>&);
};


#endif // BERNOULLI_NAIVE_BAYES_HPP
