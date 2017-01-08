/*!
 * \author Alexandre Krebs
 * \file qda.hpp
 * \brief QDA classifier
 */


#pragma once

#ifndef QDA_HPP
#define QDA_HPP


#include "../../matrix.hpp"


/*!
 * \class QDA
 * \brief class representing a classifier based on Quadratic Discriminant Analysis
 */
class QDA
{
    private :
        std::size_t nb_classes;
        std::size_t nb_features;
        Matrix<float> means;
        std::vector<Matrix<float> > covars;
        Matrix<float> priors;
        Matrix<std::size_t> unique_labels;

    public :
        QDA();
        Matrix<float> getCenters();
        std::vector<Matrix<float> > getCovariances();
        Matrix<float> getPriors();
        void fit(const Matrix<float>&, const Matrix<std::size_t>&);
        Matrix<std::size_t> fit_predict(const Matrix<float>&, const Matrix<std::size_t>&);
        Matrix<std::size_t> predict(const Matrix<float>&);
        Matrix<float> predict_proba(const Matrix<float>&);
};


#endif // QDA_HPP
