/*!
 * \author Alexandre Krebs
 * \file knn.hpp
 * \brief k-Nearest Neighbor classifier
 */


#pragma once

#ifndef KNN_HPP
#define KNN_HPP


#include "../../matrix.hpp"


/*!
 * \class KNN
 * \brief class representing a K-Nearest Neighbor classifier
 */
class KNN
{
    private :
        int k;
        int nb_class;
        Matrix<float> training_data;
        Matrix<int> labels;

    public :
        explicit KNN(int);
        void fit(const Matrix<float>&, const Matrix<int>&);
        Matrix<int> fit_predict(const Matrix<float>&, const Matrix<int>&);
        Matrix<int> predict(const Matrix<float>&);
        Matrix<float> predict_proba(const Matrix<float>&);
};


#endif // KNN_HPP
