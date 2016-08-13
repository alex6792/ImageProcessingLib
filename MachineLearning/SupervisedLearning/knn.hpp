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
        std::size_t k;
        std::size_t nb_class;
        Matrix<float> training_data;
        Matrix<std::size_t> labels;

    public :
        explicit KNN(std::size_t);
        void fit(const Matrix<float>&, const Matrix<std::size_t>&);
        Matrix<std::size_t> fit_predict(const Matrix<float>&, const Matrix<std::size_t>&);
        Matrix<std::size_t> predict(const Matrix<float>&);
        Matrix<float> predict_proba(const Matrix<float>&);
};


#endif // KNN_HPP
