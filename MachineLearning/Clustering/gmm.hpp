/*!
 * \author Alexandre Krebs
 * \file gmm.hpp
 * \brief Gaussian Mixture Model classifier
 */


#pragma once

#ifndef GMM_HPP
#define GMM_HPP


#include "../../matrix.hpp"


/*!
 * \class GMM
 * \brief class representing a gaussian mixture based classifier
 */
class GMM
{
    private :
        Matrix<float> centers;// cluster centers
        Matrix<float> vars;// cluster standard deviation
        Matrix<float> weights;// cluster weights
        Matrix<float> responsibilities;
        Matrix<std::size_t> labels;// labels
        std::size_t nb_clusters;// number of cluster
        std::size_t nb_features;// number of features
        std::size_t max_iteration;// max iteration

    public :
        explicit GMM(std::size_t);
        Matrix<float> getCenters();
        Matrix<float> getStddevs();
        Matrix<float> getWeights();
        void fit(const Matrix<float>&);
        Matrix<std::size_t> fit_predict(const Matrix<float>&);
        Matrix<std::size_t> predict(const Matrix<float>&);

    private :
        void expectation(const Matrix<float>&);
        void maximization(const Matrix<float>&);
};


#endif // GMM_HPP

