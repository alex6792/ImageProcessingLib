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
        Matrix<int> labels;// labels
        int nb_clusters;// number of cluster
        int nb_features;// number of features
        int max_iteration;// max iteration

    public :
        explicit GMM(int);
        Matrix<float> getCenters();
        Matrix<float> getStddevs();
        Matrix<float> getWeights();
        void fit(const Matrix<float>&);
        Matrix<int> fit_predict(const Matrix<float>&);
        Matrix<int> predict(const Matrix<float>&);

    private :
        void expectation(const Matrix<float>&);
        void maximization(const Matrix<float>&);
};


#endif // GMM_HPP

