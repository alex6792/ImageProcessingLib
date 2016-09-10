/*!
 * \author Alexandre Krebs
 * \file kmeans.hpp
 * \brief K-means classifier
 */


#pragma once

#ifndef KMEANS_HPP
#define KMEANS_HPP


#include "../../matrix.hpp"


/*!
 * \class Kmeans
 * \brief class representing a k-means classifier
 */
class Kmeans
{
    private :
        Matrix<float> centers;// cluster centers
        Matrix<std::size_t> labels;// labels
        std::size_t nb_clusters;// number of clusters
        std::size_t nb_features;// number of features
        std::size_t max_iteration;// max iteration

    public :
        explicit Kmeans(std::size_t);
        Matrix<float> getCenters();
        std::size_t getMaxIteration();
        void setMaxIteration(std::size_t);
        void fit(const Matrix<float>&);
        Matrix<std::size_t> fit_predict(const Matrix<float>&);
        Matrix<std::size_t> predict(const Matrix<float>&);

        static Matrix<float> KmeansPlusPlusCenters(const Matrix<float>&, std::size_t);
        static Matrix<size_t> KmeansPlusPlusCentersIndices(const Matrix<float>&, std::size_t);

    private :
        void update_centers(const Matrix<float>&);
        void init_centers(const Matrix<float>&);
};


#endif // KMEANS_HPP
