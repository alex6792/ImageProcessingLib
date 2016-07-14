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
        Matrix<int> labels;// labels
        int nb_clusters;// number of clusters
        int max_iteration;// max iteration

    public :
        explicit Kmeans(int);
        Matrix<float> getCenters();
        int getMaxIteration();
        void setMaxIteration(int);
        void fit(const Matrix<float>&);
        Matrix<int> fit_predict(const Matrix<float>&);
        Matrix<int> predict(const Matrix<float>&);

        static Matrix<float> KmeansPlusPlusInit(const Matrix<float>&, int);

    private :
        void update_centers(const Matrix<float>&);
        void init_centers(const Matrix<float>&);
};


#endif // KMEANS_HPP
