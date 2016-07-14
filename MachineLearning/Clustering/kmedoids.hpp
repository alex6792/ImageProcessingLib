/*!
 * \author Alexandre Krebs
 * \file kmedoids.hpp
 * \brief K-medoids classifier
 */


#pragma once

#ifndef KMEDOIDS_HPP
#define KMEDOIDS_HPP


#include "../../matrix.hpp"


/*!
 * \class Kmedoids
 * \brief class representing a k-medoids classifier
 */
class Kmedoids
{
    private :
        Matrix<int> centers;// cluster centers
        Matrix<int> labels;// labels
        int nb_clusters;// number of clusters
        int nb_features;
        Matrix<float> distances;
        std::vector<Matrix<float> > data;
        float old_cost, cur_cost;

    public :
        explicit Kmedoids(int);
        Matrix<float> getCenters();
        void fit(const Matrix<float>&);
        Matrix<int> fit_predict(const Matrix<float>&);
        Matrix<int> predict(const Matrix<float>&);
};


#endif // KMEDOIDS_HPP

