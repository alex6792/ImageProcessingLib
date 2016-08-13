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
        Matrix<std::size_t> centers;// cluster centers
        Matrix<std::size_t> labels;// labels
        std::size_t nb_clusters;// number of clusters
        std::size_t nb_features;
        Matrix<float> distances;
        std::vector<Matrix<float> > data;
        float old_cost, cur_cost;

    public :
        explicit Kmedoids(std::size_t);
        Matrix<float> getCenters();
        void fit(const Matrix<float>&);
        Matrix<std::size_t> fit_predict(const Matrix<float>&);
        Matrix<std::size_t> predict(const Matrix<float>&);
};


#endif // KMEDOIDS_HPP

