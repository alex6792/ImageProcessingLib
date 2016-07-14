/*!
 * \author Alexandre Krebs
 * \file hierarchicalclustering.hpp
 * \brief Hierarchical Clustering classifier
 */


#pragma once

#ifndef HIERARCHICALCLUSTERING_HPP
#define HIERARCHICALCLUSTERING_HPP


#include "../../matrix.hpp"


/*!
 * \class HierarchicalClustering
 * \brief class representing a hierarchical clustering based classifier
 */
class HierarchicalClustering
{
    private :
        Matrix<float> data;
        Matrix<int> labels;
        int nb_clusters;
        int nb_features;

    public :
        explicit HierarchicalClustering(int);
        void fit(const Matrix<float>&);
        Matrix<int> fit_predict(const Matrix<float>&);
        Matrix<int> predict(const Matrix<float>&);
};


#endif // HIERARCHICALCLUSTERING_HPP


