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
        Matrix<std::size_t> labels;
        std::size_t nb_clusters;
        std::size_t nb_features;

    public :
        explicit HierarchicalClustering(std::size_t);
        void fit(const Matrix<float>&);
        Matrix<std::size_t> fit_predict(const Matrix<float>&);
        Matrix<std::size_t> predict(const Matrix<float>&);
};


#endif // HIERARCHICALCLUSTERING_HPP


