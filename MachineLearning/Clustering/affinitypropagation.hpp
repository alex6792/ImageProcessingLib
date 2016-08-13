/*!
 * \author Alexandre Krebs
 * \file affinitypropagation.hpp
 * \brief Mean shift classifier
 */


#pragma once

#ifndef AFFINITYPROPAGATION_HPP
#define AFFINITYPROPAGATION_HPP


#include "../../matrix.hpp"


/*!
 * \class AffinityPropagation
 * \brief class representing a classifier based on the Affinity Propagation algorithm
 */
class AffinityPropagation
{
    private :
        Matrix<float> centers;
        Matrix<float> availability;
        Matrix<float> responsibility;
        Matrix<std::size_t> labels;
        std::size_t nb_clusters;
        std::size_t nb_features;
        std::size_t max_iteration;

    public :
        AffinityPropagation();
        Matrix<float> getAvailability();
        Matrix<float> getResponsibility();
        Matrix<float> getCenters();
        std::size_t getMaxIteration();
        void setMaxIteration(std::size_t);
        void fit(const Matrix<float>&);
        Matrix<std::size_t> fit_predict(const Matrix<float>&);
        Matrix<std::size_t> predict(const Matrix<float>&);
};


#endif // AFFINITYPROPAGATION_HPP

