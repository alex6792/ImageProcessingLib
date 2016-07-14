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
        Matrix<int> labels;
        int nb_clusters;
        int nb_features;
        int max_iteration;

    public :
        AffinityPropagation();
        Matrix<float> getAvailability();
        Matrix<float> getResponsibility();
        Matrix<float> getCenters();
        int getMaxIteration();
        void setMaxIteration(int);
        void fit(const Matrix<float>&);
        Matrix<int> fit_predict(const Matrix<float>&);
        Matrix<int> predict(const Matrix<float>&);
};


#endif // AFFINITYPROPAGATION_HPP

