/*!
 * \author Alexandre Krebs
 * \file meanshift.hpp
 * \brief Mean shift classifier
 */


#pragma once

#ifndef MEANSHIFT_HPP
#define MEANSHIFT_HPP


#include "../../matrix.hpp"


/*!
 * \class MeanShift
 * \brief class representing a mean shift classifier
 */
class MeanShift
{
    private :
        float (*kernel)(float);
        Matrix<float> centers;
        Matrix<float> data;
        Matrix<int> labels;
        int nb_clusters;
        int nb_features;
        int max_iteration;

    public :
        MeanShift();
        void fit(const Matrix<float>&);
        Matrix<int> fit_predict(const Matrix<float>&);
        Matrix<int> predict(const Matrix<float>&);

    private :
        Matrix<float> iterate(const Matrix<float>&);
};


#endif // MEANSHIFT_HPP

