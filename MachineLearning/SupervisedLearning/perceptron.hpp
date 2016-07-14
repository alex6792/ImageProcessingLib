/*!
 * \author Alexandre Krebs
 * \file perceptron.hpp
 * \brief Perceptron classifier
 */


#pragma once

#ifndef PERCEPTRON_HPP
#define PERCEPTRON_HPP


#include "../../matrix.hpp"


/*!
 * \class Perceptron
 * \brief class representing a Perceptron classifier
 */
class Perceptron
{
    private :
        Matrix<float> coeff;// weight
        float thresh;// intercept, first coeff
        float alpha;// learning rate
        int max_iteration;// maximum iteration allowed

    public :
        Perceptron();
        Matrix<float> get_coeff();
        float get_thresh();
        void set_alpha(float);
        void set_max_iteration(int);
        void fit(const Matrix<float>&, const Matrix<bool>&);
        Matrix<bool> fit_predict(const Matrix<float>&, const Matrix<bool>&);
        Matrix<bool> predict(const Matrix<float>&);
};


#endif // PERCEPTRON_HPP
