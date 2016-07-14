/*!
 * \author Alexandre Krebs
 * \file mlp.hpp
 * \brief Multi Layer Perceptron classifier
 */


#pragma once

#ifndef MLP_HPP
#define MLP_HPP


#include "../../matrix.hpp"


/*!
 * \class MLP
 * \brief class representing a Multi Layer Perceptron classifier
 */
class MLP
{
    private :

        float alpha;// learning rate
        int max_iteration;// maximum iteration allowed

        std::vector<Matrix<float> > A;
        std::vector<Matrix<float> > y, z;
        std::vector<Matrix<float> > biases;
        std::vector<int> hidden_layers_sizes;
        int nb_clusters;
        int nb_hidden_layers;
        int nb_features;
        float (*activation_fct)(float);
        float (*activation_fct_der)(float);

    public :
        explicit MLP(const std::vector<int>&);
        void export_graphviz(std::string);
        void set_alpha(float);
        void set_max_iteration(int);
        void fit(const Matrix<float>&, const Matrix<int>&);
        Matrix<int> fit_predict(const Matrix<float>&, const Matrix<int>&);
        Matrix<int> predict(const Matrix<float>&);
};


#endif // MLP_HPP
