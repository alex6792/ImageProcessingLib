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
        std::size_t max_iteration;// maximum iteration allowed

        std::vector<Matrix<float> > A;
        std::vector<Matrix<float> > y, z;
        std::vector<Matrix<float> > biases;
        std::vector<std::size_t> hidden_layers_sizes;
        std::size_t nb_clusters;
        std::size_t nb_hidden_layers;
        std::size_t nb_features;
        float (*activation_fct)(float);
        float (*activation_fct_der)(float);

    public :
        explicit MLP(const std::vector<std::size_t>&);
        void export_graphviz(std::string);
        void set_alpha(float);
        void set_max_iteration(std::size_t);
        void fit(const Matrix<float>&, const Matrix<std::size_t>&);
        Matrix<std::size_t> fit_predict(const Matrix<float>&, const Matrix<std::size_t>&);
        Matrix<std::size_t> predict(const Matrix<float>&);
};


#endif // MLP_HPP
