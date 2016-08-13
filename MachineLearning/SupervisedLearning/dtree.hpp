/*!
 * \author Alexandre Krebs
 * \file dtree.hpp
 * \brief Decision Tree classifier
 */


#pragma once

#ifndef DTREE_HPP
#define DTREE_HPP


#include "../../matrix.hpp"



struct DtreeNode
{
    public :
        std::size_t feature;
        float thresh;
        float gini;
        std::size_t samples;
        std::size_t label;
        Matrix<float> value;
        std::size_t left_child;
        std::size_t right_child;
        bool is_a_leaf;

        DtreeNode();
        void show();
};


/*!
 * \class Dtree
 * \brief class representing a Decision Tree classifier
 */
class Dtree
{
    private :
        std::vector<DtreeNode> tree;
        std::size_t nb_clusters;
        std::size_t nb_features;

    public :
        Dtree();
        void export_graphviz(std::string);
        void fit(const Matrix<float>&, const Matrix<std::size_t>&);
        Matrix<std::size_t> fit_predict(const Matrix<float>&, const Matrix<std::size_t>&);
        Matrix<std::size_t> predict(const Matrix<float>&);
    //private :
        std::pair<std::size_t, float> find_best_split(const Matrix<float>&, const Matrix<std::size_t>&);
        void split_node(DtreeNode&, const Matrix<float>&, const Matrix<std::size_t>&);
};


#endif // DTREE_HPP

