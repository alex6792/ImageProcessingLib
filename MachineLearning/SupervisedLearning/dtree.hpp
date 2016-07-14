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
        int feature;
        float thresh;
        float gini;
        int samples;
        int label;
        Matrix<float> value;
        int left_child;
        int right_child;

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
        int nb_clusters;
        int nb_features;

    public :
        Dtree();
        void export_graphviz(std::string);
        void fit(const Matrix<float>&, const Matrix<int>&);
        Matrix<int> fit_predict(const Matrix<float>&, const Matrix<int>&);
        Matrix<int> predict(const Matrix<float>&);
    //private :
        std::pair<int, float> find_best_split(const Matrix<float>&, const Matrix<int>&);
        void split_node(DtreeNode&, const Matrix<float>&, const Matrix<int>&);
};


#endif // DTREE_HPP

