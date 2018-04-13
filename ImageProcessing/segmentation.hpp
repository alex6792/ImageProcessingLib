/*!
 * \author Alexandre Krebs
 * \file segmentation.hpp
 * \brief Segmentation Algorithms
 */


#pragma once

#ifndef SEGMENTATION_HPP
#define SEGMENTATION_HPP


#include "../matrix.hpp"


struct Edge
{
    public :
        unsigned char w;
        std::size_t v1, v2;

        Edge();
        Edge(std::size_t, std::size_t, unsigned char);
        bool operator<(const Edge&) const;
};


Matrix<bool> hysteresis(const Matrix<unsigned char>&, unsigned char, unsigned char);
Matrix<bool> otsu(const Matrix<unsigned char>&);
Matrix<std::size_t> felzenszwalb(const Matrix<unsigned char>&);
Matrix<std::size_t> watershed(const Matrix<unsigned char>&, const Matrix<std::size_t>&);


#endif // SEGMENTATION_HPP
