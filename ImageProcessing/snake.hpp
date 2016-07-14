/*!
 * \author Alexandre Krebs
 * \file snake.hpp
 * \brief active contours
 */


#pragma once

#ifndef SNAKE_HPP
#define SNAKE_HPP


#include "../matrix.hpp"


/*!
 * \class Snake
 * \brief class representing an active contour
 */
class Snake
{
    private :
        Matrix<float> Vx; /*!< x coordinates*/
        Matrix<float> Vy; /*!< y coordinates*/
        int max_iteration; /*!< max iterations*/
        int nb_points; /*!< number of points*/
        float dt; /*!< step*/
        float lambda1, lambda2, lambda3; /*!< energy weights*/
        Matrix<float> D2, D4; /*!< derivative matrices*/
        Matrix<float> A, A_inv; /*!< transition matrices*/
        Matrix<float> Px, Py; /*!< gradient matrices*/
        Matrix<float> I; /*!< image*/

    public :
        Snake();

        std::pair<Matrix<float>, Matrix<float> > getcoordinates();
        int getMaxIteration();
        int getNbPoints();
        float getDeltaT();
        float getLambda1();
        float getLambda2();
        float getLambda3();
        Matrix<float> getgradient();

        void setCoordinates(const Matrix<float>&, const Matrix<float>&);
        void setDeltaT(float);
        void setLambda1(float);
        void setLambda2(float);
        void setLambda3(float);
        void setNbPoints(float);

        void init_image(const Matrix<unsigned char>&);
        void init_contour(float, float, float);
        void init_contour(const Matrix<float>&, const Matrix<float>&);
        void iterate();
        void run();

};

#endif // MATRIX_HPP
