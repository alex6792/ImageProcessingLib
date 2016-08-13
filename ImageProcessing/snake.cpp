#include "../linalg.hpp"
#include "../mmath.hpp"
#include "linear_filtering.hpp"
#include "snake.hpp"


Snake::Snake()
{
    lambda1 = 5;
    lambda2 = 5;
    lambda3 = 0.5;
    nb_points = 20;
    dt = 0.01;
    max_iteration = 200;
}

std::pair<Matrix<float>, Matrix<float> > Snake::getcoordinates()
{
    return std::make_pair(Vx, Vy);
}

std::size_t Snake::getMaxIteration()
{
    return max_iteration;
}

std::size_t Snake::getNbPoints()
{
    return nb_points;
}

float Snake::getDeltaT()
{
    return  dt;
}

float Snake::getLambda1()
{
    return lambda1;
}

float Snake::getLambda2()
{
    return lambda2;
}

float Snake::getLambda3()
{
    return lambda3;
}

Matrix<float> Snake::getgradient()
{
    return Px;
}

void Snake::setCoordinates(const Matrix<float>& Vx_arg, const Matrix<float>& Vy_arg)
{
    if(Vx_arg.rowNb()==Vy_arg.rowNb() && Vx_arg.colNb()==1 && Vy_arg.colNb()==1)
    {
        Vx = Vx_arg;
        Vy = Vy_arg;
        nb_points = Vx_arg.rowNb();
        D2 = -2.0f*id<float>(nb_points);
        for(std::size_t i=0;i<nb_points-1;++i)
        {
            D2(i, i+1) = 1.0f;
            D2(i+1, i) = 1.0f;
        }
        D2(0, nb_points-1) = 1.0f;
        D2(nb_points-1, 0) = 1.0f;
        D4 = dot(D2, D2);
        A = id<float>(nb_points)+dt*(-lambda1*D2+lambda2*D4);
        A_inv = inv(A);
    }
    else
        std::cout<<"invalid matrix size"<<std::endl;

}

void Snake::setDeltaT(float dt_arg)
{
    dt = dt_arg;
    A = id<float>(nb_points)+dt*(-lambda1*D2+lambda2*D4);
    A_inv = inv(A);
}

void Snake::setLambda1(float lambda1_arg)
{
    lambda1 = lambda1_arg;
    A = id<float>(nb_points)+dt*(-lambda1*D2+lambda2*D4);
    A_inv = inv(A);
}

void Snake::setLambda2(float lambda2_arg)
{
    lambda2 = lambda2_arg;
    A = id<float>(nb_points)+dt*(-lambda1*D2+lambda2*D4);
    A_inv = inv(A);
}

void Snake::setLambda3(float lambda3_arg)
{
    lambda3 = lambda3_arg;
}

void Snake::setNbPoints(float nb_points_arg)
{
    nb_points = nb_points_arg;
}

void Snake::init_image(const Matrix<unsigned char>& I_arg)
{
    I = Matrix<float>(I_arg);
    Matrix<float> Gx = conv(I, sobelx());
    Matrix<float> Gy = conv(I, sobely());
    Matrix<float> G = Gx*Gx+Gy*Gy;
    Px = conv(G, sobelx());
    Py = conv(G, sobely());
}

void Snake::init_contour(float xc, float yc, float r)
{
    setCoordinates(r*cos(arange<float>(0, nb_points)*2*PI/nb_points)+xc,
                   r*sin(arange<float>(0, nb_points)*2*PI/nb_points)+yc);
}

void Snake::init_contour(const Matrix<float>& Vx_arg, const Matrix<float>& Vy_arg)
{
    setCoordinates(Vx_arg, Vy_arg);
}

void Snake::iterate()
{
    Matrix<float> Gx(Vx.rowNb(), Vx.colNb());
    Matrix<float> Gy(Vy.rowNb(), Vy.colNb());
    for(std::size_t i=0;i<Vx.rowNb();++i)
    {
        Gx(i, 0) = Px((int)Vx(i, 0), (int)Vy(i, 0));
        Gy(i, 0) = Py((int)Vx(i, 0), (int)Vy(i, 0));
    }

    Vx = dot(A_inv, Vx+dt*lambda3*Gx);
    Vy = dot(A_inv, Vy+dt*lambda3*Gy);

    replace_if(Vx, Vx<0, 0.0f);
    replace_if(Vy, Vy<0, 0.0f);
    replace_if(Vx, Vx>=I.rowNb(), float(I.rowNb()-1));
    replace_if(Vy, Vy>=I.colNb(), float(I.colNb()-1));
}

void Snake::run()
{
    for(std::size_t i=0;i<max_iteration;++i)
    {
        iterate();
    }
}
