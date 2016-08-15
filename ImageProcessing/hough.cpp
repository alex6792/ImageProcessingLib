#include "hough.hpp"
#include "../linalg.hpp"
#include "../mmath.hpp"

Matrix<std::size_t> hough_transform(const Matrix<bool>& M)
{
    float dtheta = 2*PI/32.0f;
    float dphi = 1.0f;
    float phi_max = 2*(M.rowNb()+M.colNb())-1;//std::ceil(std::sqrt((M.rowNb()+1)*(M.rowNb()+1)+(M.colNb()+1)*(M.colNb()+1)));
    std::size_t W = std::ceil(2*PI/dtheta), H = std::ceil(phi_max/dphi);
    Matrix<std::size_t> h_transform = zeros<std::size_t>(H, W);

    Matrix<float> argw = Matrix<float>(argwhere(M));
    Matrix<float> theta = transpose(arange<float>(0, 2*PI+dtheta, dtheta));
    Matrix<float> cos_angles = cos(theta);
    Matrix<float> sin_angles = sin(theta);
    //std::cout<<cos_angles.rowNb()<<" "<<cos_angles.colNb()<<std::endl;
    //std::cout<<argw.rowNb()<<" "<<argw.colNb()<<std::endl;
    //std::cout<<h_transform.rowNb()<<" "<<h_transform.colNb()<<std::endl;
    Matrix<float> X = argw.getCol(1);
    Matrix<float> Y = M.rowNb()-1.0f-argw.getCol(0);
    Matrix<float> phi = dot(X, cos_angles)+dot(Y, sin_angles);
    theta = dot(ones<float>(argw.rowNb(), 1), theta);
    //std::cout<<phi.rowNb()<<" "<<phi.colNb()<<std::endl;
    //std::cout<<theta.rowNb()<<" "<<theta.colNb()<<std::endl;

    theta/=dtheta;
    phi+=phi_max/2.0f;
    //phi/=2.0f;
    phi = apply<float, float>(phi, std::ceil);
    //std::cout<<phi<<std::endl;
    //std::cout<<phi<<std::endl;
    //auto it_theta = theta.cbegin(), it_phi = phi.cbegin();
    for(auto it_theta = theta.cbegin(), it_end = theta.cend(), it_phi = phi.cbegin();it_theta<it_end;++it_theta, ++it_phi)
    {
        ++h_transform(*it_phi, *it_theta);
    }

    return h_transform;
}

