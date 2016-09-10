#include "shape_features.hpp"
#include "morphology.hpp"
#include "../mmath.hpp"
#include "../statistics.hpp"


float circularity(const Matrix<bool>& M)
{
    float n1 = count(M, true);
    Matrix<bool> G = inner_gradient(M);
    float n2 = count(G, true);
    return 4*PI*n1/(n2*n2);
}


float major_axis(const Matrix<bool> & M)
{
    Matrix<float> M_copy(M);
    float result = 2*central_moment(M_copy,1,1)/(central_moment(M_copy,2,0)-central_moment(M_copy,0,2));
    return 0.5f*std::atan(result);
}

Zernike::Zernike(int m_arg, int n_arg)
{
    m = m_arg;
    n = n_arg;

    if(n==0 && m==0)
    {
        rho_fct = [](float rho){return 1.0f;};
    }
    else if(n==1 && (m==1||m==-1))
    {
        rho_fct = [](float rho){return rho;};
    }
    else if(n==2 && m==0)
    {
        rho_fct = [](float rho){return 2*rho*rho-1;};
    }
    else if(n==2 && (m==2||m==-2))
    {
        rho_fct = [](float rho){return rho*rho;};
    }
    else if(n==3 && (m==1||m==-1))
    {
        rho_fct = [](float rho){return rho*(3*rho*rho-2);};
    }
    else if(n==3 && (m==3||m==-3))
    {
        rho_fct = [](float rho){return rho*rho*rho;};
    }
    else if(n==4 && m==0)
    {
        rho_fct = [](float rho){return 6.0f*std::pow(rho, 4.0f)-6.0f*std::pow(rho, 2.0f)+1.0f;};
    }
    else if(n==4 && (m==2||m==-2))
    {
        rho_fct = [](float rho){return 4.0f*std::pow(rho, 4.0f)-3.0f*std::pow(rho, 2.0f);};
    }
    else if(n==4 && (m==4||m==-4))
    {
        rho_fct = [](float rho){return std::pow(rho, 4.0f);};
    }
    else if(n==5 && (m==1||m==-1))
    {
        rho_fct = [](float rho){return 10.0f*std::pow(rho, 5.0f)-12.0f*std::pow(rho, 3.0f)+3.0f*rho;};
    }
    else if(n==5 && (m==3||m==-3))
    {
        rho_fct = [](float rho){return 5.0f*std::pow(rho, 5.0f)-4.0f*std::pow(rho, 3.0f);};
    }
    else if(n==5 && (m==5||m==-5))
    {
        rho_fct = [](float rho){return std::pow(rho, 5.0f);};
    }
    else
    {
        std::cout<<"error, unknown function"<<std::endl;
    }
}

std::complex<float> Zernike::polynom(float rho, float phi)
{
    return std::exp(std::complex<float>(0, m*phi))*rho_fct(rho);
}
