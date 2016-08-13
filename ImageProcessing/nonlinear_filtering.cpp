#include <float.h>
#include "../statistics.hpp"
#include "linear_filtering.hpp"
#include "nonlinear_filtering.hpp"


Matrix<unsigned char> bilateral(Matrix<unsigned char> M, std::size_t filtersize, float range_var, float spatial_var)
{
    Matrix<unsigned char> filtered_img(M.rowNb(), M.colNb());
    for(std::size_t i=0, I=M.rowNb();i<I;++i)
    {
        for(std::size_t j=0, J=M.colNb();j<J;++j)
        {
            float weighted_sum = 0.0f;
            float normalization = 0.0f;
            for(std::size_t k=0;k<filtersize;++k)
            {
                for(std::size_t l=0;l<filtersize;++l)
                {

                    if(i+k>=filtersize/2 && j+l>=filtersize/2 && i+k<I+filtersize/2 && j+l<J+filtersize/2)
                    {
                        std::size_t dist = filtersize*filtersize/2+k*k+l*l-2*(k+l)*filtersize;
                        std::size_t x = i+k-filtersize/2;
                        std::size_t y = j+l-filtersize/2;
                        float w = std::exp(-((float)M(i, j)-(float)M(x, y))*((float)M(i, j)-(float)M(x, y))/(2.0f*range_var));
                        w*=std::exp(-(dist)/(2.0f*spatial_var));
                        weighted_sum+=w*M(x, y);
                        normalization+=w;
                    }
                }
            }
            filtered_img(i, j) = std::round(weighted_sum/normalization);
        }
    }
    return filtered_img;
}

Matrix<unsigned char> despeckle(Matrix<unsigned char> M, std::size_t filtersize)
{
    Matrix<float> M_copy = Matrix<float>(M);
    Matrix<float> f = ones<float>(filtersize);
    f(filtersize/2, filtersize/2) = 0.0f;
    f/=sum(f);
    Matrix<float> mean_img = conv(M_copy, f);
    Matrix<float> dif_img = M_copy-mean_img;
    Matrix<float> stddev_img = conv(M_copy*M_copy, f)-mean_img*mean_img;
    return where(dif_img*dif_img>stddev_img, Matrix<unsigned char>(mean_img), M);
}

Matrix<unsigned char> nagao(Matrix<unsigned char> M, std::size_t filtersize)
{
    Matrix<float> M_copy = Matrix<float>(M);
    Matrix<float> f = average(filtersize);
    Matrix<float> mean_img = conv(M_copy, f);
    Matrix<float> var_img = conv(M_copy*M_copy, f)-mean_img*mean_img;

    Matrix<unsigned char> filtered_img(M.rowNb(), M.colNb());
    for(std::size_t i=0, I=M.rowNb();i<I;++i)
    {
        for(std::size_t j=0, J=M.colNb();j<J;++j)
        {
            float min_var = FLT_MAX;
            for(std::size_t k=0, K=f.rowNb();k<K;++k)
            {
                for(std::size_t l=0, L=f.colNb();l<L;++l)
                {
                    if(i+k>=K/2 && j+l>=L/2 && i+k<I+K/2 && j+l<J+L/2)
                    {
                        std::size_t x = i+k-K/2;
                        std::size_t y = j+l-L/2;
                        if(var_img(x, y)<min_var)
                        {
                            min_var = var_img(x, y);
                            filtered_img(i, j) = std::round(mean_img(x, y));
                        }
                    }
                }
            }
        }
    }
    return filtered_img;
}
