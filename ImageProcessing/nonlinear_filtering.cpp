#include <float.h>
#include "../statistics.hpp"
#include "linear_filtering.hpp"
#include "nonlinear_filtering.hpp"


Matrix<unsigned char> bilateral(Matrix<unsigned char> M, int filtersize, float range_var, float spatial_var)
{
    Matrix<unsigned char> filtered_img(M.rowNb(), M.colNb());
    for(int i=0, I=M.rowNb();i<I;++i)
    {
        for(int j=0, J=M.colNb();j<J;++j)
        {
            float weighted_sum = 0.0f;
            float normalization = 0.0f;
            for(int k=0;k<filtersize;++k)
            {
                for(int l=0;l<filtersize;++l)
                {
                    int di = k-filtersize/2;
                    int dj = l-filtersize/2;
                    int x = i+di;
                    int y = j+dj;
                    if(x>=0 && y>=0 && x<I && y<J)
                    {
                        float w = std::exp(-((float)M(i, j)-(float)M(x, y))*((float)M(i, j)-(float)M(x, y))/(2.0f*range_var));
                        w*=std::exp(-(di*di+dj*dj)/(2.0f*spatial_var));
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

Matrix<unsigned char> despeckle(Matrix<unsigned char> M, int filtersize)
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

Matrix<unsigned char> nagao(Matrix<unsigned char> M, int filtersize)
{
    Matrix<float> M_copy = Matrix<float>(M);
    Matrix<float> f = average(filtersize);
    Matrix<float> mean_img = conv(M_copy, f);
    Matrix<float> var_img = conv(M_copy*M_copy, f)-mean_img*mean_img;

    Matrix<unsigned char> filtered_img(M.rowNb(), M.colNb());
    for(int i=0;i<M.rowNb();++i)
    {
        for(int j=0;j<M.colNb();++j)
        {
            float min_var = FLT_MAX;
            for(int k=0;k<f.rowNb();++k)
            {
                for(int l=0;l<f.colNb();++l)
                {
                    int x = i+k-f.rowNb()/2;
                    int y = j+l-f.colNb()/2;
                    if(x>=0 && y>=0 && x<M_copy.rowNb() && y<M_copy.colNb() && var_img(x, y)<min_var)
                    {
                        min_var = var_img(x, y);
                        filtered_img(i, j) = std::round(mean_img(x, y));
                    }
                }
            }
        }
    }
    return filtered_img;
}
