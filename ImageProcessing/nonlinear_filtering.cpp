#include <float.h>
#include "../statistics.hpp"
#include "linear_filtering.hpp"
#include "nonlinear_filtering.hpp"


Matrix<float> variance_filter(const Matrix<float>& M, std::size_t filtersize)
{
    Matrix<float> M_copy(M);
    Matrix<float> f = average(filtersize);
    Matrix<float> mean_img = conv(M_copy, f);
    Matrix<float> var_img = conv(M_copy*M_copy, f)-mean_img*mean_img;
    return var_img;
}

Matrix<float> geometric_mean_filter(const Matrix<float>& M, std::size_t filtersize)
{
    Matrix<float> filtered_img = log(M);
    Matrix<float> f = average(filtersize);
    filtered_img = conv(filtered_img, f);
    filtered_img = exp(filtered_img);
    return filtered_img;
}

Matrix<float> harmonic_mean_filter(const Matrix<float>& M, std::size_t filtersize)
{
    Matrix<float> filtered_img = 1.0f/M;
    Matrix<float> f = average(filtersize);
    filtered_img = conv(filtered_img, f);
    filtered_img = 1.0f/filtered_img;
    return filtered_img;
}

Matrix<float> quadratic_mean_filter(const Matrix<float>& M, std::size_t filtersize)
{
    Matrix<float> filtered_img = M*M;
    Matrix<float> f = average(filtersize);
    filtered_img = conv(filtered_img, f);
    filtered_img = sqrt(filtered_img);
    return filtered_img;
}

Matrix<float> imguided_filtering(const Matrix<float>& I, const Matrix<float>& G, std::size_t filtersize, float epsilon)
{
    Matrix<float> f = average(filtersize);
    Matrix<float> meanG = conv(G,f);
    Matrix<float> meanI = conv(I,f);
    Matrix<float> corrG = conv(G*G,f);
    Matrix<float> corrIG = conv(I*G,f);
    Matrix<float> varG = corrG-meanG*meanG;
    Matrix<float> covIG = corrIG-meanI*meanG;
    Matrix<float> a = covIG/(varG+epsilon);
    Matrix<float> b = meanI-a*meanG;
    Matrix<float> meana = conv(a,f);
    Matrix<float> meanb = conv(b,f);
    return meana*G+meanb;
}

Matrix<float> imguided_filtering(const Matrix<float>& I, std::size_t filtersize, float epsilon)
{
    Matrix<float> f = average(filtersize);
    Matrix<float> meanI = conv(I,f);
    Matrix<float> corrI = conv(I*I,f);
    Matrix<float> varI = corrI-meanI*meanI;
    Matrix<float> a = 1.0f-epsilon/(I+epsilon);
    Matrix<float> b = (1.0f-a)*meanI;
    Matrix<float> meana = conv(a,f);
    Matrix<float> meanb = conv(b,f);
    return meana*I+meanb;
}

Matrix<unsigned char> bilateral(const Matrix<unsigned char>& M, std::size_t filtersize, float range_var, float spatial_var)
{
    Matrix<unsigned char> filtered_img(M.rowNb(), M.colNb());
    Matrix<float> M_copy(M);

    std::pair<Matrix<float>, Matrix<float> > XY = meshgrid<float>(filtersize);
    Matrix<float>& X = XY.first;
    Matrix<float>& Y = XY.second;
    X-=filtersize/2.0f-0.5f;
    Y-=filtersize/2.0f-0.5f;
    Matrix<float> W_spatial = exp(-X*X-Y*Y)/(2.0f*spatial_var);

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

                        std::size_t x = i+k-filtersize/2;
                        std::size_t y = j+l-filtersize/2;
                        float w = std::exp(-(M_copy(i, j)-M_copy(x, y))*(M_copy(i, j)-M_copy(x, y))/(2.0f*range_var));
                        w*=W_spatial(k, l);
                        weighted_sum+=w*M_copy(x, y);
                        normalization+=w;
                    }
                }
            }
            filtered_img(i, j) = std::round(weighted_sum/normalization);
        }
    }
    return filtered_img;
}

Matrix<unsigned char> despeckle(const Matrix<unsigned char>& M, std::size_t filtersize)
{
    Matrix<float> M_copy = Matrix<float>(M);
    Matrix<float> f = ones<float>(filtersize);
    f(filtersize/2, filtersize/2) = 0.0f;
    f/=sum(f);
    Matrix<float> mean_img = conv(M_copy, f);
    Matrix<float> dif_img = M_copy-mean_img;
    Matrix<float> var_img = conv(M_copy*M_copy, f)-mean_img*mean_img;
    return where(dif_img*dif_img>var_img, Matrix<unsigned char>(mean_img), M);
}

Matrix<unsigned char> nagao(const Matrix<unsigned char>& M, std::size_t filtersize)
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
