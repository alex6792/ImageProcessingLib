#include "fourier.hpp"
#include "../mmath.hpp"
#include "../statistics.hpp"


Matrix<std::complex<float> > FFT(const Matrix<float>& M)
{
    std::complex<float> cste = std::complex<float>(0.0f, -2.0f*PI);
    Matrix<std::complex<float> > M_copy = Matrix<std::complex<float> >(M);
    std::pair<Matrix<std::size_t>, Matrix<std::size_t> > XY = meshgrid<std::size_t>(M.rowNb(), M.colNb());
    Matrix<std::complex<float> > X = Matrix<std::complex<float> >(XY.first);
    Matrix<std::complex<float> > Y = Matrix<std::complex<float> >(XY.second);
    Matrix<std::complex<float> > fft(M.rowNb(), M.colNb());

    for(float i=0, I=M.rowNb();i<I;++i)
    {
        for(float j=0, J=M.colNb();j<J;++j)
        {
            std::complex<float> u = i/I;
            std::complex<float> v = j/J;

            Matrix<std::complex<float> > Data = cste*(X*u+Y*v);
            std::transform(Data.begin(), Data.end(), Data.begin(), [](std::complex<float> arg){return std::exp(arg);});
            Data*=M_copy;
            fft(i, j) = sum(Data);
        }
    }
    return fft/std::complex<float>(fft.size(), 0.0f);
}

Matrix<std::complex<float> > iFFT(const Matrix<std::complex<float> >& M)
{
    std::complex<float> cste = std::complex<float>(0.0f, 2.0f*PI);
    std::pair<Matrix<std::size_t>, Matrix<std::size_t> > XY = meshgrid<std::size_t>(M.rowNb(), M.colNb());
    Matrix<std::complex<float> > X = Matrix<std::complex<float> >(XY.first);
    Matrix<std::complex<float> > Y = Matrix<std::complex<float> >(XY.second);
    Matrix<std::complex<float> > fft(M.rowNb(), M.colNb());

    for(float i=0, I=M.rowNb();i<I;++i)
    {
        for(float j=0, J=M.colNb();j<J;++j)
        {
            std::complex<float> u = i/I;
            std::complex<float> v = j/J;

            Matrix<std::complex<float> > Data = cste*(X*u+Y*v);
            std::transform(Data.begin(), Data.end(), Data.begin(), [](std::complex<float> arg){return std::exp(arg);});
            Data*=M;
            fft(i, j) = sum(Data);
        }
    }
    return fft;
}
