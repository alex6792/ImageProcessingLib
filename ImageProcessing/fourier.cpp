#include "fourier.hpp"
#include "../mmath.hpp"
#include "../statistics.hpp"


Matrix<std::complex<float> > FFT(const Matrix<float>& M)
{
    std::size_t H = M.rowNb(), W = M.colNb();
    std::complex<float> cste = std::complex<float>(0.0f, -2.0f*PI);
    Matrix<std::complex<float> > M_copy = Matrix<std::complex<float> >(M);
    std::pair<Matrix<std::size_t>, Matrix<std::size_t> > XY = meshgrid<std::size_t>(H, W);
    Matrix<std::complex<float> > X = Matrix<std::complex<float> >(XY.first);
    Matrix<std::complex<float> > Y = Matrix<std::complex<float> >(XY.second);
    Matrix<std::complex<float> > fft(H, W);
    X*=cste;
    Y*=cste;

    for(float i=0;i<H;++i)
    {
        for(float j=0;j<W;++j)
        {
            std::complex<float> u = i/H;
            std::complex<float> v = j/W;
            Matrix<std::complex<float> > Data = exp(X*u+Y*v);
            Data*=M_copy;
            fft(i, j) = sum(Data);
        }
    }
    return fft/std::complex<float>(fft.size(), 0.0f);
}

Matrix<std::complex<float> > iFFT(const Matrix<std::complex<float> >& M)
{
    std::size_t H = M.rowNb(), W = M.colNb();
    std::complex<float> cste = std::complex<float>(0.0f, 2.0f*PI);
    std::pair<Matrix<std::size_t>, Matrix<std::size_t> > XY = meshgrid<std::size_t>(M.rowNb(), M.colNb());
    Matrix<std::complex<float> > X = Matrix<std::complex<float> >(XY.first);
    Matrix<std::complex<float> > Y = Matrix<std::complex<float> >(XY.second);
    Matrix<std::complex<float> > fft(H, W);
    X*=cste;
    Y*=cste;

    for(float i=0;i<H;++i)
    {
        for(float j=0;j<W;++j)
        {
            std::complex<float> u = i/H;
            std::complex<float> v = j/W;
            Matrix<std::complex<float> > Data = exp(X*u+Y*v);
            Data*=M;
            fft(i, j) = sum(Data);
        }
    }
    return fft;
}
