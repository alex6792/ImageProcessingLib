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


Matrix<std::complex<float> > FFTshift(const Matrix<std::complex<float> >& M)
{
    std::size_t H = M.rowNb(), W = M.colNb();
    Matrix<std::complex<float>> res(H, W);

    if(H>1 && W>1)
        res.setSubmat(0, 0, M.getSubmat((H+1)/2, H, (W+1)/2, W));
    if(W>1)
        res.setSubmat(H-(H+1)/2, 0, M.getSubmat(0, (H+1)/2, (W+1)/2, W));
    if(H>1)
        res.setSubmat(0, W-(W+1)/2, M.getSubmat((H+1)/2, H, 0, (W+1)/2));

        res.setSubmat(H-(H+1)/2, W-(W+1)/2, M.getSubmat(0, (H+1)/2, 0, (W+1)/2));

    return res;
}

Matrix<float> FFTfreq(std::size_t n, float d)
{
    Matrix<float> freqs(n, 1);
    if(n%2==0)
    {
        freqs.setSubmat(0, 0, arange<float>(0, n/2));
        freqs.setSubmat(n/2, 0, arange<float>(-float(n/2), 0.0f));
    }
    else
    {
        freqs.setSubmat(0, 0, arange<float>(0, (n+1)/2));
        freqs.setSubmat((n+1)/2, 0, arange<float>(-float((n-1)/2), 0.0f));
    }
    return freqs/(d*n);
}
