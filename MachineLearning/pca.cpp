#include "pca.hpp"
#include "../mmath.hpp"
#include "../linalg.hpp"
#include "../statistics.hpp"


PCA::PCA()
{
}

void PCA::fit(const Matrix<float>& M)
{
    mean = axismean(M, 2);
    Matrix<float> centered_mat = M-dot(ones<float>(M.rowNb(), mean.rowNb()), mean);
    covar = centered_mat;
    covar.transpose();
    covar = dot(covar, centered_mat);
    covar/=M.rowNb();
    std::pair<Matrix<float>, Matrix<float> > DP = Jacobi(covar);
    eigenvalues = DP.first;
    eigenvectors = DP.second;
}

Matrix<float> PCA::fit_transform(const Matrix<float>& M)
{
    fit(M);
    return transform(M);
}

Matrix<float> PCA::inverse_transform(const Matrix<float>& M)
{
    Matrix<float> D = eigenvalues;
    for(std::size_t i = 0, I = D.rowNb();i<I;++i)
    {
        if(std::abs(D(i, i))>10e-7)
            D(i,i)=1.0/D(i,i);
    }
    eigenvectors.transpose();
    Matrix<float> temp = dot(M, dot(sqrt(D),eigenvectors));
    eigenvectors.transpose();
    temp+=dot(ones<float>(M.rowNb(), mean.rowNb()), mean);
    return temp;
}

Matrix<float> PCA::transform(const Matrix<float>& M)
{
    Matrix<float> temp = dot(M-dot(ones<float>(M.rowNb(), mean.rowNb()), mean), eigenvectors);
    return dot(temp , sqrt(eigenvalues));
}
