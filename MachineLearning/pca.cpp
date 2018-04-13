#include "pca.hpp"
#include "../mmath.hpp"
#include "../linalg.hpp"
#include "../statistics.hpp"


PCA::PCA(std::size_t nb_components_arg)
{
    nb_components = nb_components_arg;
}

void PCA::fit(const Matrix<float>& M)
{
    if(nb_components==0)
        nb_components = M.colNb();
    StdSc = StandardScaler(true, true, 2);
    Matrix<float> scaled_mat = StdSc.fit_transform(M);
    covar = dot(transpose(scaled_mat), scaled_mat);
    std::pair<Matrix<float>, Matrix<float> > DP = jacobi(covar);
    eigenvalues = DP.first.getSubmat(0,nb_components,0,nb_components);
    eigenvectors = DP.second.getCols(0, nb_components);
}

Matrix<float> PCA::fit_transform(const Matrix<float>& M)
{
    fit(M);
    return transform(M);
}

Matrix<float> PCA::inverse_transform(const Matrix<float>& M)
{
    Matrix<float> temp = transpose(eigenvectors);
    temp = dot(M, temp);
    return StdSc.inverse_transform(temp);
}

Matrix<float> PCA::transform(const Matrix<float>& M)
{
    Matrix<float> centered_data = StdSc.transform(M);
    Matrix<float> temp = eigenvectors;
    return dot(centered_data, temp);
}
