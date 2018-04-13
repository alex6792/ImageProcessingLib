#include "nnmf.hpp"
#include "../mmath.hpp"
#include "../linalg.hpp"
#include "../optimization.hpp"
#include "../statistics.hpp"
#include <cfloat>


NNMF::NNMF(std::size_t nb_components_arg)
{
    nb_components = nb_components_arg;
}

void NNMF::fit(const Matrix<float>& M)
{
    if(nb_components==0)
        nb_components = M.colNb();
    W = rand<float>(M.rowNb(), nb_components);
    H = rand<float>(nb_components, M.colNb());
    float alpha = sum(M*M)/sum(dot(W,H));
    W*=std::sqrt(std::abs(alpha));
    H*=alpha/std::sqrt(std::abs(alpha));

    Matrix<float> Zh = zeros<float>(H.rowNb(),H.colNb());
    Matrix<float> Ubh = full<float>(H.rowNb(),H.colNb(),FLT_MAX);
    Matrix<float> Zw = zeros<float>(W.colNb(),W.rowNb());
    Matrix<float> Ubw = full<float>(W.colNb(),W.rowNb(),FLT_MAX);
    std::size_t cpt = 0;
    while(cpt<100)
    {
        ++cpt;
        Matrix<float> err = dot(W,H)-M;
        Matrix<float> gradW = dot(err,transpose(H));
        if(sum(gradW*gradW)<10e-10)
            break;
        Matrix<float> gradH = dot(transpose(W),err);
        if(sum(gradH*gradH)<10e-10)
            break;

        W = transpose(solve(transpose(H), transpose(M), transpose(W), Zw, Ubw));//
        H = solve(W,M,H,Zh,Ubh);//
    }
}

Matrix<float> NNMF::fit_transform(const Matrix<float>& M)
{
    fit(M);
    return W;
}

Matrix<float> NNMF::inverse_transform(const Matrix<float>& M)
{
    return dot(M,H);
}

Matrix<float> NNMF::transform(const Matrix<float>& M)
{
    //find W such that M = WH ie M' = H'W'
    Matrix<float> Z = zeros<float>(H.rowNb(),M.rowNb());
    Matrix<float> Ub = full<float>(H.rowNb(),M.rowNb(),FLT_MAX);
    return transpose(solve(transpose(H), transpose(M), Z, Z, Ub));
}

