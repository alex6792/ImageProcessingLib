#include "gaussian_naive_bayes.hpp"
#include "../../mmath.hpp"
#include "../../linalg.hpp"
#include "../../statistics.hpp"


GaussianNaiveBayes::GaussianNaiveBayes()
{
    nb_class = 0;
}

Matrix<float> GaussianNaiveBayes::getCenters()
{
    return means;
}

Matrix<float> GaussianNaiveBayes::getStddevs()
{
    return stdevs;
}


Matrix<float> GaussianNaiveBayes::getPriors()
{
    return priors;
}

void GaussianNaiveBayes::fit(const Matrix<float>& M, const Matrix<std::size_t>& label)
{
    nb_class = max(label)+1;
    means = Matrix<float>(nb_class, M.colNb());
    stdevs = Matrix<float>(nb_class, M.colNb());
    priors = Matrix<float>(nb_class, 1);
    for(std::size_t i=0;i<M.rowNb();++i)
    {
        for(std::size_t j=0, J=M.colNb();j<J;++j)
        {
            means(label(i, 0), j)+=M(i, j);
            stdevs(label(i, 0), j)+=M(i, j)*M(i, j);
        }
        ++priors(label(i, 0), 0);
    }
    for(std::size_t k=0;k<nb_class;++k)
    {
        means.setRow(k, means.getRow(k)/priors(k, 0));
        stdevs.setRow(k, sqrt(stdevs.getRow(k)/priors(k, 0)-means.getRow(k)*means.getRow(k)));
    }

    priors/=sum(priors);
}

Matrix<std::size_t> GaussianNaiveBayes::fit_predict(const Matrix<float>& M, const Matrix<std::size_t>& label)
{
    fit(M, label);
    return predict(M);
}

Matrix<std::size_t> GaussianNaiveBayes::predict(const Matrix<float>& M)
{
    Matrix<float> results_proba = predict_proba(M);
    Matrix<std::size_t> results(M.rowNb(), 1);
    for(std::size_t i=0;i<M.rowNb();++i)
        results(i, 0) = argmax(results_proba.getRow(i))(0, 1);
    return results;
}

Matrix<float> GaussianNaiveBayes::predict_proba(const Matrix<float>& M)
{
    Matrix<float> results = ones<float>(M.rowNb(), nb_class);
    for(std::size_t i=0; i<M.rowNb();++i)
    {
        for(std::size_t j=0, J=M.colNb();j<J;++j)
        {
            for(std::size_t k=0;k<nb_class;++k)
            {
                results(i, k)*=std::exp(-0.5*std::pow(((M(i, j)-means(k, j))/stdevs(k, j)),2.0))/stdevs(k, j);
            }
        }
    }
    /*for(std::size_t k=0;k<nb_class;++k)
        results.setCol(k, results.getCol(k)*priors(k, 0));*/
    return results;
}
