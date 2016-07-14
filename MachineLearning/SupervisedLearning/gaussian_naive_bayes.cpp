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

void GaussianNaiveBayes::fit(const Matrix<float>& M, const Matrix<int>& label)
{
    nb_class = max(label)+1;
    means = Matrix<float>(nb_class, M.colNb());
    stdevs = Matrix<float>(nb_class, M.colNb());
    priors = Matrix<float>(nb_class, 1);
    for(int i=0;i<M.rowNb();++i)
    {
        for(int j=0;j<M.colNb();++j)
        {
            means(label(i, 0), j)+=M(i, j);
            stdevs(label(i, 0), j)+=M(i, j)*M(i, j);
        }
        ++priors(label(i, 0), 0);
    }
    for(int k=0;k<nb_class;++k)
    {
        means.setRow(k, means.getRow(k)/priors(k, 0));
        stdevs.setRow(k, sqrt(stdevs.getRow(k)/priors(k, 0)-means.getRow(k)*means.getRow(k)));
    }

    priors/=sum(priors);
}

Matrix<int> GaussianNaiveBayes::fit_predict(const Matrix<float>& M, const Matrix<int>& label)
{
    fit(M, label);
    return predict(M);
}

Matrix<int> GaussianNaiveBayes::predict(const Matrix<float>& M)
{
    Matrix<float> results_proba = predict_proba(M);
    Matrix<int> results(M.rowNb(), 1);
    for(int i=0;i<M.rowNb();++i)
        results(i, 0) = argmax(results_proba.getRow(i))(0, 1);
    return results;
}

Matrix<float> GaussianNaiveBayes::predict_proba(const Matrix<float>& M)
{
    Matrix<float> results = ones<float>(M.rowNb(), nb_class);
    for(int i=0; i<M.rowNb();++i)
    {
        for(int j=0;j<M.colNb();++j)
        {
            for(int k=0;k<nb_class;++k)
            {
                results(i, k)*=std::exp(-0.5*std::pow(((M(i, j)-means(k, j))/stdevs(k, j)),2.0))/stdevs(k, j);
            }
        }
    }
    /*for(int k=0;k<nb_class;++k)
        results.setCol(k, results.getCol(k)*priors(k, 0));*/
    return results;
}
