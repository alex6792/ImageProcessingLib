#include "bernoulli_naive_bayes.hpp"
#include "../../linalg.hpp"
#include "../../statistics.hpp"


BernoulliNaiveBayes::BernoulliNaiveBayes()
{
    nb_class = 0;
}

void BernoulliNaiveBayes::fit(const Matrix<bool>& M, const Matrix<std::size_t>& label)
{
    nb_class = max(label)+1;
    proba = Matrix<float>(nb_class, M.colNb());
    priors = Matrix<float>(nb_class, 1);
    for(std::size_t i=0; i<M.rowNb();++i)
    {
        for(std::size_t j=0;j<M.colNb();++j)
        {
            if(M(i, j))
                ++proba(label(i, 0), j);
            ++priors(label(i, 0), 0);
        }
    }
    for(std::size_t k=0;k<nb_class;++k)
        proba.setRow(k, proba.getRow(k)*M.colNb()/priors(k, 0));
    priors/=sum(priors);
}

Matrix<std::size_t> BernoulliNaiveBayes::fit_predict(const Matrix<bool>& M, const Matrix<std::size_t>& label)
{
    fit(M, label);
    return predict(M);
}

Matrix<std::size_t> BernoulliNaiveBayes::predict(const Matrix<bool>& M)
{
    Matrix<float> results_proba = predict_proba(M);
    Matrix<std::size_t> results(M.rowNb(), 1);
    for(std::size_t i=0;i<M.rowNb();++i)
        results(i, 0) = argmax(results_proba.getRow(i))(0, 1);
    return results;
}

Matrix<float> BernoulliNaiveBayes::predict_proba(const Matrix<bool>& M)
{
    Matrix<float> results = ones<float>(M.rowNb(), nb_class);
    for(std::size_t i=0; i<M.rowNb();++i)
    {
        for(std::size_t j=0;j<M.colNb();++j)
        {
            for(std::size_t k=0;k<nb_class;++k)
            {
                if(M(i, j))
                    results(i, k)*=proba(k, j);
                else
                    results(i, k)*=1.0f-proba(k, j);
            }
        }
    }
    /*for(std::size_t k=0;k<nb_class;++k)
        results.setCol(k, results.getCol(k)*priors(k, 0));*/
    return results;
}
