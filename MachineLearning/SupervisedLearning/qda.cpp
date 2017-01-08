#include "qda.hpp"
#include "../../mmath.hpp"
#include "../../linalg.hpp"
#include "../../statistics.hpp"


QDA::QDA()
{
    nb_classes = 1;
}

Matrix<float> QDA::getCenters()
{
    return means;
}

std::vector<Matrix<float> > QDA::getCovariances()
{
    return covars;
}


Matrix<float> QDA::getPriors()
{
    return priors;
}

void QDA::fit(const Matrix<float>& M, const Matrix<std::size_t>& label)
{
    unique_labels = unique(label);
    nb_classes = unique_labels.size();
    nb_features = M.colNb();
    means = zeros<float>(nb_classes, nb_features);
    covars = std::vector<Matrix<float> >(nb_classes, zeros<float>(nb_features));
    priors = zeros<float>(nb_classes, 1);
    std::size_t nb_samples = M.rowNb();

    for(std::size_t i=0;i<nb_samples;++i)
    {
        Matrix<float> cur_sample = M.getRow(i);
        std::size_t cur_label = label(i, 0);
        means.setRow(cur_label, means.getRow(cur_label)+cur_sample);
        ++priors(label(i, 0), 0);
    }
    for(std::size_t k=0;k<nb_classes;++k)
    {
        means.setRow(k, means.getRow(k)/priors(k, 0));
    }

    for(std::size_t i=0;i<nb_samples;++i)
    {
        Matrix<float> cur_sample = M.getRow(i);
        std::size_t cur_label = label(i, 0);
        cur_sample-=means.getRow(cur_label);
        covars[cur_label]+=dot(transpose(cur_sample), cur_sample);
    }
    for(std::size_t k=0;k<nb_classes;++k)
    {
        covars[k]/=priors(k, 0);
    }

    priors/=nb_samples;
}

Matrix<std::size_t> QDA::fit_predict(const Matrix<float>& M, const Matrix<std::size_t>& label)
{
    fit(M, label);
    return predict(M);
}

Matrix<std::size_t> QDA::predict(const Matrix<float>& M)
{
    std::size_t nb_samples = M.rowNb();
    Matrix<float> results_proba = predict_proba(M);
    Matrix<std::size_t> results(nb_samples, 1);
    for(std::size_t i=0;i<nb_samples;++i)
        results(i, 0) = argmax(results_proba.getRow(i))(0, 1);
    return results;
}

Matrix<float> QDA::predict_proba(const Matrix<float>& M)
{
    std::size_t nb_samples = M.rowNb();
    Matrix<float> results(nb_samples, nb_classes);

    for(std::size_t i=0; i<nb_samples;++i)
    {
        Matrix<float> cur_sample = M.getRow(i);
        for(std::size_t j=0;j<nb_classes;++j)
        {
            Matrix<float> inv_covar = pinv(covars[j]);
            Matrix<float> cur_mean = means.getRow(j);
            Matrix<float> X = cur_sample-cur_mean;
            float dist = exp(-0.5f*dot(dot(X, inv_covar), transpose(X)))(0, 0);
            dist/=std::sqrt(det(covars[j]));
            results(i, j) = dist;
        }
    }
    for(std::size_t k=0;k<nb_classes;++k)
        results.setCol(k, results.getCol(k)*priors(k, 0));
    return results;
}
