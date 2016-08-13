#include "knn.hpp"
#include "../../statistics.hpp"
#include "../../linalg.hpp"


KNN::KNN(std::size_t k_arg)
{
    nb_class = 0;
    k = k_arg;
}

void KNN::fit(const Matrix<float>& M, const Matrix<std::size_t>& label)
{
    training_data = M;
    labels = label;
    nb_class = max(label)+1;
}

Matrix<std::size_t> KNN::fit_predict(const Matrix<float>& M, const Matrix<std::size_t>& label)
{
    fit(M, label);
    return predict(M);
}

Matrix<std::size_t> KNN::predict(const Matrix<float>& M)
{
    Matrix<float> results_proba = predict_proba(M);
    Matrix<std::size_t> results(M.rowNb(), 1);
    for(std::size_t i=0;i<M.rowNb();++i)
        results(i, 0) = argmax(results_proba.getRow(i))(0, 1);
    return results;
}

Matrix<float> KNN::predict_proba(const Matrix<float>& M)
{
    Matrix<float> bestresults(M.rowNb(), nb_class);
    Matrix<float> results_on_training_data(training_data.rowNb(), 1);
    for(std::size_t i=0; i<M.rowNb();++i)
    {
        Matrix<float> cur_trainingData = M.getRow(i);
        for(std::size_t j=0;j<training_data.rowNb();++j)
            results_on_training_data(j, 0) = norm(cur_trainingData-training_data.getRow(j));
        Matrix<std::size_t> argmaxmat = argsort(results_on_training_data);
        for(std::size_t j=0;j<k;++j)
            ++bestresults(i, labels(argmaxmat(j, 0), 0));
    }
    return bestresults;
}
