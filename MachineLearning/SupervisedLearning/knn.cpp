#include "knn.hpp"
#include "../../statistics.hpp"
#include "../../linalg.hpp"


KNN::KNN(int k_arg)
{
    nb_class = 0;
    k = k_arg;
}

void KNN::fit(const Matrix<float>& M, const Matrix<int>& label)
{
    training_data = M;
    labels = label;
    nb_class = max(label)+1;
}

Matrix<int> KNN::fit_predict(const Matrix<float>& M, const Matrix<int>& label)
{
    fit(M, label);
    return predict(M);
}

Matrix<int> KNN::predict(const Matrix<float>& M)
{
    Matrix<float> results_proba = predict_proba(M);
    Matrix<int> results(M.rowNb(), 1);
    for(int i=0;i<M.rowNb();++i)
        results(i, 0) = argmax(results_proba.getRow(i))(0, 1);
    return results;
}

Matrix<float> KNN::predict_proba(const Matrix<float>& M)
{
    Matrix<float> bestresults(M.rowNb(), nb_class);
    Matrix<float> results_on_training_data(training_data.rowNb(), 1);
    for(int i=0; i<M.rowNb();++i)
    {
        Matrix<float> cur_trainingData = M.getRow(i);
        for(int j=0;j<training_data.rowNb();++j)
            results_on_training_data(j, 0) = norm(cur_trainingData-training_data.getRow(j));
        Matrix<int> argmaxmat = argsort(results_on_training_data);
        for(int j=0;j<k;++j)
            ++bestresults(i, labels(argmaxmat(j, 0), 0));
    }
    return bestresults;
}
