#include <cfloat>
#include "../../linalg.hpp"
#include "../../mmath.hpp"
#include "../../statistics.hpp"
#include "gmm.hpp"
#include "kmeans.hpp"


GMM::GMM(int nb_clusters_arg)
{
    nb_clusters = nb_clusters_arg;
    nb_features = -1;
    max_iteration = 40;
}

Matrix<float> GMM::getCenters()
{
    return centers;
}

Matrix<float> GMM::getStddevs()
{
    return sqrt(vars);
}

Matrix<float> GMM::getWeights()
{
    return weights;
}

void GMM::fit(const Matrix<float>& M)
{
    nb_features = M.colNb();
    vars = full<float>(nb_clusters, nb_features, 1000.0f);
    centers = zeros<float>(nb_clusters, nb_features);
    //Matrix<int> seeds = rand<int>(nb_clusters, 1);
    //for(int i=0;i<nb_clusters;++i)
        //centers.setRow(i, M.getRow(seeds(i, 0)%M.rowNb()));
    centers = Kmeans::KmeansPlusPlusInit(M, nb_clusters);

    weights = ones<float>(nb_clusters, 1)/nb_clusters;
    labels = zeros<int>(M.rowNb(), 1);
    responsibilities = ones<float>(M.rowNb(), nb_clusters)/nb_clusters;

    for(int i=0;i<max_iteration;++i)
    {
        expectation(M);
        maximization(M);
    }
}

Matrix<int> GMM::fit_predict(const Matrix<float>& M)
{
    fit(M);
    return predict(M);
}

Matrix<int> GMM::predict(const Matrix<float>& M)
{
    labels = zeros<int>(M.rowNb(), 1);
    for(int i=0;i<M.rowNb();++i)
    {
        float distance = 0.0f;
        Matrix<float> cur_row = M.getRow(i);
        for(int j=0;j<nb_clusters;++j)
        {
            Matrix<float> cur_mat = cur_row-centers.getRow(j);
            cur_mat*=cur_mat;
            for(int k=0;k<nb_features;++k)
                cur_mat(0, k)/=vars(j, k);

            float dist = weights(j, 0)*std::exp(-0.5*sum(cur_mat))/std::sqrt(prod(vars.getRow(j)));

            if(dist>distance)
            {
                labels(i, 0) = j;
                distance = dist;
            }
        }
    }
    return labels;
}

void GMM::expectation(const Matrix<float>& M)
{
    for(int i=0;i<M.rowNb();++i)
    {
        Matrix<float> cur_row = M.getRow(i);
        for(int j=0;j<nb_clusters;++j)
        {
            Matrix<float> cur_mat = cur_row-centers.getRow(j);
            cur_mat*=cur_mat;
            for(int k=0;k<M.colNb();++k)
                cur_mat(0, k)/=vars(j, k);
            responsibilities(i, j) = weights(j, 0)*std::exp(-0.5*sum(cur_mat))/std::sqrt(prod(vars.getRow(j)));
        }
    }
    Matrix<float> part_sum = axissum(responsibilities, 1);
    for(int i=0;i<M.rowNb();++i)
    {
        for(int j=0;j<nb_clusters;++j)
        {
            responsibilities(i, j)/=part_sum(i, 0);
        }
    }
}

void GMM::maximization(const Matrix<float>& M)
{
    centers = zeros<float>(nb_clusters, nb_features);
    vars = zeros<float>(nb_clusters, nb_features);
    weights = zeros<float>(nb_clusters, 1);
    for(int i=0;i<M.rowNb();++i)
    {
        Matrix<float> cur_mat = M.getRow(i);
        for(int k=0;k<nb_clusters;++k)
        {
            centers.setRow(k, centers.getRow(k)+cur_mat*responsibilities(i, k));
            weights(k, 0) += responsibilities(i, k);
        }
    }
    for(int i=0;i<nb_clusters;++i)
    {
        for(int j=0;j<M.colNb();++j)
            centers(i, j)/=weights(i, 0);
    }
    for(int i=0;i<M.rowNb();++i)
    {
        Matrix<float> cur_row = M.getRow(i);
        for(int k=0;k<nb_clusters;++k)
        {
            Matrix<float> cur_mat = cur_row-centers.getRow(k);
            cur_mat*=cur_mat;
            vars.setRow(k, vars.getRow(k)+cur_mat*responsibilities(i, k));
        }
    }
    for(int i=0;i<nb_clusters;++i)
    {
        for(int j=0;j<nb_features;++j)
            vars(i, j)/=weights(i, 0);
    }
    weights/=sum(weights);
}
