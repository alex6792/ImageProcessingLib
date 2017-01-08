#include <cfloat>
#include "../../linalg.hpp"
#include "../../mmath.hpp"
#include "../../statistics.hpp"
#include "gmm.hpp"
#include "kmeans.hpp"


GMM::GMM(std::size_t nb_clusters_arg, std::string covar_type_arg)
{
    nb_clusters = nb_clusters_arg;
    covar_type = covar_type_arg;
    nb_features = -1;
    max_iteration = 100;
}

Matrix<float> GMM::getCenters()
{
    return centers;
}

std::vector<Matrix<float> > GMM::getCovariances()
{
    return covars;
}

Matrix<float> GMM::getWeights()
{
    return weights;
}

void GMM::setMaxIteration(std::size_t new_max_iteration)
{
    max_iteration = new_max_iteration;
}

void GMM::fit(const Matrix<float>& M)
{
    nb_features = M.colNb();

    if(!covar_type.compare("full"))
        covars = std::vector<Matrix<float> >(nb_clusters, 10.0f*id<float>(nb_features, nb_features));
    else if(!covar_type.compare("tied"))
        covars = std::vector<Matrix<float> >(1, 10.0f*id<float>(nb_features, nb_features));
    else if(!covar_type.compare("diag"))
        covars = std::vector<Matrix<float> >(nb_clusters, full<float>(1, nb_features, 10.0f));
    else if(!covar_type.compare("spherical"))
        covars = std::vector<Matrix<float> >(nb_clusters, full<float>(1, 1, 10.0f));

    centers = Kmeans::KmeansPlusPlusCenters(M, nb_clusters);
    weights = ones<float>(nb_clusters, 1)/nb_clusters;// mixture weights
    labels = zeros<std::size_t>(M.rowNb(), 1);
    responsibilities = ones<float>(M.rowNb(), nb_clusters)/nb_clusters;//membership weights

    float old_log_likelihood = 0.0f;
    float new_log_likelihood = -10.0f;
    std::size_t it = 0;
    while((it<max_iteration)&&(std::abs(new_log_likelihood-old_log_likelihood)>10e-9))
    {
        expectation(M);
        maximization(M);
        old_log_likelihood = new_log_likelihood;
        new_log_likelihood = compute_log_likelihood(M);
    }
}

Matrix<std::size_t> GMM::fit_predict(const Matrix<float>& M)
{
    fit(M);
    return predict(M);
}

Matrix<std::size_t> GMM::predict(const Matrix<float>& M)
{
    std::size_t nb_samples = M.rowNb();
    labels = zeros<std::size_t>(nb_samples, 1);
    for(std::size_t i=0;i<nb_samples;++i)
    {
        float max_proba = 0.0f;
        Matrix<float> cur_row = M.getRow(i);
        for(std::size_t j=0;j<nb_clusters;++j)
        {
            Matrix<float> cur_mat = cur_row-centers.getRow(j);

            Matrix<float> cur_covar;
            if(!covar_type.compare("full"))
            {
                cur_covar = covars[j];
            }
            else if(!covar_type.compare("tied"))
            {
                cur_covar = covars[0];
            }
            else if(!covar_type.compare("diag"))
            {
                cur_covar = Matrix<float>(nb_features);
                cur_covar.setDiag(covars[j]);
            }
            else if(!covar_type.compare("spherical"))
            {
                cur_covar = covars[j](0,0)*id<float>(nb_features);
            }
            Matrix<float> inv_covar = pinv(cur_covar);

            float proba = weights(j, 0)*exp(-0.5f*dot(dot(cur_mat, inv_covar), transpose(cur_mat)))(0, 0);
            proba/=std::sqrt(det(cur_covar));

            if(proba>max_proba)
            {
                labels(i, 0) = j;
                max_proba = proba;
            }
        }
    }
    return labels;
}

void GMM::expectation(const Matrix<float>& M)
{
    std::size_t nb_samples = M.rowNb();
    Matrix<float> cur_covar;
    for(std::size_t j=0;j<nb_clusters;++j)
    {
        Matrix<float> cur_center = centers.getRow(j);

        if(!covar_type.compare("full"))
        {
            cur_covar = covars[j];
        }
        else if(!covar_type.compare("tied"))
        {
            cur_covar = covars[0];
        }
        else if(!covar_type.compare("diag"))
        {
            cur_covar = zeros<float>(nb_features, nb_features);
            cur_covar.setDiag(covars[j]);
        }
        else if(!covar_type.compare("spherical"))
        {
            cur_covar = covars[j](0,0)*id<float>(nb_features);
        }
        Matrix<float> inv_covar = pinv(cur_covar);

        for(std::size_t i=0;i<nb_samples;++i)
        {
            Matrix<float> cur_mat = M.getRow(i)-cur_center;
            float dist = weights(j, 0)*exp(-0.5f*dot(dot(cur_mat, inv_covar), transpose(cur_mat)))(0, 0);
            dist/=std::sqrt(det(cur_covar));
            responsibilities(i, j) = dist;
        }
    }
    Matrix<float> part_sum = axissum(responsibilities, 1);
    for(std::size_t i=0;i<nb_samples;++i)
    {
        for(std::size_t j=0;j<nb_clusters;++j)
        {
            responsibilities(i, j)/=part_sum(i, 0);
        }
    }
}

void GMM::maximization(const Matrix<float>& M)
{
    std::size_t nb_samples = M.rowNb();
    centers = zeros<float>(nb_clusters, nb_features);
    weights = zeros<float>(nb_clusters, 1);
    for(std::size_t i=0;i<nb_samples;++i)
    {
        Matrix<float> cur_mat = M.getRow(i);
        for(std::size_t k=0;k<nb_clusters;++k)
        {
            centers.setRow(k, centers.getRow(k)+cur_mat*responsibilities(i, k));
            weights(k, 0) += responsibilities(i, k);
        }
    }
    for(std::size_t i=0;i<nb_clusters;++i)
    {
        for(std::size_t j=0;j<nb_features;++j)
            centers(i, j)/=weights(i, 0);
    }

    if(!covar_type.compare("full"))
    {
        covars = std::vector<Matrix<float> >(nb_clusters, zeros<float>(nb_features));

        for(std::size_t i=0;i<nb_samples;++i)
        {
            Matrix<float> cur_row = M.getRow(i);
            for(std::size_t k=0;k<nb_clusters;++k)
            {
                Matrix<float> cur_mat = cur_row-centers.getRow(k);
                cur_mat=dot(transpose(cur_mat), cur_mat);
                covars[k]+=cur_mat*responsibilities(i, k);
            }
        }

        for(std::size_t i=0;i<nb_clusters;++i)
        {
            covars[i]/=weights(i, 0);
        }
    }
    else if(!covar_type.compare("tied"))
    {
        covars = std::vector<Matrix<float> >(1, zeros<float>(nb_features));

        for(std::size_t i=0;i<nb_samples;++i)
        {
            Matrix<float> cur_row = M.getRow(i);
            for(std::size_t k=0;k<nb_clusters;++k)
            {
                Matrix<float> cur_mat = cur_row-centers.getRow(k);
                cur_mat=dot(transpose(cur_mat), cur_mat);
                covars[0]+=cur_mat*responsibilities(i, k);
            }
        }

        covars[0]/=sum(weights);
    }
    else if(!covar_type.compare("diag"))
    {
        covars = std::vector<Matrix<float> >(nb_clusters, zeros<float>(1, nb_features));

        for(std::size_t i=0;i<nb_samples;++i)
        {
            Matrix<float> cur_row = M.getRow(i);
            for(std::size_t k=0;k<nb_clusters;++k)
            {
                Matrix<float> cur_mat = cur_row-centers.getRow(k);
                cur_mat*=cur_mat;
                covars[k]+=cur_mat*responsibilities(i, k);
            }
        }
        for(std::size_t i=0;i<nb_clusters;++i)
        {
            covars[i]/=weights(i, 0);
        }
    }
    else if(!covar_type.compare("spherical"))
    {
        covars = std::vector<Matrix<float> >(nb_clusters, zeros<float>(1, nb_features));

        for(std::size_t i=0;i<nb_samples;++i)
        {
            Matrix<float> cur_row = M.getRow(i);
            for(std::size_t k=0;k<nb_clusters;++k)
            {
                Matrix<float> cur_mat = cur_row-centers.getRow(k);
                cur_mat*=cur_mat;
                covars[k]+=cur_mat*responsibilities(i, k);
            }
        }
        for(std::size_t i=0;i<nb_clusters;++i)
        {
            covars[i]/=weights(i, 0);
            covars[i] = geometric_mean(covars[i])*id<float>(1);
        }
    }

    weights/=sum(weights);
}

float GMM::compute_log_likelihood(const Matrix<float>& M)
{
    float log_likelihood = 0.0f;
    std::size_t nb_samples = M.rowNb();
    Matrix<float> cur_covar;
    for(std::size_t j=0;j<nb_clusters;++j)
    {
        Matrix<float> cur_center = centers.getRow(j);

        if(!covar_type.compare("full"))
        {
            cur_covar = covars[j];
        }
        else if(!covar_type.compare("tied"))
        {
            cur_covar = covars[0];
        }
        else if(!covar_type.compare("diag"))
        {
            cur_covar = zeros<float>(nb_features, nb_features);
            cur_covar.setDiag(covars[j]);
        }
        else if(!covar_type.compare("spherical"))
        {
            cur_covar = covars[j](0,0)*id<float>(nb_features);
        }
        Matrix<float> inv_covar = pinv(cur_covar);

        for(std::size_t i=0;i<nb_samples;++i)
        {
            Matrix<float> cur_mat = M.getRow(i)-cur_center;
            float dist = std::log(weights(j, 0))-0.5f*dot(dot(cur_mat, inv_covar), transpose(cur_mat))(0, 0);
            dist-=std::log(std::sqrt(det(cur_covar)));
            log_likelihood += dist;
        }
    }
    return log_likelihood;
}
