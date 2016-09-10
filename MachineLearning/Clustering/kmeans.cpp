#include <cfloat>
#include <set>
#include "kmeans.hpp"
#include "../../logical_operators.hpp"
#include "../../mmath.hpp"
#include "../../statistics.hpp"
#include "../../linalg.hpp"


Kmeans::Kmeans(std::size_t nb_clusters_arg)
{
    nb_clusters = nb_clusters_arg;
    max_iteration = 20;
}

Matrix<float> Kmeans::getCenters()
{
    return centers;
}

std::size_t Kmeans::getMaxIteration()
{
    return max_iteration;
}

void Kmeans::setMaxIteration(std::size_t new_max_iteration)
{
    max_iteration = new_max_iteration;
}

void Kmeans::fit(const Matrix<float>& M)
{
    std::size_t nb_samples = M.rowNb();
    nb_features = M.colNb();
    labels = zeros<std::size_t>(nb_samples, 1);
    //init_centers(M);
    centers = KmeansPlusPlusCenters(M, nb_clusters);
    std::size_t cpt = 0;
    Matrix<std::size_t> old_labels = ones<std::size_t>(nb_samples, 1);
    while(cpt<max_iteration && any(old_labels!=labels))
    {
        old_labels = labels;
        predict(M);
        update_centers(M);
        ++cpt;
    }
}

Matrix<std::size_t> Kmeans::fit_predict(const Matrix<float>& M)
{
    fit(M);
    return labels;
}

Matrix<std::size_t> Kmeans::predict(const Matrix<float>& M)
{
    std::size_t nb_samples = M.rowNb();
    labels = Matrix<std::size_t>(nb_samples, 1);
    for(std::size_t i=0;i<nb_samples;++i)
    {
        Matrix<float> distances = dot(ones<float>(nb_clusters, 1), M.getRow(i))-centers;
        distances *= distances;
        distances = axissum(distances, 1);
        labels(i, 0) = argmin(distances)(0, 0);
    }
    return labels;
}

Matrix<float> Kmeans::KmeansPlusPlusCenters(const Matrix<float>& M, std::size_t nb_clusters)
{
    std::size_t n_samples = M.rowNb();
    std::size_t n_features = M.colNb();
    Matrix<float> centers = Matrix<float>(nb_clusters, n_features);
    std::vector<std::size_t> init_indices;
    init_indices.clear();
    init_indices.push_back(rand()%n_samples);

    Matrix<float> distances(n_samples, 1);

    Matrix<float> cur_center = M.getRow(init_indices[0]);
    for(std::size_t i=0;i<n_samples;++i)
        distances(i, 0) = sum(pow(cur_center-M.getRow(i), 2.0f));

    for(std::size_t i=1;i<nb_clusters;++i)
    {
        Matrix<float> cum_dist = cumsum(distances);
        float sum_dist = cum_dist(n_samples-1, 0);
        float rand_nb = rand()%1000;
        rand_nb/=1000.0f/sum_dist;
        for(std::size_t j=0;j<n_samples;++j)
        {
            if(rand_nb<cum_dist(j, 0))
            {
                init_indices.push_back(j);
                cur_center = M.getRow(j);

                for(std::size_t k=0;k<n_samples;++k)
                    distances(k, 0) = std::min(sum(pow(cur_center-M.getRow(k), 2.0f)), distances(k, 0));
                break;
            }
        }
    }

    for(std::size_t i=0, I=init_indices.size();i<I;++i)
    {
        centers.setRow(i, M.getRow(init_indices[i]));
    }

    return centers;
}

Matrix<size_t> Kmeans::KmeansPlusPlusCentersIndices(const Matrix<float>& M, std::size_t nb_clusters)
{
    std::size_t n_samples = M.rowNb();
    Matrix<size_t> centers = Matrix<size_t>(nb_clusters, 1);
    std::vector<std::size_t> init_indices;
    init_indices.clear();
    init_indices.push_back(rand()%n_samples);

    Matrix<float> distances(n_samples, 1);

    Matrix<float> cur_center = M.getRow(init_indices[0]);
    for(std::size_t i=0;i<n_samples;++i)
        distances(i, 0) = sum(pow(cur_center-M.getRow(i), 2.0f));

    for(std::size_t i=1;i<nb_clusters;++i)
    {
        Matrix<float> cum_dist = cumsum(distances);
        float sum_dist = cum_dist(n_samples-1, 0);
        float rand_nb = rand()%1000;
        rand_nb/=1000.0f/sum_dist;
        for(std::size_t j=0;j<n_samples;++j)
        {
            if(rand_nb<cum_dist(j, 0))
            {
                init_indices.push_back(j);
                cur_center = M.getRow(j);

                for(std::size_t k=0;k<n_samples;++k)
                    distances(k, 0) = std::min(sum(pow(cur_center-M.getRow(k), 2.0f)), distances(k, 0));
                break;
            }
        }
    }

    for(std::size_t i=0, I=init_indices.size();i<I;++i)
    {
        //centers.setRow(i, M.getRow(init_indices[i]));
        centers(i, 0) = init_indices[i];
    }

    return centers;
}

void Kmeans::update_centers(const Matrix<float>& M)
{
    centers = zeros<float>(nb_clusters, nb_features);
    Matrix<float> elem_counter = zeros<float>(nb_clusters, 1);
    for(std::size_t i=0;i<M.rowNb();++i)
    {
        std::size_t cur_label = labels(i, 0);
        centers.setRow(cur_label, centers.getRow(cur_label)+M.getRow(i));
        ++elem_counter(cur_label, 0);
    }
    for(std::size_t k=0;k<nb_clusters;++k)
        centers.setRow(k, centers.getRow(k)/elem_counter(k, 0));
}

void Kmeans::init_centers(const Matrix<float>& M)
{
    centers = Matrix<float>(nb_clusters, nb_features);
    std::set<std::size_t> init_indices;
    init_indices.clear();
    std::size_t H = M.rowNb();
    while(init_indices.size()<nb_clusters)
    {
        std::size_t k = rand()%H;
        init_indices.insert(k);
    }
    std::size_t i = 0;
    for(std::set<std::size_t>::iterator it=init_indices.begin(); it!=init_indices.end(); ++it)
    {
        centers.setRow(i, M.getRow(*it));
        ++i;
    }
}
