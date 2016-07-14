#include <cfloat>
#include <set>
#include "kmeans.hpp"
#include "../../logical_operators.hpp"
#include "../../mmath.hpp"
#include "../../statistics.hpp"
#include "../../linalg.hpp"


Kmeans::Kmeans(int nb_clusters_arg)
{
    nb_clusters = nb_clusters_arg;
    max_iteration = 20;
}

Matrix<float> Kmeans::getCenters()
{
    return centers;
}

int Kmeans::getMaxIteration()
{
    return max_iteration;
}

void Kmeans::setMaxIteration(int new_max_iteration)
{
    max_iteration = new_max_iteration;
}

void Kmeans::fit(const Matrix<float>& M)
{
    labels = zeros<int>(M.rowNb(), 1);
    //init_centers(M);
    centers = KmeansPlusPlusInit(M, nb_clusters);
    int cpt = 0;
    Matrix<int> old_labels = ones<int>(M.rowNb(), 1);
    while(cpt<max_iteration && any(old_labels!=labels))
    {
        old_labels = labels;
        predict(M);
        update_centers(M);
        ++cpt;
    }
}

Matrix<int> Kmeans::fit_predict(const Matrix<float>& M)
{
    fit(M);
    return labels;
}

Matrix<int> Kmeans::predict(const Matrix<float>& M)
{
    int nb_samples = M.rowNb();
    labels = Matrix<int>(nb_samples, 1);
    for(int i=0;i<nb_samples;++i)
    {
        Matrix<float> distances = dot(ones<float>(nb_clusters, 1), M.getRow(i))-centers;
        distances *= distances;
        distances = axissum(distances, 1);
        labels(i, 0) = argmin(distances)(0, 0);
    }
    return labels;
}

Matrix<float> Kmeans::KmeansPlusPlusInit(const Matrix<float>& M, int nb_clusters)
{
    int n_samples = M.rowNb();
    int n_features = M.colNb();
    Matrix<float> centers = Matrix<float>(nb_clusters, n_features);
    std::vector<int> init_indices;
    init_indices.clear();
    init_indices.push_back(rand()%n_samples);

    Matrix<float> distances(n_samples, 1);

    Matrix<float> cur_center = M.getRow(init_indices[0]);
    for(int i=0;i<n_samples;++i)
        distances(i, 0) = sum(pow(cur_center-M.getRow(i), 2.0f));

    for(int i=1;i<nb_clusters;++i)
    {
        Matrix<float> cum_dist = cumsum(distances);
        float sum_dist = cum_dist(n_samples-1, 0);
        float rand_nb = rand()%1000;
        rand_nb/=1000.0f/sum_dist;
        for(int j=0;j<n_samples;++j)
        {
            if(rand_nb<cum_dist(j, 0))
            {
                init_indices.push_back(j);
                cur_center = M.getRow(j);

                for(int k=0;k<n_samples;++k)
                    distances(k, 0) = std::min(sum(pow(cur_center-M.getRow(k), 2.0f)), distances(k, 0));
                break;
            }
        }
    }

    for(int i=0, I=init_indices.size();i<I;++i)
    {
        centers.setRow(i, M.getRow(init_indices[i]));
    }

    return centers;
}

void Kmeans::update_centers(const Matrix<float>& M)
{
    centers = zeros<float>(nb_clusters, M.colNb());
    Matrix<float> elem_counter = zeros<float>(nb_clusters, 1);
    for(int i=0;i<M.rowNb();++i)
    {
        int cur_label = labels(i, 0);
        centers.setRow(cur_label, centers.getRow(cur_label)+M.getRow(i));
        ++elem_counter(cur_label, 0);
    }
    for(int k=0;k<nb_clusters;++k)
        centers.setRow(k, centers.getRow(k)/elem_counter(k, 0));
}

void Kmeans::init_centers(const Matrix<float>& M)
{
    centers = Matrix<float>(nb_clusters, M.colNb());
    std::set<int> init_indices;
    init_indices.clear();
    while(init_indices.size()<nb_clusters)
    {
        int k = rand()%M.rowNb();
        init_indices.insert(k);
    }
    int i = 0;
    for(std::set<int>::iterator it=init_indices.begin(); it!=init_indices.end(); ++it)
    {
        centers.setRow(i, M.getRow(*it));
        ++i;
    }
}
