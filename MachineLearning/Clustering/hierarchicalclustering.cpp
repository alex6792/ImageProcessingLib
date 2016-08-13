#include <cfloat>
#include "hierarchicalclustering.hpp"
#include "../../statistics.hpp"


HierarchicalClustering::HierarchicalClustering(std::size_t nb_cluster_arg)
{
    nb_clusters = nb_cluster_arg;
    nb_features = 0;
}



void HierarchicalClustering::fit(const Matrix<float>& M)
{
    data = M;
    std::size_t nb_sample = M.rowNb();
    nb_features = M.colNb();
    std::vector<Matrix<float> > rows(nb_sample);
    for(std::size_t i=0;i<nb_sample;++i)
        rows[i] = M.getRow(i);
    Matrix<float> dist = full(nb_sample, nb_sample, FLT_MAX);
    for(std::size_t i=0;i<nb_sample-1;++i)
    {
        for(std::size_t j=i+1;j<nb_sample;++j)
            dist(i, j) = sum(pow(rows[i]-rows[j], 2.0f));
    }
    labels = arange<std::size_t>(0, nb_sample);

    Matrix<std::size_t> sorted_dist = argsort(dist);

    std::size_t cur_nb_clusters = nb_sample;
    std::size_t cur_idx = 0;
    while(cur_nb_clusters>nb_clusters)
    {
        std::size_t cluster1_idx = sorted_dist(cur_idx, 0);
        std::size_t cluster2_idx = sorted_dist(cur_idx, 1);
        std::size_t cluster1 = labels(cluster1_idx, 0);
        std::size_t cluster2 = labels(cluster2_idx, 0);
        if(cluster1!=cluster2)
        {
            replace(labels, cluster1, cluster2);
            --cur_nb_clusters;
        }
        ++cur_idx;
    }
    Matrix<std::size_t> unique_labels = unique(labels);
    for(std::size_t i=0;i<unique_labels.size();++i)
        replace(labels, unique_labels(i, 0), i);
}

Matrix<std::size_t> HierarchicalClustering::fit_predict(const Matrix<float>& M)
{
    fit(M);
    return labels;
}

Matrix<std::size_t> HierarchicalClustering::predict(const Matrix<float>& M)
{
    std::size_t nb_samples = M.rowNb();
    Matrix<std::size_t> returned_labels(nb_samples, 1);
    for(std::size_t i=0;i<nb_samples;++i)
    {
        float min_dist = FLT_MAX;
        Matrix<float> cur_sample = M.getRow(i);
        for(std::size_t j=0;j<data.rowNb();++j)
        {
            Matrix<float> dif = cur_sample-data.getRow(j);
            dif*=dif;
            float dist = sum(dif);
            if(dist<min_dist)
            {
                min_dist = dist;
                returned_labels(i, 0) = labels(j, 0);
            }
        }
    }
    return returned_labels;
}
