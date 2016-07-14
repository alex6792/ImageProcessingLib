#include <cfloat>
#include "hierarchicalclustering.hpp"
#include "../../statistics.hpp"


HierarchicalClustering::HierarchicalClustering(int nb_cluster_arg)
{
    nb_clusters = nb_cluster_arg;
}



void HierarchicalClustering::fit(const Matrix<float>& M)
{
    data = M;
    int nb_sample = M.rowNb();
    nb_features = M.colNb();
    std::vector<Matrix<float> > rows(nb_sample);
    for(int i=0;i<nb_sample;++i)
        rows[i] = M.getRow(i);
    Matrix<float> dist = full(nb_sample, nb_sample, FLT_MAX);
    for(int i=0;i<nb_sample-1;++i)
    {
        for(int j=i+1;j<nb_sample;++j)
            dist(i, j) = sum(pow(rows[i]-rows[j], 2.0f));
    }
    labels = arange<int>(0, nb_sample);

    Matrix<int> sorted_dist = argsort(dist);

    int cur_nb_clusters = nb_sample;
    int cur_idx = 0;
    while(cur_nb_clusters>nb_clusters)
    {
        int cluster1_idx = sorted_dist(cur_idx, 0);
        int cluster2_idx = sorted_dist(cur_idx, 1);
        int cluster1 = labels(cluster1_idx, 0);
        int cluster2 = labels(cluster2_idx, 0);
        if(cluster1!=cluster2)
        {
            replace(labels, cluster1, cluster2);
            --cur_nb_clusters;
        }
        ++cur_idx;
    }
    Matrix<int> unique_labels = unique(labels);
    for(int i=0;i<unique_labels.size();++i)
        replace(labels, unique_labels(i, 0), i);
}

Matrix<int> HierarchicalClustering::fit_predict(const Matrix<float>& M)
{
    fit(M);
    return labels;
}

Matrix<int> HierarchicalClustering::predict(const Matrix<float>& M)
{
    int nb_samples = M.rowNb();
    Matrix<int> returned_labels(nb_samples, 1);
    for(int i=0;i<nb_samples;++i)
    {
        float min_dist = FLT_MAX;
        Matrix<float> cur_sample = M.getRow(i);
        for(int j=0;j<data.rowNb();++j)
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
