#include <cfloat>
#include "kmedoids.hpp"
#include "../../statistics.hpp"


Kmedoids::Kmedoids(std::size_t nb_clusters_arg)
{
    nb_clusters = nb_clusters_arg;
    nb_features = -1;
    old_cost = 0.0f;
    cur_cost = 0.0f;
}

Matrix<float> Kmedoids::getCenters()
{
    Matrix<float> centers_values(nb_clusters, nb_features);
    for(std::size_t i=0;i<nb_clusters;++i)
        centers_values.setRow(i, data[centers(i, 0)]);
    return centers_values;
}


void Kmedoids::fit(const Matrix<float>& M)
{

    std::size_t nb_samples = M.rowNb();
    nb_features = M.colNb();
    data = std::vector<Matrix<float> >(nb_samples);
    for(std::size_t i=0;i<nb_samples;++i)
        data[i] = M.getRow(i);
    distances = zeros<float>(nb_samples, nb_samples);
    for(std::size_t i=0;i<nb_samples-1;++i)
    {
        for(std::size_t j=i+1;j<nb_samples;++j)
        {
            distances(i, j) = sum(pow(data[i]-data[j], 2.0f));
            distances(j, i) = distances(i, j);
        }
    }

    centers = rand<std::size_t>(nb_clusters, 1)%nb_samples;

    labels = zeros<std::size_t>(M.rowNb(), 1);
    Matrix<float> min_dist = distances.getCol(centers(0, 0));
    for(std::size_t i=1;i<nb_clusters;++i)
    {
        Matrix<float> cur_dist = distances.getCol(centers(i, 0));
        replace_if(labels, cur_dist<min_dist, i);
        min_dist = min(min_dist, cur_dist);
    }


    old_cost = sum(min_dist);
    bool stop = false;

    while(!stop)
    {
        stop = true;
        for(std::size_t i=0;i<nb_clusters;++i)
        {

            std::size_t temp = centers(i, 0);
            for(std::size_t j=0;j<nb_samples;++j)
            {
                if(labels(j, 0)==i && j!=temp)
                {
                    centers(i, 0) = j;

                    labels = zeros<std::size_t>(nb_samples, 1);
                    Matrix<float> min_dist = distances.getCol(centers(0, 0));
                    for(std::size_t k=1;k<nb_clusters;++k)
                    {
                        Matrix<float> cur_dist = distances.getCol(centers(k, 0));
                        Matrix<bool> cdt = cur_dist<min_dist;
                        replace_if(labels, cdt, k);
                        min_dist = where(cdt, cur_dist, min_dist);
                    }
                    cur_cost = sum(min_dist);
                    if(cur_cost<old_cost)
                    {
                        old_cost = cur_cost;
                        stop = false;
                        break;
                    }
                    else
                        centers(i, 0) = temp;

                }
            }
        }
    }
}

Matrix<std::size_t> Kmedoids::fit_predict(const Matrix<float>& M)
{
    fit(M);
    return predict(M);
}

Matrix<std::size_t> Kmedoids::predict(const Matrix<float>& M)
{
    Matrix<float> centers_values = getCenters();
    std::size_t nb_samples = M.rowNb();
    labels = zeros<std::size_t>(nb_samples, 1);
    for(std::size_t i=0;i<nb_samples;++i)
    {
        Matrix<float> cur_row = M.getRow(i);
        float dist = FLT_MAX;
        for(std::size_t j=0;j<nb_clusters;++j)
        {
            Matrix<float> cur_cl = centers_values.getRow(j);
            float cur_dist = sum(pow(cur_cl-cur_row, 2.0f));
            if(cur_dist<dist)
            {
                dist = cur_dist;
                labels(i, 0) = j;
            }
        }
    }
    return labels;
}
