#include <cfloat>
#include "affinitypropagation.hpp"
#include "../../linalg.hpp"
#include "../../statistics.hpp"


AffinityPropagation::AffinityPropagation()
{
    max_iteration = 40;
    nb_clusters = -1;
    nb_features = -1;
}

Matrix<float> AffinityPropagation::getAvailability()
{
    return availability;
}

Matrix<float> AffinityPropagation::getResponsibility()
{
    return responsibility;
}

Matrix<float> AffinityPropagation::getCenters()
{
    return centers;
}

std::size_t AffinityPropagation::getMaxIteration()
{
    return max_iteration;
}

void AffinityPropagation::setMaxIteration(std::size_t max_iteration_arg)
{
    max_iteration = max_iteration_arg;
}

void AffinityPropagation::fit(const Matrix<float>& M)
{
    std::size_t nb_samples = M.rowNb();
    nb_features = M.colNb();

    std::vector<Matrix<float> > rows(nb_samples);
    for(std::size_t i=0;i<nb_samples;++i)
        rows[i] = M.getRow(i);
    Matrix<float> dist = zeros<float>(nb_samples, nb_samples);
    std::vector<float> distances;
    distances.clear();
    for(std::size_t i=0;i<nb_samples-1;++i)
    {
        for(std::size_t j=i+1;j<nb_samples;++j)
        {
            dist(i, j) = -sum(pow(rows[i]-rows[j], 2.0f));
            dist(j, i) = dist(i, j);
            distances.push_back(dist(i, j));
        }
    }
    std::sort(distances.begin(), distances.end());
    float median;
    std::size_t size = distances.size();
    median= distances[(size-1)/2];
    for(std::size_t i=0;i<nb_samples;++i)
        dist(i, i) = median;


    labels = arange<std::size_t>(0, nb_samples);
    float damping_factor = 0.7;
    responsibility = zeros<float>(nb_samples, nb_samples);
    availability = zeros<float>(nb_samples, nb_samples);
    Matrix<float> E;

    std::size_t it=0;
    while(it++<max_iteration || nb_clusters==0)
    {
        Matrix<float> old_R = responsibility;
        Matrix<float> maximum = full<float>(nb_samples, nb_samples, -FLT_MAX);
        Matrix<float> temp = availability+dist;
        float switch_value = -FLT_MAX;
        for(std::size_t i=0;i<nb_samples;++i)
        {
            Matrix<float> temp_row = temp.getRow(i);
            Matrix<std::size_t> arg_max = argmax(temp_row);
            std::size_t& x = arg_max(0, 0), y = arg_max(0, 1);
            maximum.setRow(i, full<float>(1, nb_samples, temp_row(x, y)));
            temp_row(x, y) = -FLT_MAX;
            maximum(i, y) = max(temp_row);
        }
        responsibility = dist - maximum;
        responsibility = responsibility*(1-damping_factor)+damping_factor*old_R;

        Matrix<float> old_A = availability;
        Matrix<float> max_respo_zero = max(responsibility, zeros<float>(nb_samples));
        for(std::size_t k=0;k<nb_samples;++k)
        {
            float sum_m = sum(max_respo_zero.getCol(k))-max_respo_zero(k, k);
            for(std::size_t i=0;i<nb_samples;++i)
            {
                if(i!=k)
                    availability(i, k) = std::min(0.0f, responsibility(k, k)+sum_m-max_respo_zero(i, k));
                else
                    availability(k, k) = sum_m;
            }
        }
        availability = availability*(1-damping_factor)+damping_factor*old_A;
        E = availability.getDiag()+responsibility.getDiag();
        nb_clusters = count_nonzero(E>0);
    }

    Matrix<std::size_t> cluster_idx = argwhere(E>0);
    centers = Matrix<float>(nb_clusters, nb_features);
    for(std::size_t i=0;i<nb_clusters;++i)
        centers.setRow(i, M.getRow(cluster_idx(i, 0)));

}

Matrix<std::size_t> AffinityPropagation::fit_predict(const Matrix<float>& M)
{
    fit(M);
    return predict(M);
}

Matrix<std::size_t> AffinityPropagation::predict(const Matrix<float>& M)
{
    std::size_t nb_sampless = M.rowNb();
    labels = Matrix<std::size_t>(nb_sampless, 1);
    for(std::size_t i=0;i<M.rowNb();++i)
    {
        Matrix<float> distances = dot(ones<float>(nb_clusters, 1), M.getRow(i))-centers;
        distances *= distances;
        distances = axissum(distances, 1);
        labels(i, 0) = argmin(distances)(0, 0);
    }
    return labels;
}

