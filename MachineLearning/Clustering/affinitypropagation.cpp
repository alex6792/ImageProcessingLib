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

int AffinityPropagation::getMaxIteration()
{
    return max_iteration;
}

void AffinityPropagation::setMaxIteration(int max_iteration_arg)
{
    max_iteration = max_iteration_arg;
}

void AffinityPropagation::fit(const Matrix<float>& M)
{
    int nb_sample = M.rowNb();
    nb_features = M.colNb();

    std::vector<Matrix<float> > rows(nb_sample);
    for(int i=0;i<nb_sample;++i)
        rows[i] = M.getRow(i);
    Matrix<float> dist = zeros<float>(nb_sample, nb_sample);
    std::vector<float> distances;
    distances.clear();
    for(int i=0;i<nb_sample-1;++i)
    {
        for(int j=i+1;j<nb_sample;++j)
        {
            dist(i, j) = -sum(pow(rows[i]-rows[j], 2.0f));
            dist(j, i) = dist(i, j);
            distances.push_back(dist(i, j));
        }
    }
    std::sort(distances.begin(), distances.end());
    float median;
    int size = distances.size();
    median= distances[(size-1)/2];
    for(int i=0;i<nb_sample;++i)
        dist(i, i) = median;


    labels = arange<int>(0, nb_sample);
    float damping_factor = 0.8;
    responsibility = zeros<float>(nb_sample, nb_sample);
    availability = zeros<float>(nb_sample, nb_sample);
    Matrix<float> E;

    int it=0;
    while(it++<max_iteration || nb_clusters<=0)
    {

        //update responsibility
        Matrix<float> old_R = responsibility;
        for(int i=0;i<nb_sample;++i)
        {
            for(int k=0;k<nb_sample;++k)
            {
                float maximum = -FLT_MAX;
                for(int kp=0;kp<nb_sample;++kp)
                {
                    float temp = availability(i, kp)+dist(i, kp);
                    if(temp>maximum && kp!=k)
                    {
                        maximum = temp;
                    }
                }
               responsibility(i, k) = dist(i, k) - maximum;
            }
        }
        responsibility = responsibility*(1-damping_factor)+damping_factor*old_R;


        //update availability
        Matrix<float> old_A = availability;
        for(int i=0;i<nb_sample;++i)
        {
            for(int k=0;k<nb_sample;++k)
            {
                if(i!=k)
                {
                    float sum = 0.0f;
                    for(int ip=0;ip<nb_sample;++ip)
                    {
                        if(ip!=k && ip!=i)
                        {
                            sum+=std::max(responsibility(ip, k), 0.0f);
                        }
                    }
                    availability(i, k) = std::min(0.0f, responsibility(k, k)+sum);
                }
                else
                {
                    float sum = 0.0f;
                    for(int ip=0;ip<nb_sample;++ip)
                    {
                        if(ip!=k)
                        {
                            sum+=std::max(responsibility(ip, k), 0.0f);
                        }
                    }
                    availability(k, k) = sum;
                }
            }
        }
        availability = availability*(1-damping_factor)+damping_factor*old_A;
        E = diag(availability)+diag(responsibility);
        nb_clusters = count_nonzero(E>0);
    }

    Matrix<int> cluster_idx = argwhere(E>0);
    centers = Matrix<float>(nb_clusters, nb_features);
    for(int i=0;i<nb_clusters;++i)
        centers.setRow(i, M.getRow(cluster_idx(i, 0)));

}

Matrix<int> AffinityPropagation::fit_predict(const Matrix<float>& M)
{
    fit(M);
    return predict(M);
}

Matrix<int> AffinityPropagation::predict(const Matrix<float>& M)
{
    int nb_samples = M.rowNb();
    labels = Matrix<int>(nb_samples, 1);
    for(int i=0;i<M.rowNb();++i)
    {
        Matrix<float> distances = dot(ones<float>(nb_clusters, 1), M.getRow(i))-centers;
        distances *= distances;
        distances = axissum(distances, 1);
        labels(i, 0) = argmin(distances)(0, 0);
    }
    return labels;
}

