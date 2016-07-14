#include "meanshift.hpp"
#include "../../linalg.hpp"
#include "../../mmath.hpp"


MeanShift::MeanShift()
{
    //kernel = [](float x){return std::exp(-0.5f*x*x/1.0f);};
    kernel = [](float x){return float(x<3.0f);};
    max_iteration = 20;
}

void MeanShift::fit(const Matrix<float>& M)
{
    data = M;
    nb_features = M.colNb();
    int nb_samples = M.rowNb();

    Matrix<float> M_copy = M;
    int it = 0;
    M_copy = iterate(M_copy);
    while(it<max_iteration && max(pow(M_copy-data,2.0f)>1e-9))
    {
        data = M_copy;
        M_copy = iterate(M_copy);
        ++it;
    }

    labels = arange(0, nb_samples);
    for(int i=1;i<nb_samples;++i)
    {
        Matrix<float> A = M_copy.getRow(i);
        for(int j=0;j<i;++j)
        {
            Matrix<float> B = M_copy.getRow(j);
            if(max(pow(A-B, 2.0f))<1e-9)
            {
                labels(i, 0) = labels(j, 0);
                break;
            }
        }
    }
    Matrix<int> unique_labels = unique(labels);
    nb_clusters = unique_labels.size();
    for(int i=0;i<nb_clusters;++i)
        replace(labels, unique_labels(i, 0), i);

    centers = Matrix<float>(nb_clusters, nb_features);
    for(int i=0;i<nb_clusters;++i)
    {
        Matrix<int> indices = argwhere(labels==i);
        centers.setRow(i, M_copy.getRow(indices(0, 0)));
    }
    std::cout<<centers<<std::endl;
}

Matrix<int> MeanShift::fit_predict(const Matrix<float>& M)
{
    fit(M);
    return labels;
}

Matrix<int> MeanShift::predict(const Matrix<float>& M)
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

Matrix<float> MeanShift::iterate(const Matrix<float>& M)
{
    int nb_samples = M.rowNb();
    nb_features = M.colNb();
    Matrix<float> M_copy = zeros<float>(nb_samples, nb_features);

    for(int i=0;i<nb_samples;++i)
    {
        float sum_w = 0.0f;
        Matrix<float> input_row = M.getRow(i);
        Matrix<float> output_row = zeros<float>(1, nb_features);
        for(int j=0, J=data.rowNb();j<J;++j)
        {
            Matrix<float> cur_sample = data.getRow(j);
            float w = kernel(sqrt(sum(pow(cur_sample-input_row, 2.0f))));
            sum_w+=w;
            output_row+=w*cur_sample;
        }
        M_copy.setRow(i, output_row/sum_w);
    }
    return M_copy;
}
