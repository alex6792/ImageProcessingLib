#include "scaler.hpp"
#include "../linalg.hpp"
#include "../statistics.hpp"


MinMaxScaler::MinMaxScaler(float range_min_arg, float range_max_arg, int axis_arg)
{
    range_min = range_min_arg;
    range_max = range_max_arg;
    axis = axis_arg;
}

void MinMaxScaler::fit(const Matrix<float>& M)
{
    data_min = axismin(M, axis);
    data_max = axismax(M, axis);
    data_range = data_max-data_min;
    data_range_inv = 1.0f/data_range;
}

Matrix<float> MinMaxScaler::fit_transform(const Matrix<float>& M)
{
    fit(M);
    return transform(M);
}

Matrix<float> MinMaxScaler::inverse_transform(const Matrix<float>& M)
{
    if(axis==1)
        return (M-range_min)*dot(data_range, ones<float>(1, M.colNb()))/(range_max-range_min)+dot(data_min, ones<float>(1, M.colNb()));
    else if(axis==2)
        return (M-range_min)*dot(ones<float>(M.rowNb(), 1), data_range)/(range_max-range_min)+dot(ones<float>(M.rowNb(), 1), data_min);
    else
        return (M-range_min)*data_range(0, 0)/(range_max-range_min)+data_min(0, 0);
}

Matrix<float> MinMaxScaler::transform(const Matrix<float>& M)
{
    if(axis==1)
        return (range_max-range_min)*(M-dot(data_min, ones<float>(1, M.colNb())))*dot(data_range_inv, ones<float>(1, M.colNb()))+range_min;
    else if(axis==2)
        return (range_max-range_min)*(M-dot(ones<float>(M.rowNb(), 1), data_min))*dot(ones<float>(M.rowNb(), 1), data_range_inv)+range_min;
    else
        return (range_max-range_min)*(M-data_min(0, 0))*data_range_inv(0, 0)+range_min;
}



StandardScaler::StandardScaler(bool with_mean_arg, bool with_std_arg, int axis_arg)
{
    with_mean = with_mean_arg;
    with_std = with_std_arg;
    axis = axis_arg;
}

void StandardScaler::fit(const Matrix<float>& M)
{
    data_mean = axismean(M, axis);
    data_std = axisstdev(M, axis);
    data_std_inv = 1.0f/data_std;
}

Matrix<float> StandardScaler::fit_transform(const Matrix<float>& M)
{
    fit(M);
    return transform(M);
}

Matrix<float> StandardScaler::inverse_transform(const Matrix<float>& M)
{
    Matrix<float> result = M;
    if(with_std)
        result*=data_std;
    if(with_mean)
        result+=data_mean;
    return result;
}

Matrix<float> StandardScaler::transform(const Matrix<float>& M)
{
    Matrix<float> result = M;
    if(axis==1)
    {
        if(with_mean)
            result-=dot(data_mean, ones<float>(1, M.colNb()));
        if(with_std)
            result*=dot(data_std_inv, ones<float>(1, M.colNb()));
    }
    else if(axis==2)
    {
        if(with_mean)
            result-=dot(ones<float>(M.rowNb(), 1), data_mean);
        if(with_std)
            result*=dot(ones<float>(M.rowNb(), 1), data_std_inv);
    }
    else
    {
        if(with_mean)
            result-=data_mean(0, 0);
        if(with_std)
            result*=data_std_inv(0, 0);
    }

    return result;
}


