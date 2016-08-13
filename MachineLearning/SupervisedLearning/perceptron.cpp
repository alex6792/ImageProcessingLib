#include "../../linalg.hpp"
#include "perceptron.hpp"


Perceptron::Perceptron()
{
    alpha = 0.1;
    thresh = 0;
    max_iteration = 20;
}

Matrix<float> Perceptron::get_coeff()
{
    return coeff;
}

float Perceptron::get_thresh()
{
    return thresh;
}

void Perceptron::set_alpha(float a)
{
    alpha = a;
}

void Perceptron::set_max_iteration(std::size_t new_max_iteration)
{
    max_iteration = new_max_iteration;
}

void Perceptron::fit(const Matrix<float>& M, const Matrix<bool>& label)
{
    coeff = zeros<float>(M.colNb(), 1);
    std::size_t cpt = 0;
    bool coeff_changed = true;
    while(coeff_changed && cpt<max_iteration)
    {
        coeff_changed = false;
        for(std::size_t i=0;i<M.rowNb();++i)
        {
            Matrix<float> X = M.getRow(i);
            bool res = predict(X)(0, 0);
            if(res!=label(i, 0))
            {
                X.transpose();
                coeff_changed = true;
                if(label(i, 0))
                {
                    coeff+=alpha*X;
                    thresh-=alpha;
                }
                else
                {
                    coeff-=alpha*X;
                    thresh+=alpha;
                }
            }
        }
        ++cpt;
    }
}

Matrix<bool> Perceptron::fit_predict(const Matrix<float>& M, const Matrix<bool>& label)
{
    fit(M, label);
    return predict(M);
}

Matrix<bool> Perceptron::predict(const Matrix<float>& M)
{
    return dot(M, coeff)>thresh;
}
