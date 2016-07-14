#include <fstream>
#include <sstream>
#include <string>
#include "../../linalg.hpp"
#include "mlp.hpp"


MLP::MLP(const std::vector<int>& hidden_layers_sizes_arg)
{
    alpha = 0.003;
    max_iteration = 20000;
    nb_clusters = -1;
    nb_features = -1;

    hidden_layers_sizes = hidden_layers_sizes_arg;
    nb_hidden_layers = hidden_layers_sizes.size()+1;
    //activation_fct = [](float value){return 1.0f/(1.0f+std::exp(-value));};
    //activation_fct_der = [](float value){float temp = 1.0f/(1.0f+std::exp(-value)); return temp*(1.0f-temp);};

    activation_fct = std::tanh;
    activation_fct_der = [](float value){float temp = std::tanh(value); return 1.0f-temp*temp;};
}

void MLP::export_graphviz(std::string filename)
{
    std::ofstream myfile(filename, std::ios::out | std::ios::trunc);
    if(myfile)
    {
        myfile<<"digraph G {"<<std::endl;
        hidden_layers_sizes.insert(hidden_layers_sizes.begin(), nb_features);
        for(int i=0;i<nb_hidden_layers+1;++i)
        {
            myfile<<'\t'<<"subgraph cluster_"<<i<<" {"<<std::endl;
            for(int j=0;j<hidden_layers_sizes[i];++j)
                myfile<<'\t'<<'\t'<<"L"<<i<<"_"<<j<<";"<<std::endl;
            if(i!=nb_hidden_layers)
                myfile<<'\t'<<'\t'<<"L"<<i<<"_b"<<";"<<std::endl;
            myfile<<'\t'<<"}"<<std::endl;
        }
        for(int i=0;i<nb_hidden_layers;++i)
        {
            for(int k=0;k<hidden_layers_sizes[i+1];++k)
            {
                for(int j=0;j<hidden_layers_sizes[i];++j)
                    myfile<<'\t'<<"L"<<i<<"_"<<j<<"->"<<"L"<<i+1<<"_"<<k<<" [label = \""<<"A"<<j<<","<<k<<" : "<<A[i](j, k)<<"\" ];"<<std::endl;
                myfile<<'\t'<<"L"<<i<<"_b"<<"->"<<"L"<<i+1<<"_"<<k<<" [label = \""<<"b"<<k<<" : "<<biases[i](0, k)<<"\" ];"<<std::endl;
            }
        }
        hidden_layers_sizes.erase(hidden_layers_sizes.begin());
        myfile<<"}"<<std::endl;
        myfile.close();
    }
    else
        std::cout<<"unable to write in the file : "<<filename<<std::endl;
}

void MLP::set_alpha(float a)
{
    alpha = a;
}

void MLP::set_max_iteration(int new_max_iteration)
{
    max_iteration = new_max_iteration;
}

void MLP::fit(const Matrix<float>& M, const Matrix<int>& label)
{
    int nb_samples = M.rowNb();
    nb_features = M.colNb();
    nb_clusters = max(label)+1;
    hidden_layers_sizes.push_back(nb_clusters);
    nb_hidden_layers = hidden_layers_sizes.size();
    Matrix<float> true_labels = zeros<float>(nb_samples, nb_clusters);
    for(int i=0;i<nb_samples;++i)
        true_labels(i, label(i, 0)) = 1.0f;
    y = std::vector<Matrix<float> >(nb_hidden_layers);
    z = std::vector<Matrix<float> >(nb_hidden_layers+1);
    A = std::vector<Matrix<float> >(nb_hidden_layers);
    biases = std::vector<Matrix<float> >(nb_hidden_layers);
    A[0] = Matrix<float>(rand<int>(nb_features, hidden_layers_sizes[0])%10)/10.0f;
    for(int i=1;i<nb_hidden_layers;++i)
        A[i] = Matrix<float>(rand<int>(hidden_layers_sizes[i-1], hidden_layers_sizes[i])%10)/10.0f;
    for(int i=0;i<nb_hidden_layers;++i)
        biases[i] = zeros<float>(1, hidden_layers_sizes[i]);

    //A[0] = {{-1, 1},{-1, 1}};
    //biases[0] = {{0.5, -0.5}};
    for(int it=0;it<max_iteration;++it)
    {
        for(int i=0;i<nb_samples;++i)
        {

            Matrix<float> cur_label = true_labels.getRow(i);
            Matrix<float> cur_sample = M.getRow(i);
            Matrix<int> cur_prediction = predict(cur_sample);

            std::vector<Matrix<float> > y_der(nb_hidden_layers);
            std::vector<Matrix<float> > z_der(nb_hidden_layers+1);

            z_der[nb_hidden_layers] = z[nb_hidden_layers]-cur_label;
            y_der[nb_hidden_layers-1] = apply<float, float>(y[nb_hidden_layers-1], activation_fct_der)*z_der[nb_hidden_layers];
            for(int j=nb_hidden_layers-2;j>=0;--j)
            {
                z_der[j+1] = dot(y_der[j+1], transpose(A[j+1]));
                y_der[j] = apply<float, float>(y[j], activation_fct_der)*z_der[j+1];
            }
            for(int j=0;j<nb_hidden_layers;++j)
            {
                A[j]-=alpha*dot(transpose(z[j]), y_der[j]);
                biases[j]-=alpha*y_der[j];
            }
        }
    }
}

Matrix<int> MLP::fit_predict(const Matrix<float>& M, const Matrix<int>& label)
{
    fit(M, label);
    return predict(M);
}

Matrix<int> MLP::predict(const Matrix<float>& M)
{
    int nb_samples = M.rowNb();
    Matrix<int> labels(nb_samples, 1);
    for(int i=0;i<nb_samples;++i)
    {
        z[0] = M.getRow(i);
        for(int j=0;j<nb_hidden_layers;++j)
        {
            y[j] = dot(z[j], A[j])+biases[j];
            z[j+1] = apply<float, float>(y[j], activation_fct);
        }
        labels(i, 0) = argmax(z[nb_hidden_layers])(0, 1);
    }
    return labels;
}
