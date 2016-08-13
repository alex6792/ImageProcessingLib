#include <fstream>
#include <sstream>
#include <string>
#include "../../linalg.hpp"
#include "../../ImageProcessing/imgconverter.hpp"
#include "dtree.hpp"


DtreeNode::DtreeNode()
{
    feature = 0;
    thresh = 0.0f;
    gini = 1.0f;
    samples = 0;
    label = 0;
    value = Matrix<float>();
    left_child = 0;
    right_child = 0;
    is_a_leaf = true;
}

void DtreeNode::show()
{
    std::cout<<feature<<" ";
    std::cout<<thresh<<" ";
    std::cout<<gini<<" ";
    std::cout<<samples<<" ";
    std::cout<<label<<" ";
    //std::cout<<value<<" ";
    std::cout<<left_child<<" ";
    std::cout<<right_child<<std::endl;
}


Dtree::Dtree()
{
    nb_clusters = 0;
    nb_features = 0;
    tree = std::vector<DtreeNode>(1);
}

void Dtree::export_graphviz(std::string filename)
{
    std::ofstream myfile(filename, std::ios::out | std::ios::trunc);
    if(myfile)
    {
        myfile<<"digraph Tree {"<<std::endl;
        myfile<<"node [shape=box]"<<std::endl;
        std::size_t n = tree.size();
        for(std::size_t i=0;i<n;++i)
        {
            myfile<<i;
            myfile<<" [label=\"";
            if(!tree[i].is_a_leaf)
            {
                myfile<<"X["<<tree[i].feature<<"]";
                myfile<<" <= "<<tree[i].thresh<<"\\n";
            }
            myfile<<"gini = "<<tree[i].gini<<"\\n";
            myfile<<"samples = "<<tree[i].samples<<"\\n";
            //myfile<<"value = "<<tree[i].value;
            myfile<<"\"] ;";
            myfile<<std::endl;

            if(!tree[i].is_a_leaf)
            {
                myfile<<i<<" -> "<<tree[i].left_child;
                //if(i==0)
                    //myfile<<" [labeldistance=2.5, labelangle=-45, headlabel=\"False\"]";
                myfile<<" ;"<<std::endl;

                myfile<<i<<" -> "<<tree[i].right_child;
                //if(i==0)
                    //myfile<<" [labeldistance=2.5, labelangle=45, headlabel=\"True\"]";
                myfile<<" ;"<<std::endl;
            }
        }
        myfile<<"}"<<std::endl;
        myfile.close();
    }
    else
        std::cout<<"unable to write in the file : "<<filename<<std::endl;
}


void Dtree::fit(const Matrix<float>& M, const Matrix<std::size_t>& label)
{
    nb_features = M.colNb();
    nb_clusters = max(label)+1;

    split_node(tree[0], M, label);
}

Matrix<std::size_t> Dtree::fit_predict(const Matrix<float>& M, const Matrix<std::size_t>& label)
{
    fit(M, label);
    return predict(M);
}

Matrix<std::size_t> Dtree::predict(const Matrix<float>& M)
{
    std::size_t nb_samples = M.rowNb();
    Matrix<std::size_t> labels(nb_samples, 1);
    for(std::size_t i=0;i<nb_samples;++i)
    {
        const DtreeNode* cur_node = tree.data();
        while(!cur_node->is_a_leaf)
        {
            if(M(i, cur_node->feature) > cur_node->thresh)
                cur_node = &tree[cur_node->left_child];
            else
                cur_node = &tree[cur_node->right_child];
        }
        labels(i, 0) = cur_node->label;
    }
    return labels;
}

std::pair<std::size_t, float> Dtree::find_best_split(const Matrix<float>& M, const Matrix<std::size_t>& label)
{
    Matrix<std::size_t> unique_labels = unique(label);
    Matrix<float> priors(unique_labels.rowNb(), 1);
    for(std::size_t i=0;i<unique_labels.rowNb();++i)
        priors(i, 0) = count(label, unique_labels(i, 0));
    priors/=sum(priors);
    float gini = prod(priors);
    std::size_t best_feature = 0;
    float best_thresh = 0.0f;
    float best_gini = 0.0f;
    for(std::size_t i=0;i<nb_features;++i)
    {
        Matrix<float> cur_feature = M.getCol(i);
        Matrix<float> sorted_feature = unique(cur_feature);
        Matrix<float> threshes(sorted_feature.rowNb(), 1);
        std::adjacent_difference(sorted_feature.cbegin(), sorted_feature.cend(), threshes.begin(), std::plus<float>());
        threshes/=2.0f;
        for(std::size_t j=1;j<threshes.rowNb();++j)
        {
            Matrix<bool> A = cur_feature>threshes(j, 0);
            Matrix<bool> NOTA = NOT(A);
            Matrix<float> cur_hist = Matrix<float>(histogram(A));
            Matrix<float> proba_0(unique_labels.rowNb(), 1);
            Matrix<float> proba_1(unique_labels.rowNb(), 1);
            for(std::size_t k=0;k<unique_labels.rowNb();++k)
            {
                proba_0(k, 0) = count_nonzero(NOTA*(label==unique_labels(k, 0)));
                proba_1(k, 0) = count_nonzero(A*(label==unique_labels(k, 0)));
            }
            cur_hist/=sum(cur_hist);
            proba_0/=sum(proba_0);
            proba_1/=sum(proba_1);
            float cur_gini = gini-cur_hist(0, 0)*prod(proba_0)-cur_hist(0, 1)*prod(proba_1);
            if(cur_gini>best_gini)
            {
                best_feature = i;
                best_gini = cur_gini;
                best_thresh = threshes(j, 0);
            }
        }
    }
    return std::make_pair(best_feature, best_thresh);
}

void Dtree::split_node(DtreeNode& node, const Matrix<float>& M, const Matrix<std::size_t>& label)
{
    Matrix<std::size_t> unique_label = unique(label);
    if(unique_label.size()>1)
    {
        std::pair<std::size_t, float> split = find_best_split(M, label);
        std::size_t n = tree.size();
        node.feature = split.first;
        node.thresh = split.second;
        node.samples = M.rowNb();
        node.left_child = n;
        node.right_child = n+1;
        node.is_a_leaf = false;
        Matrix<float> H = Matrix<float>(histogram(label));
        node.value = H;
        H/=node.samples;
        node.gini = sum(H*(1.0f-H));
        Matrix<float> cur_feature = M.getCol(node.feature);
        Matrix<bool> cdt = cur_feature > node.thresh;
        Matrix<bool> not_cdt = NOT(cdt);
        Matrix<float> A = compress(cdt, M, 1);
        Matrix<float> B = compress(not_cdt, M, 1);
        Matrix<std::size_t> C = compress(cdt, label, 1);
        Matrix<std::size_t> D = compress(not_cdt, label, 1);
        tree.push_back(DtreeNode());
        tree.push_back(DtreeNode());
        split_node(tree[n], A, C);
        split_node(tree[n+1], B, D);
    }
    else
    {
        node.gini = 0.0f;
        node.samples = M.rowNb();
        node.value = Matrix<float>(histogram(label));
        node.label = unique_label(0, 0);
    }
}
