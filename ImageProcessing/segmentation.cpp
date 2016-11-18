#include "../mmath.hpp"
#include "../statistics.hpp"
#include "imgconverter.hpp"
#include "morphology.hpp"
#include "segmentation.hpp"

#include <algorithm>
#include <cfloat>
#include <queue>
#include <vector>


Edge::Edge()
{
}

Edge::Edge(std::size_t v1_arg, std::size_t v2_arg, unsigned char w_arg)
{
    w = w_arg;
    v1 = v1_arg;
    v2 = v2_arg;
}

bool Edge::operator<(const Edge& e) const
{
    return w<e.w;
}


Matrix<bool> hysteresis(const Matrix<unsigned char>& img, unsigned char low_thresh, unsigned char high_thresh)
{
    return reconstruct(img>=low_thresh, img>=high_thresh);
}

Matrix<bool> otsu(const Matrix<unsigned char>& img)
{
    Matrix<std::size_t> hist = histogram(img);

    // BG
    Matrix<float> hist_normalized = Matrix<float>(hist)/float(img.size());
    Matrix<float> cum_hist = cumsum(hist_normalized);
    auto XY = meshgrid<float>(hist.rowNb(), hist.colNb());
    Matrix<float>& Y = XY.second;
    Matrix<float> mean_BG = hist_normalized*Y;
    Matrix<float> var_BG = mean_BG*Y;
    mean_BG = cumsum(mean_BG);
    var_BG = cumsum(var_BG);
    std::transform(mean_BG.cbegin(), mean_BG.cend(), cum_hist.cbegin(), mean_BG.begin(), std::divides<float>());
    std::transform(var_BG.cbegin(), var_BG.cend(), cum_hist.cbegin(), var_BG.begin(), std::divides<float>());
    var_BG-=mean_BG*mean_BG;

    // FG
    hist_normalized.fliplr();
    Matrix<float> cum_hist_inv = cumsum(hist_normalized);
    Y.fliplr();
    Matrix<float> mean_FG = hist_normalized*Y;
    Matrix<float> var_FG = mean_FG*Y;
    mean_FG = cumsum(mean_FG);
    var_FG = cumsum(var_FG);
    std::transform(mean_FG.cbegin(), mean_FG.cend(), cum_hist_inv.cbegin(), mean_FG.begin(), std::divides<float>());
    std::transform(var_FG.cbegin(), var_FG.cend(), cum_hist_inv.cbegin(), var_FG.begin(), std::divides<float>());
    var_FG-=mean_FG*mean_FG;
    cum_hist_inv.fliplr();
    mean_FG.fliplr();
    var_FG.fliplr();
    cum_hist.delCol(cum_hist.colNb()-1);
    var_BG.delCol(var_BG.colNb()-1);
    cum_hist_inv.delCol(0);
    var_FG.delCol(0);

    Matrix<float> weigthed_sum = cum_hist*var_BG+cum_hist_inv*var_FG;
    std::replace_if(weigthed_sum.begin(), weigthed_sum.end(), [](float value){return std::isnan(value);}, FLT_MAX);
    Matrix<std::size_t> results = argmin(weigthed_sum);
    return img>results(0, 1);
}

Matrix<std::size_t> felzenszwalb(const Matrix<unsigned char>& M)
{
    std::size_t H = M.rowNb(), W = M.colNb();
    Matrix<float> M_copy_l(M);
    Matrix<float> M_copy_r(M);
    Matrix<float> M_copy_u(M);
    Matrix<float> M_copy_d(M);
    M_copy_l.delCol(0);
    M_copy_r.delCol(W-1);
    M_copy_u.delRow(0);
    M_copy_d.delRow(H-1);
    Matrix<unsigned char> diff_img_h(abs(M_copy_l-M_copy_r));
    Matrix<unsigned char> diff_img_v(abs(M_copy_u-M_copy_d));
    Matrix<std::size_t> labeled_img = arange<std::size_t>(0, H*W);
    Matrix<std::vector<std::size_t> >cluster_elements(1, H*W);
    auto it_m = cluster_elements.begin();
    for(std::size_t i=0;i<H*W;++i)
    {
        (*it_m).clear();
        (*it_m).push_back(i);
        ++it_m;
    }
    Matrix<unsigned char> cluster_Int = zeros<unsigned char>(1, H*W);

    std::vector<Edge> edges(W*(H-1)+H*(W-1));
    auto it = diff_img_v.cbegin();
    for(std::size_t i=0;i<W*(H-1);++i)
    {
        std::size_t x = i/W;
        std::size_t y = i%W;
        edges[i] = Edge(x*W+y,(x+1)*W+y,*it++);
    }
    it = diff_img_h.cbegin();
    for(std::size_t i=0;i<(W-1)*H;++i)
    {
        std::size_t x = i/(W-1);
        std::size_t y = i%(W-1);
        edges[i+W*(H-1)] = Edge(x*W+y,x*W+y+1,*it++);
    }
    std::sort(edges.begin(), edges.end(), std::less<Edge>());
    std::size_t K = 10000;
    for(std::size_t i=0;i<W*(H-1)+H*(W-1);++i)
    {
        //std::cout<<i<<"/"<<W*(H-1)+H*(W-1)<<std::endl;
        std::size_t l1 = labeled_img(edges[i].v1, 0);
        std::size_t l2 = labeled_img(edges[i].v2, 0);

        if(l1!=l2)
        {
            std::vector<std::size_t>& C1 = cluster_elements(0, l1);
            std::vector<std::size_t>& C2 = cluster_elements(0, l2);
            std::size_t n1 = C1.size();
            std::size_t n2 = C2.size();
            unsigned char& w = edges[i].w;
            unsigned char& w1 = cluster_Int(0, l1);
            unsigned char& w2 = cluster_Int(0, l2);

            if((w-w1)*n1<=K && (w-w2)*n2<=K)
            {
                std::for_each(C2.cbegin(), C2.cend(), [&labeled_img, l1](std::size_t idx){labeled_img(idx, 0)=l1;});
                w1 = w;
                C1.insert(C1.end(), C2.cbegin(), C2.cend());
                C2.clear();
            }
        }
    }

    labeled_img.reshape(H, W);
    return labeled_img;
}

Matrix<std::size_t> watershed(const Matrix<unsigned char>& M, const Matrix<std::size_t>& markers)
{
    std::size_t H = M.rowNb(), W = M.colNb();
    Matrix<std::size_t> labeled_img = markers;
    Matrix<std::size_t> seeds = argwhere(markers>0);
    std::vector<std::queue<Matrix<std::size_t> > > Hqueues(256);
    for(std::size_t i=0, I=seeds.rowNb();i<I;++i)
    {
        Hqueues[M(seeds(i, 0), seeds(i, 1))].push(seeds.getRow(i));
    }

    while(!std::all_of(Hqueues.cbegin(), Hqueues.cend(), [](const std::queue<Matrix<std::size_t> >& q){return q.empty();}))
    {
        unsigned char c = 0;
        while(Hqueues[c].empty())
        {
            ++c;
        }
        Matrix<std::size_t> x = Hqueues[c].front();
        Hqueues[c].pop();
        std::size_t cur_label = labeled_img(x(0, 0), x(0, 1));
        for(std::size_t i=0;i<3;++i)
        {
            for(std::size_t j=0;j<3;++j)
            {
                if(x(0, 0)+i>0 && x(0,1)+j>0 && x(0, 0)+i<=H && x(0,1)+j<=W && (!(j==1 && i==1)))
                {
                    Matrix<std::size_t> y = {{x(0, 0)+i-1, x(0,1)+j-1}};
                    unsigned char v = M(y(0, 0), y(0, 1));
                    if(labeled_img(y(0, 0), y(0, 1))==0)
                    {
                        labeled_img(y(0, 0), y(0, 1)) = cur_label;
                        Hqueues[v].push(y);
                    }
                }
            }
        }
    }
    return labeled_img;
}
