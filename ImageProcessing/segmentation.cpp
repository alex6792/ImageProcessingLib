#include "segmentation.hpp"
#include <algorithm>
#include <queue>
#include <vector>

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
