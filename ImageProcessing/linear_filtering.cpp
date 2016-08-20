#include "linear_filtering.hpp"
#include "../statistics.hpp"


Matrix<float> average(std::size_t filtersize)
{
    return ones<float>(filtersize)/float(filtersize*filtersize);
}

Matrix<float> disk(std::size_t radius)
{
    Matrix<float> new_disk = zeros<float>(2*radius+1);
    for(std::size_t i=0;i<new_disk.rowNb();++i)
    {
        for(std::size_t j=0;j<new_disk.colNb();++j)
        {
            if((i-radius)*(i-radius)+(j-radius)*(j-radius)<=radius*radius)
                new_disk(i, j) = 1.0f;
        }
    }
    return new_disk/sum(new_disk);
}

Matrix<float> gaussian(std::size_t filtersize, float stddev)
{
    std::pair<Matrix<float>, Matrix<float> > XY = meshgrid<float>(filtersize);
    Matrix<float>& X = XY.first;
    Matrix<float>& Y = XY.second;
    X-=filtersize/2.0f-0.5f;
    Y-=filtersize/2.0f-0.5f;
    Matrix<float> F = exp(-(X*X+Y*Y)/(2*stddev*stddev));
    return F/sum(F);
}

Matrix<float> isotropicx()
{
    Matrix<float> iso({{-1.0f,0.0f,1.0f},
                    {-std::sqrt(2.0f),0.0f,sqrt(2.0f)},
                    {-1.0f,0.0f,1.0f}});
    return iso/(2.0f+sqrt(2.0f));
}

Matrix<float> isotropicy()
{
    Matrix<float> iso({{-1.0f,-sqrt(2.0f),-1.0f},
                    {0.0f,0.0f,0.0f},
                    {1.0f,sqrt(2.0f),1.0f}});
    return iso/(2.0f+sqrt(2.0f));
}

Matrix<float> kirch()
{
    Matrix<float> kir({{-3.0f,-3.0f,-3.0f},
                    {-3.0f,0.0f,-3.0f},
                    {5.0f,5.0f,5.0f}});
    return kir/15.0f;
}

Matrix<float> laplacian(float alpha)
{
    Matrix<float> lap({{alpha,1.0f-alpha,alpha},
                    {1.0f-alpha,-4.0f,1.0f-alpha},
                    {alpha,1.0f-alpha,alpha}});
    return lap/4.0f;
}

Matrix<float> log(std::size_t filtersize, float stddev)
{
    std::pair<Matrix<float>, Matrix<float> > XY = meshgrid<float>(filtersize);
    Matrix<float>& X = XY.first;
    Matrix<float>& Y = XY.second;
    X-=filtersize/2.0f-0.5f;
    Y-=filtersize/2.0f-0.5f;
    X*=X;
    Y*=Y;
    float var = stddev*stddev;
    Matrix<float> F = (X+Y-2*var)*exp(-(X+Y)/(2*var));
    return F;
}

Matrix<float> prewittx()
{
    Matrix<float> pre({{-1.0f,0.0f,1.0f},
                    {-1.0f,0.0f,1.0f},
                    {-1.0f,0.0f,1.0f}});
    return pre/3.0f;
}

Matrix<float> prewitty()
{
    Matrix<float> pre({{-1.0f,-1.0f,-1.0f},
                {0.0f,0.0f,0.0f},
                {1.0f,1.0f,1.0f}});
    return pre/3.0f;
}

Matrix<float> robertsx()
{
    Matrix<float> rob = zeros<float>(3);
    rob(1, 2) = 1.0f;
    rob(2, 1) = -1.0f;
    return rob;
}

Matrix<float> robertsy()
{
    Matrix<float> rob = zeros<float>(3);
    rob(1, 1) = 1.0f;
    rob(2, 2) = -1.0f;
    return rob;
}

Matrix<float> robinson()
{
    Matrix<float> rob({{-1.0f,-1.0f,-1.0f},
                    {1.0f,-2.0f,1.0f},
                    {1.0f,1.0f,1.0f}});
    return rob/5.0f;
}


Matrix<float> scharrx()
{
    Matrix<float> sch({{-3.0f,0.0f,3.0f},
                    {-10.0f,0.0f,10.0f},
                    {-3.0f,0.0f,3.0f}});
    return sch/16.0f;
}

Matrix<float> scharry()
{
    Matrix<float> sch({{-3.0f,-10.0f,-3.0f},
                    {0.0f,0.0f,0.0f},
                    {3.0f,10.0f,3.0f}});
    return sch/16.0f;
}

Matrix<float> sobelx()
{
    Matrix<float> sob({{-1.0f,0.0f,1.0f},
                    {-2.0f,0.0f,2.0f},
                    {-1.0f,0.0f,1.0f}});
    return sob/4.0f;
}

Matrix<float> sobely()
{
    Matrix<float> sob({{-1.0f,-2.0f,-1.0f},
                    {0.0f,0.0f,0.0f},
                    {1.0f,2.0f,1.0f}});
    return sob/4.0f;
}

Matrix<float> unsharp(float alpha)
{
    Matrix<float> unsharp({{alpha,1.0f-alpha,alpha},
                    {1.0f-alpha,-4.0f,1.0f-alpha},
                    {alpha,1.0f-alpha,alpha}});
    return unsharp/4.0f;
}

Matrix<unsigned char> gradient(Matrix<unsigned char> M, std::string method)
{
    Matrix<float> M_copy = Matrix<float>(M);
    Matrix<float> fx, fy;
    if(!method.compare("sobel"))
    {
        fx = sobelx();
        fy = sobely();
    }
    else if(!method.compare("scharr"))
    {
        fx = scharrx();
        fy = scharry();
    }
    else if(!method.compare("isotropic"))
    {
        fx = isotropicx();
        fy = isotropicy();
    }

    Matrix<float> Gx = conv(M_copy, fx);
    Matrix<float> Gy = conv(M_copy, fy);
    return Matrix<unsigned char>(sqrt((Gx*Gx+Gy*Gy)/2.0f));
}

Matrix<unsigned char> filter(Matrix<unsigned char> M, Matrix<float> f, int mode)
{
    Matrix<unsigned char> filtered_img(M.rowNb(), M.colNb());
    for(std::size_t i=0, I=M.rowNb();i<I;++i)
    {
        for(std::size_t j=0, J=M.colNb();j<J;++j)
        {
            float weighted_sum = 0;
            for(std::size_t k=0, K=f.rowNb();k<K;++k)
            {
                for(std::size_t l=0, L=f.colNb();l<L;++l)
                {
                    std::size_t x = i+k>K/2?i+k-K/2:0;
                    std::size_t y = j+l>L/2?j+l-L/2:0;

                    if(x>=I)
                        x = I-1;
                    if(y>=J)
                        y = J-1;
                    weighted_sum+=f(k, l)*M(x, y);
                }
            }
            filtered_img(i, j) = round(weighted_sum);
        }
    }
    return filtered_img;
}

Matrix<float> conv(Matrix<float> M, Matrix<float> f, int mode)
{
    std::size_t H = f.rowNb(), W = f.colNb();
    Matrix<float> filtered_img = zeros<float>(M.rowNb(), M.colNb());
    for(std::size_t i=0;i<H;++i)
    {
        for(std::size_t j=0;j<W;++j)
        {
            Matrix<float> temp = f(i, j)*M;
            for(std::size_t k=H/2;k<i;++k)
            {
                temp.newRow();
                temp.setRow(temp.rowNb()-1, temp.getRow(temp.rowNb()-2));
                temp.delRow(0);
            }
            for(std::size_t k=i;k<H/2;++k)
            {
                temp.newRow(0);
                temp.setRow(0, temp.getRow(1));
                temp.delRow(temp.rowNb()-1);
            }
            for(std::size_t k=W/2;k<j;++k)
            {
                temp.newCol();
                temp.setCol(temp.colNb()-1, temp.getCol(temp.colNb()-2));
                temp.delCol(0);
            }
            for(std::size_t k=j;k<W/2;++k)
            {
                temp.newCol(0);
                temp.setCol(0, temp.getCol(1));
                temp.delCol(temp.colNb()-1);
            }
            filtered_img+=temp;
        }
    }
    return filtered_img;
}
/*
binomial

1 4 6 4 1
4 16 24 16 4
6 24 36 24 6
4 16 24 16 4
1 4 6 4 1

pyramidal
1 2 3 2 1
2 4 6 4 2
3 6 9 6 3
2 4 6 4 2
1 2 3 2 1
*/
