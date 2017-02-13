#include "linear_filtering.hpp"
#include "../linalg.hpp"
#include "../statistics.hpp"


Matrix<float> average(std::size_t filtersize)
{
    return ones<float>(filtersize)/float(filtersize*filtersize);
}

Matrix<float> binomial(std::size_t filtersize)
{
    Matrix<float> v(1, filtersize);
    for(std::size_t i=0;i<filtersize;++i)
        v(0, i) = std::tgamma(filtersize)/std::tgamma(filtersize-i)/std::tgamma(i+1);
    v = dot(transpose(v), v);
    return v/sum(v);
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
    Matrix<float> iso({{-1.0f,-std::sqrt(2.0f),-1.0f},
                    {0.0f,0.0f,0.0f},
                    {1.0f,std::sqrt(2.0f),1.0f}});
    return iso/(2.0f+std::sqrt(2.0f));
}

Matrix<float> isotropicy()
{
    Matrix<float> iso({{-1.0f,0.0f,1.0f},
                    {-std::sqrt(2.0f),0.0f,sqrt(2.0f)},
                    {-1.0f,0.0f,1.0f}});
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

Matrix<float> mdifx(std::size_t filtersize, float stddev)
{
    std::pair<Matrix<float>, Matrix<float> > XY = meshgrid<float>(filtersize);
    Matrix<float>& X = XY.first;
    Matrix<float>& Y = XY.second;
    X-=filtersize/2.0f-0.5f;
    Y-=filtersize/2.0f-0.5f;
    float variance = stddev*stddev;
    Matrix<float> F = X*exp(-(X*X+Y*Y)/(2*variance));
    float alpha = 2.0f/sum(abs(F));
    return alpha*F;
}

Matrix<float> mdify(std::size_t filtersize, float stddev)
{
    std::pair<Matrix<float>, Matrix<float> > XY = meshgrid<float>(filtersize);
    Matrix<float>& X = XY.first;
    Matrix<float>& Y = XY.second;
    X-=filtersize/2.0f-0.5f;
    Y-=filtersize/2.0f-0.5f;
    float variance = stddev*stddev;
    Matrix<float> F = Y*exp(-(X*X+Y*Y)/(2*variance));
    float alpha = 2.0f/sum(abs(F));
    return alpha*F;
}

Matrix<float> mdifxx(std::size_t filtersize, float stddev)
{
    std::pair<Matrix<float>, Matrix<float> > XY = meshgrid<float>(filtersize);
    Matrix<float>& X = XY.first;
    Matrix<float>& Y = XY.second;
    X-=filtersize/2.0f-0.5f;
    Y-=filtersize/2.0f-0.5f;
    float variance = stddev*stddev;
    Matrix<float> F = (X*X-variance)*exp(-(X*X+Y*Y)/(2*variance));
    float alpha = 2.0f/sum(abs(F));
    return alpha*F;
}

Matrix<float> mdifxy(std::size_t filtersize, float stddev)
{
    std::pair<Matrix<float>, Matrix<float> > XY = meshgrid<float>(filtersize);
    Matrix<float>& X = XY.first;
    Matrix<float>& Y = XY.second;
    X-=filtersize/2.0f-0.5f;
    Y-=filtersize/2.0f-0.5f;
    float variance = stddev*stddev;
    Matrix<float> F = X*Y*exp(-(X*X+Y*Y)/(2*variance));
    float alpha = 2.0f/sum(abs(F));
    return alpha*F;
}

Matrix<float> mdifyy(std::size_t filtersize, float stddev)
{
    std::pair<Matrix<float>, Matrix<float> > XY = meshgrid<float>(filtersize);
    Matrix<float>& X = XY.first;
    Matrix<float>& Y = XY.second;
    X-=filtersize/2.0f-0.5f;
    Y-=filtersize/2.0f-0.5f;
    float variance = stddev*stddev;
    Matrix<float> F = (Y*Y-variance)*exp(-(X*X+Y*Y)/(2*variance));
    float alpha = 2.0f/sum(abs(F));
    return alpha*F;
}

Matrix<float> prewittx()
{
    Matrix<float> pre({{-1.0f,-1.0f,-1.0f},
                    {0.0f,0.0f,0.0f},
                    {1.0f,1.0f,1.0f}});
    return pre/3.0f;
}

Matrix<float> prewitty()
{
    Matrix<float> pre({{-1.0f,0.0f,1.0f},
                    {-1.0f,0.0f,1.0f},
                    {-1.0f,0.0f,1.0f}});
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


Matrix<float> savgol(std::size_t filtersize, std::size_t degree, std::size_t derivx, std::size_t derivy)
{
    Matrix<float> J(filtersize*filtersize, (degree+1)*(degree+2)/2);

    std::pair<Matrix<float>, Matrix<float> > XY = meshgrid<float>(filtersize);
    Matrix<float>& X = XY.first;
    Matrix<float>& Y = XY.second;
    X-=filtersize/2.0f-0.5f;
    Y-=filtersize/2.0f-0.5f;
    X.reshape(filtersize*filtersize,1);
    Y.reshape(filtersize*filtersize,1);

    std::size_t cur_col = 0;
    for(std::size_t i=0;i<=degree;++i)
    {
        for(std::size_t j=0;j<=i;++j)
        {
            J.setCol(cur_col, pow(X, float(i-j))*pow(Y, float(j)));
            ++cur_col;
        }
    }

    Matrix<float> F = pinv(J);
    if(derivx==0 && derivy==0)
        F = F.getRow(0);
    else if(derivx==1 && derivy==0)
        F = F.getRow(1);
    else if(derivx==0 && derivy==1)
        F = F.getRow(2);
    else if(derivx==2 && derivy==0)
        F = F.getRow(3)*2.0f;
    else if(derivx==1 && derivy==1)
        F = F.getRow(4);
    else if(derivx==0 && derivy==2)
        F = F.getRow(5)*2.0f;

    F.reshape(filtersize, filtersize);
    return F;
}

Matrix<float> scharrx()
{
    Matrix<float> sch({{-3.0f,-10.0f,-3.0f},
                    {0.0f,0.0f,0.0f},
                    {3.0f,10.0f,3.0f}});
    return sch/16.0f;
}

Matrix<float> scharry()
{
    Matrix<float> sch({{-3.0f,0.0f,3.0f},
                    {-10.0f,0.0f,10.0f},
                    {-3.0f,0.0f,3.0f}});
    return sch/16.0f;
}

Matrix<float> sobelx()
{
    Matrix<float> sob({{-1.0f,-2.0f,-1.0f},
                    {0.0f,0.0f,0.0f},
                    {1.0f,2.0f,1.0f}});
    return sob/4.0f;
}

Matrix<float> sobely()
{
    Matrix<float> sob({{-1.0f,0.0f,1.0f},
                    {-2.0f,0.0f,2.0f},
                    {-1.0f,0.0f,1.0f}});
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
    else if(!method.compare("mdif"))
    {
        fx = mdifx();
        fy = mdify();
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

Matrix<float> conv(Matrix<float> M, Matrix<float> f, std::string mode)
{
    std::size_t H = f.rowNb(), W = f.colNb();
    Matrix<float> filtered_img = zeros<float>(M.rowNb(), M.colNb());
    for(std::size_t i=0;i<H;++i)
    {
        for(std::size_t j=0;j<W;++j)
        {
            Matrix<float> temp = f(H-i, W-j)*M;
            for(std::size_t k=H/2;k<i;++k)
            {
                if(mode.compare("valid"))
                {
                    temp.newRow();
                    temp.setRow(temp.rowNb()-1, temp.getRow(temp.rowNb()-2));
                }
                if(mode.compare("full"))
                    temp.delRow(0);
            }
            for(std::size_t k=i;k<H/2;++k)
            {
                if(mode.compare("valid"))
                {
                    temp.newRow(0);
                    temp.setRow(0, temp.getRow(1));
                }
                if(mode.compare("full"))
                    temp.delRow(temp.rowNb()-1);
            }
            for(std::size_t k=W/2;k<j;++k)
            {
                if(mode.compare("valid"))
                {
                    temp.newCol();
                    temp.setCol(temp.colNb()-1, temp.getCol(temp.colNb()-2));
                }
                if(mode.compare("full"))
                    temp.delCol(0);
            }
            for(std::size_t k=j;k<W/2;++k)
            {
                if(mode.compare("valid"))
                {
                    temp.newCol(0);
                    temp.setCol(0, temp.getCol(1));
                }
                if(mode.compare("full"))
                    temp.delCol(temp.colNb()-1);
            }
            filtered_img+=temp;
        }
    }
    return filtered_img;
}

Matrix<float> xcorr(Matrix<float> M, Matrix<float> f, std::string mode)
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
                if(mode.compare("valid"))
                {
                    temp.newRow();
                    temp.setRow(temp.rowNb()-1, temp.getRow(temp.rowNb()-2));
                }
                if(mode.compare("full"))
                    temp.delRow(0);
            }
            for(std::size_t k=i;k<H/2;++k)
            {
                if(mode.compare("valid"))
                {
                    temp.newRow(0);
                    temp.setRow(0, temp.getRow(1));
                }
                if(mode.compare("full"))
                    temp.delRow(temp.rowNb()-1);
            }
            for(std::size_t k=W/2;k<j;++k)
            {
                if(mode.compare("valid"))
                {
                    temp.newCol();
                    temp.setCol(temp.colNb()-1, temp.getCol(temp.colNb()-2));
                }
                if(mode.compare("full"))
                    temp.delCol(0);
            }
            for(std::size_t k=j;k<W/2;++k)
            {
                if(mode.compare("valid"))
                {
                    temp.newCol(0);
                    temp.setCol(0, temp.getCol(1));
                }
                if(mode.compare("full"))
                    temp.delCol(temp.colNb()-1);
            }
            filtered_img+=temp;
        }
    }
    return filtered_img;
}


Matrix<float> interp1(Matrix<float> x, Matrix<float> v, Matrix<float> xq, std::string mode)
{
    std::size_t nb_elem_x = x.colNb();
    std::size_t nb_elem_xq = xq.colNb();
    Matrix<float> vq(1, nb_elem_xq);
    if(!mode.compare("linear"))
    {
        for(std::size_t i=0;i<nb_elem_xq;++i)
        {
            if(xq(0, i)<=x(0, 0))
                vq(0, i) = v(0, 0);
            else if(xq(0, i)>=x(0, nb_elem_x-1))
                vq(0, i) = v(0, nb_elem_x-1);
            else
            {
                std::size_t cur_idx = 1;
                while(xq(0, i)>=x(0, cur_idx))
                {
                    ++cur_idx;
                }
                vq(0, i) = v(0, cur_idx-1)+(xq(0,i)-x(0, cur_idx-1))*((v(0, cur_idx)-v(0, cur_idx-1))/(x(0, cur_idx)-x(0, cur_idx-1)));
            }
        }
    }
    else if(!mode.compare("nearest"))
    {
        for(std::size_t i=0;i<nb_elem_xq;++i)
        {
            Matrix<float> dist = abs(xq(0,i)-x);
            Matrix<std::size_t> idx = argmin(dist);
            vq(0,i) = v(0, idx(0,1));
        }
    }
    else if(!mode.compare("cubic"))
    {

    }
    return vq;
}

/*
pyramidal
1 2 3 2 1
2 4 6 4 2
3 6 9 6 3
2 4 6 4 2
1 2 3 2 1
*/
