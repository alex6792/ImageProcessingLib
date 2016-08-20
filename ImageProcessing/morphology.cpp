#include <cmath>
#include <vector>
#include "../mmath.hpp"
#include "../logical_operators.hpp"
#include "morphology.hpp"


Mask circle(std::size_t r)
{
    return ellipse(r, r);
}

Mask cross(std::size_t a)
{
    Mask new_cross = zeros<bool>(2*a+1);
    new_cross.setCol(a, ones<bool>(new_cross.rowNb(), 1));
    new_cross.setRow(a, ones<bool>(1, new_cross.colNb()));
    return new_cross;
}

Mask ellipse(std::size_t a, std::size_t b)
{
    Mask new_ellipse = zeros<bool>(2*a+1, 2*b+1);
    for(std::size_t i=0;i<new_ellipse.rowNb();++i)
    {
        for(std::size_t j=0;j<new_ellipse.colNb();++j)
        {
            if(b*b*(i*i+a*a-2*i*a)+a*a*(j*j+b*b-2*j*b)<=a*a*b*b)
                new_ellipse(i, j) = true;
        }
    }
    return new_ellipse;
}

Mask diamond(std::size_t a)
{
    Mask new_diamond = zeros<bool>(a);
    for(std::size_t i=0;i<(a-1)/2;++i)
    {
        for(std::size_t j=-i+(a-1)/2;j<=i+a/2;++j)
            new_diamond(i, j) = true;
    }
    for(std::size_t i=a/2;i<a;++i)
    {
        for(std::size_t j=i-a/2;j<-i+3*a/2;++j)
            new_diamond(i, j) = true;
    }
    return new_diamond;
}

Mask rect(std::size_t a, std::size_t b)
{
    return ones<bool>(a, b);
}

Mask square(std::size_t a)
{
    return rect(a, a);
}


Matrix<bool> BTH(const Matrix<bool>& M, const Mask& mask)
{
    return close(M, mask)-M;
}

Matrix<bool> close(const Matrix<bool>& M, const Mask& mask)
{
    Matrix<bool> I = dilate(M, mask);
    Mask mask2 = mask;
    mask2.rot180();
    return erode(I, mask2);
}

Matrix<bool> closing_by_reconstruction(const Matrix<bool>& M, const Mask& mask)
{
    std::vector<Matrix<bool> > I;
    I.clear();
    I.push_back(dilate(M, mask));
    I.push_back(max(erode(I[0], mask), M));
    while(any(I[I.size()-1]!=I[I.size()-2]))
        I.push_back(max(erode(I[I.size()-1], mask), M));
    return I[I.size()-1];
}

Matrix<bool> conservative_smoothing(const Matrix<bool>& M)
{
    Mask N= {{1, 1, 1},
            {1, 0, 1},
            {1, 1, 1}};

    Matrix<bool> dil = dilate(M, N);
    Matrix<bool> ero = erode(M, N);
    Matrix<bool> I = max(min(M, dil), ero);
    return I;
}

Matrix<bool> dilate(const Matrix<bool>& M, const Mask& mask)
{
    std::size_t I = M.rowNb(), J = M.colNb();
    Matrix<bool> res = zeros<bool>(I, J);
    for(std::size_t k=0, K=mask.rowNb();k<K;++k)
    {
        for(std::size_t l=0, L=mask.colNb();l<L;++l)
        {
            if(mask(k, l))
            {
                std::size_t sup_bound_i = std::min(I, I+K/2-k);
                std::size_t inf_bound_i = sup_bound_i+k-I-K/2;
                std::size_t sup_bound_j = std::min(J, J+L/2-l);
                std::size_t inf_bound_j = sup_bound_j+l-J-L/2;
                for(std::size_t i=inf_bound_i;i<sup_bound_i;++i)
                {
                    for(std::size_t j=inf_bound_j;j<sup_bound_j;++j)
                    {
                        std::size_t x = i+k-K/2;
                        std::size_t y = j+l-L/2;
                        if(M(x, y))
                            res(i, j) = true;
                    }
                }
            }
        }
    }
    return res;
}

Matrix<bool> erode(const Matrix<bool>& M, const Mask& mask)
{
    std::size_t I = M.rowNb(), J = M.colNb();
    Matrix<bool> res = ones<bool>(I, J);
    for(std::size_t k=0, K=mask.rowNb();k<K;++k)
    {
        for(std::size_t l=0, L=mask.colNb();l<L;++l)
        {
            if(mask(k, l))
            {
                std::size_t sup_bound_i = std::min(I, I+K/2-k);
                std::size_t inf_bound_i = sup_bound_i+k-I-K/2;
                std::size_t sup_bound_j = std::min(J, J+L/2-l);
                std::size_t inf_bound_j = sup_bound_j+l-J-L/2;
                for(std::size_t i=inf_bound_i;i<sup_bound_i;++i)
                {
                    for(std::size_t j=inf_bound_j;j<sup_bound_j;++j)
                    {
                        std::size_t x = i+k-K/2;
                        std::size_t y = j+l-L/2;
                        if(!M(x, y))
                            res(i, j) = false;
                    }
                }
            }
        }
    }
    return res;
}

Matrix<bool> gradient(const Matrix<bool>& M, const Mask& mask)
{
    return dilate(M, mask)-erode(M, mask);
}

Matrix<bool> HitOrMiss(const Matrix<bool>& M, const Mask& maskFG, const Mask& maskBG)
{
    return AND(erode(M, maskFG), erode(NOT(M), maskBG));
}

Matrix<bool> inner_gradient(const Matrix<bool>& M, const Mask& mask)
{
    return M-erode(M, mask);
}

Matrix<bool> median_filter(const Matrix<bool>& M, const Mask& mask)
{
    Matrix<bool> res = ones<bool>(M.rowNb(), M.colNb());
    for(std::size_t i=0, I=M.rowNb();i<I;++i)
    {
        for(std::size_t j=0, J=M.colNb();j<J;++j)
        {
            std::size_t cpt_white = 0;
            std::size_t cpt_black = 0;
            for(std::size_t k=0, K=mask.rowNb();k<K;++k)
            {
                for(std::size_t l=0, L=mask.colNb();l<L;++l)
                {
                    if(k+i>=K/2 && k+i<I+K/2 && l+j>=L/2 && l+j<J+L/2)
                    {
                        std::size_t x = i+k-K/2;
                        std::size_t y = j+l-L/2;
                        if(M(x, y))
                            ++cpt_white;
                        else
                            ++cpt_black;
                    }
                }
            }
            res(i, j) = cpt_white>=cpt_black;
        }
    }
    return res;
}

Matrix<bool> open(const Matrix<bool>& M, const Mask& mask)
{
    Matrix<bool> I = erode(M, mask);
    Mask mask2 = mask;
    mask2.rot180();
    return dilate(I, mask2);
}

Matrix<bool> opening_by_reconstruction(const Matrix<bool>& M, const Mask& mask)
{
    std::vector<Matrix<bool> > I;
    I.clear();
    I.push_back(erode(M, mask));
    I.push_back(min(dilate(I[0], mask), M));
    while(any(I[I.size()-1]!=I[I.size()-2]))
        I.push_back(min(dilate(I[I.size()-1], mask), M));
    return I[I.size()-1];
}

Matrix<bool> outer_gradient(const Matrix<bool>& M, const Mask& mask)
{
    return dilate(M, mask)-M;
}

Matrix<bool> reconstruct(const Matrix<bool>& M, const Matrix<bool>& marker, const Mask& mask)
{
    if(M.rowNb()==marker.rowNb() && M.colNb()==marker.colNb())
    {
        std::vector<Matrix<bool> > cur_marker;
        cur_marker.clear();
        cur_marker.push_back(AND(marker, M));
        cur_marker.push_back(AND(dilate(cur_marker[cur_marker.size()-1], mask), M));
        while(any(cur_marker[cur_marker.size()-1]!=cur_marker[cur_marker.size()-2]))
            cur_marker.push_back(AND(dilate(cur_marker[cur_marker.size()-1], mask), M));
        return cur_marker[cur_marker.size()-1];
    }
    std::cout<<"dimension mismatch"<<std::endl;
    return M;
}

Matrix<bool> skeleton(const Matrix<bool>& M, const Mask& mask)
{
    std::vector<Matrix<bool> > list_erosion;
    list_erosion.clear();
    list_erosion.push_back(M);
    list_erosion.push_back(erode(M, mask));
    while(any((list_erosion[list_erosion.size()-1]!=list_erosion[list_erosion.size()-2])))
        list_erosion.push_back(erode(list_erosion[list_erosion.size()-1], mask));

    Mask mask2 = mask;
    mask2.rot180();
    Matrix<bool> S = M-dilate(list_erosion[1], mask2);
    for(std::size_t i=1;i<list_erosion.size()-1;++i)
        S+=list_erosion[i]-dilate(list_erosion[i+1], mask2);
    return S;
}

Matrix<bool> WTH(const Matrix<bool>& M, const Mask& mask)
{
    return M-open(M, mask);
}

Matrix<std::size_t> label(const Matrix<bool>& M, const Mask& mask)
{
    Matrix<std::size_t> labeledmap = zeros<std::size_t>(M.rowNb(), M.colNb());
    std::size_t curlabel = 0;
    for(std::size_t i=0;i<M.rowNb();++i)
    {
        for(std::size_t j=0;j<M.colNb();++j)
        {
            if(M(i, j) && labeledmap(i, j)==0)
            {
                ++curlabel;
                Matrix<bool> marker = zeros<bool>(M.rowNb(), M.colNb());
                marker(i, j) = true;
                marker = reconstruct(M, marker, mask);
                labeledmap = where(marker, full(M.rowNb(), M.colNb(), curlabel), labeledmap);
            }
        }
    }
    return labeledmap;
}

Matrix<unsigned char> BTH(const Matrix<unsigned char>& M, const Mask& mask)
{
    return close(M, mask)-M;
}

Matrix<unsigned char> close(const Matrix<unsigned char>& M, const Mask& mask)
{
    Matrix<unsigned char> I = dilate(M, mask);
    Mask mask2 = mask;
    mask2.rot180();
    return erode(I, mask2);
}

Matrix<unsigned char> closing_by_reconstruction(const Matrix<unsigned char>& M, const Mask& mask)
{
    std::vector<Matrix<unsigned char> > I;
    I.clear();
    I.push_back(dilate(M, mask));
    I.push_back(max(erode(I[0], mask), M));
    while(any(I[I.size()-1]!=I[I.size()-2]))
        I.push_back(max(erode(I[I.size()-1], mask), M));
    return I[I.size()-1];
}

Matrix<unsigned char> conservative_smoothing(const Matrix<unsigned char>& M)
{
    Mask N= {{1, 1, 1},
            {1, 0, 1},
            {1, 1, 1}};

    Matrix<unsigned char> dil = dilate(M, N);
    Matrix<unsigned char> ero = erode(M, N);
    Matrix<unsigned char> I = max(min(M, dil), ero);
    return I;
}

Matrix<unsigned char> erode(const Matrix<unsigned char>& M, const Mask& mask)
{
    std::size_t I = M.rowNb(), J = M.colNb();
    Matrix<unsigned char> res = full<unsigned char>(I, J, 255);
    for(std::size_t k=0, K=mask.rowNb();k<K;++k)
    {
        for(std::size_t l=0, L=mask.colNb();l<L;++l)
        {
            if(mask(k, l))
            {
                std::size_t sup_bound_i = std::min(I, I+K/2-k);
                std::size_t inf_bound_i = sup_bound_i+k-I-K/2;
                std::size_t sup_bound_j = std::min(J, J+L/2-l);
                std::size_t inf_bound_j = sup_bound_j+l-J-L/2;
                for(std::size_t i=inf_bound_i;i<sup_bound_i;++i)
                {
                    for(std::size_t j=inf_bound_j;j<sup_bound_j;++j)
                    {
                        std::size_t x = i+k-K/2;
                        std::size_t y = j+l-L/2;
                        if(M(x, y)<res(i, j))
                            res(i, j) = M(x, y);
                    }
                }
            }
        }
    }
    return res;
}

Matrix<unsigned char> dilate(const Matrix<unsigned char>& M, const Mask& mask)
{
    std::size_t I = M.rowNb(), J = M.colNb();
    Matrix<unsigned char> res = zeros<unsigned char>(I, J);
    for(std::size_t k=0, K=mask.rowNb();k<K;++k)
    {
        for(std::size_t l=0, L=mask.colNb();l<L;++l)
        {
            if(mask(k, l))
            {
                std::size_t sup_bound_i = std::min(I, I+K/2-k);
                std::size_t inf_bound_i = sup_bound_i+k-I-K/2;
                std::size_t sup_bound_j = std::min(J, J+L/2-l);
                std::size_t inf_bound_j = sup_bound_j+l-J-L/2;
                for(std::size_t i=inf_bound_i;i<sup_bound_i;++i)
                {
                    for(std::size_t j=inf_bound_j;j<sup_bound_j;++j)
                    {
                        std::size_t x = i+k-K/2;
                        std::size_t y = j+l-L/2;
                        if(M(x, y)>res(i, j))
                            res(i, j) = M(x, y);
                    }
                }
            }
        }
    }
    return res;
}

Matrix<unsigned char> gradient(const Matrix<unsigned char>& M, const Mask& mask)
{
    return dilate(M, mask)-erode(M, mask);
}

Matrix<unsigned char> HitOrMiss(const Matrix<unsigned char>& M, const Mask& maskFG, const Mask& maskBG)
{
    return min(erode(M, maskFG), erode((unsigned char)255-M, maskBG));
}

Matrix<unsigned char> inner_gradient(const Matrix<unsigned char>& M, const Mask& mask)
{
    return M-erode(M, mask);
}

Matrix<unsigned char> median_filter(const Matrix<unsigned char>& M, const Mask& mask)
{
    Matrix<unsigned char> res(M.rowNb(), M.colNb());
    for(std::size_t i=0, I=M.rowNb();i<I;++i)
    {
        for(std::size_t j=0, J=M.colNb();j<J;++j)
        {
            std::vector<unsigned char> vec;
            vec.clear();
            for(std::size_t k=0, K=mask.rowNb();k<K;++k)
            {
                for(std::size_t l=0, L=mask.colNb();l<L;++l)
                {
                    if(k+i>=K/2 && k+i<I+K/2 && l+j>=L/2 && l+j<J+L/2)
                    {
                        std::size_t x = i+k-K/2;
                        std::size_t y = j+l-L/2;
                        vec.push_back(M(x, y));
                    }
                }
            }
            std::sort(vec.begin(), vec.end());
            if(vec.size()%2==1)
                res(i, j) = vec[(vec.size()-1)/2];
            else
                res(i, j) = ((int)vec[vec.size()/2-1]+(int)vec[vec.size()/2])/2;
        }
    }
    return res;
}

Matrix<unsigned char> open(const Matrix<unsigned char>& M, const Mask& mask)
{
    Matrix<unsigned char> I = erode(M, mask);
    Mask mask2 = mask;
    mask2.rot180();
    return dilate(I, mask2);
}

Matrix<unsigned char> opening_by_reconstruction(const Matrix<unsigned char>& M, const Mask& mask)
{
    std::vector<Matrix<unsigned char> > I;
    I.clear();
    I.push_back(erode(M, mask));
    I.push_back(min(dilate(I[0], mask), M));
    while(any(I[I.size()-1]!=I[I.size()-2]))
        I.push_back(min(dilate(I[I.size()-1], mask), M));
    return I[I.size()-1];
}

Matrix<unsigned char> outer_gradient(const Matrix<unsigned char>& M, const Mask& mask)
{
    return dilate(M, mask)-M;
}

Matrix<unsigned char> reconstruct(const Matrix<unsigned char>& M, const Matrix<unsigned char>& marker, const Mask& mask)
{
    if(M.rowNb()==marker.rowNb() && M.colNb()==marker.colNb())
    {
        std::vector<Matrix<unsigned char> > cur_marker;
        cur_marker.clear();
        cur_marker.push_back(min(marker, M));
        cur_marker.push_back(min(dilate(cur_marker[cur_marker.size()-1], mask), M));
        while(any(cur_marker[cur_marker.size()-1]!=cur_marker[cur_marker.size()-2]))
            cur_marker.push_back(min(dilate(cur_marker[cur_marker.size()-1], mask), M));
        return cur_marker[cur_marker.size()-1];
    }
    std::cout<<"dimension mismatch"<<std::endl;
    return M;
}

Matrix<unsigned char> skeleton(const Matrix<unsigned char>& M, const Mask& mask)
{
    std::vector<Matrix<unsigned char> > list_erosion;
    list_erosion.clear();
    list_erosion.push_back(M);
    list_erosion.push_back(erode(M, mask));
    while(any((list_erosion[list_erosion.size()-1]!=list_erosion[list_erosion.size()-2])))
        list_erosion.push_back(erode(list_erosion[list_erosion.size()-1], mask));

    Mask mask2 = mask;
    mask2.rot180();
    Matrix<unsigned char> S = M-dilate(list_erosion[1], mask2);
    for(std::size_t i=1;i<list_erosion.size()-1;++i)
        S+=list_erosion[i]-dilate(list_erosion[i+1], mask2);
    return S;
}

Matrix<unsigned char> WTH(const Matrix<unsigned char>& M, const Mask& mask)
{
    return M-open(M, mask);
}

Matrix<std::size_t> label(const Matrix<unsigned char>& M, const Mask& mask)
{
    Matrix<std::size_t> labeledmap = zeros<std::size_t>(M.rowNb(), M.colNb());
    std::size_t curlabel = 0;
    for(std::size_t i=0;i<M.rowNb();++i)
    {
        for(std::size_t j=0;j<M.colNb();++j)
        {
            if(M(i, j) && labeledmap(i, j)==0)
            {
                ++curlabel;
                Matrix<unsigned char> marker = zeros<unsigned char>(M.rowNb(), M.colNb());
                marker(i, j) = 255;
                marker = reconstruct(M, marker, mask);
                labeledmap = where(marker>0, full<std::size_t>(M.rowNb(), M.colNb(), curlabel), labeledmap);
            }
        }
    }
    return labeledmap;
}
