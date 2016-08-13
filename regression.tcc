#include "linalg.hpp"
#include "statistics.hpp"


template <class T> Matrix<T> poly_regression(const Matrix<T>& Y, const Matrix<T>& X, std::size_t degree)
{
    if(X.colNb()==1 && X.rowNb()==Y.rowNb() && X.colNb()==Y.colNb())
    {
        if(degree>=0)
        {
            if(X.rowNb()>degree)
            {
                Matrix<T> Sn(degree+1);
                Sn(0, 0) = X.rowNb();
                Matrix<T> copy_matX = X;
                for(std::size_t i=1;i<=2*degree;++i)
                {
                    T cur_s = sum(copy_matX);
                    for(std::size_t j=0;j<degree+1;++j)
                    {
                        if(j<Sn.rowNb() && i-j<Sn.colNb() && i>=j)
                            Sn(j, i-j) = cur_s;
                    }
                    copy_matX*=X;
                }
                Matrix<T> Tn(degree+1, 1);
                Tn(0, 0) = sum(Y);
                Matrix<T> copy_matXY = Y;
                for(std::size_t j=1;j<degree+1;++j)
                {
                    copy_matXY*=X;
                    Tn(j, 0) = sum(copy_matXY);
                }

                std::tuple<Matrix<T>, Matrix<T>, Matrix<T> > LUP = lu(Sn);
                Matrix<T>& L = std::get<0>(LUP);
                Matrix<T>& U = std::get<1>(LUP);
                Matrix<T>& P = std::get<2>(LUP);
                Matrix<T> y = fwdsub(L, dot(P,Tn));
                return bwdsub(U, y);

                /*std::pair<Matrix<T>, Matrix<T> > DP = Jacobi(Sn);
                Matrix<T> D = DP.first;
                for(int k=0;k<D.rowNb();++k)
                    D(k, k) = 1/D(k, k);
                Matrix<T> P = DP.second;
                Matrix<T> Pt = DP.second;
                Pt.transpose();
                return dot(dot(P, D), dot(Pt,Tn));*/

                /*std::pair<Matrix<T>, Matrix<T> > LD = crout(Sn);
                Matrix<T>& L = LD.first;
                Matrix<T>& D = LD.second;
                Matrix<T> Lt = L;
                Lt.transpose();
                Matrix<T> y = fwdsub(L, Tn);
                return bwdsub(dot(D, Lt), y);*/

                /*std::pair<Matrix<T>, Matrix<T> > QR = qr(Sn);
                Matrix<T>& Q = QR.first;
                Matrix<T>& R = QR.second;
                Q.transpose();
                Matrix<T> y = dot(Q, Tn);
                return bwdsub(R, y);*/
            }
            else
                std::cout<<"not enough elements"<<std::endl;
        }
        else
            std::cout<<"invalid degree"<<std::endl;
    }
    else
        std::cout<<"invalid matrix size"<<std::endl;
    return Y;
}
