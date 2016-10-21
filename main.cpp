#include "test_functions.hpp"
#include "matrix.hpp"
#include "statistics.hpp"
#include <numeric>
#include <algorithm>
#include <functional>
#include <deque>
#include "linalg.hpp"

float entropy(const Matrix<bool>& M)
{
    std::size_t s = count_nonzero(M);
    float p = float(s)/M.size();
    return -p*std::log2(p)-(1-p)*std::log2(1-p);
}

float entropy(const Matrix<unsigned char>& M)
{
    Matrix<std::size_t> hist = zeros<std::size_t>(1, 256);
    std::for_each(M.cbegin(), M.cend(), [&hist](unsigned char value){++hist(0, (std::size_t)value);});
    hist = compress(hist>0, hist, 2);
    Matrix<float> p(hist);
    p/=M.size();
    Matrix<float> logp = -log2(p);
    return sum(logp*p);
}



int main()
{
    /*Matrix<float> A = {
    {0.8147,    0.0975,    0.1576},
    {0.9058,    0.2785,    0.9706},
    {0.1270,    0.5469,    0.9572},
    {0.9134,    0.9575,    0.4854},
    {0.6324,    0.9649,    0.8003}};
    std::cout<<A<<std::endl;
    std::cout<<pinv(A)<<std::endl;
    std::cout<<dot(pinv(A), A)<<std::endl;*/

    Matrix<bool> M = {
                    { 1, 0, 0, 0, 0},
                    { 0, 0, 0, 0, 0}
                    };
    std::cout<<entropy(M)<<std::endl;;

    Matrix<unsigned char> N = {
                            { 1, 2, 3, 4, 5},
                            {1, 2, 3, 4, 5}
                            };
    std::cout<<entropy(N)<<std::endl;
    return 0;
}
