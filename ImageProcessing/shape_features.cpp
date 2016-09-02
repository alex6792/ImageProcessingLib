#include "shape_features.hpp"
#include "morphology.hpp"
#include "../mmath.hpp"
#include "../statistics.hpp"


float circularity(const Matrix<bool>& M)
{
    float n1 = count(M, true);
    Matrix<bool> G = inner_gradient(M);
    float n2 = count(G, true);
    return 4*PI*n1/(n2*n2);
}


float major_axis(const Matrix<bool> & M)
{
    Matrix<float> M_copy(M);
    float result = 2*central_moment(M_copy,1,1)/(central_moment(M_copy,2,0)-central_moment(M_copy,0,2));
    return 0.5f*std::atan(result);
}
