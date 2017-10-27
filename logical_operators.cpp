#include "matrix.hpp"
#include "logical_operators.hpp"


// all, any and none functions
bool all(const Matrix<bool>& M)
{
    return std::all_of(M.cbegin(), M.cend(), [](bool i){return i;});
}

bool any(const Matrix<bool>& M)
{
    return std::any_of(M.cbegin(), M.cend(), [](bool i){return i;});
}

bool none(const Matrix<bool>& M)
{
    return std::none_of(M.cbegin(), M.cend(), [](bool i){return i;});
}


// logical operators
bool NOT(bool b)
{
    return !b;
}

bool AND(bool a, bool b)
{
    return a&&b;
}

bool CDT(bool a, bool b)
{
    return (!a)||b;
}

bool EQ(bool a, bool b)
{
    return a==b;
}

bool NAND(bool a, bool b)
{
    return !(a&&b);
}

bool NOR(bool a, bool b)
{
    return !(a||b);
}

bool OR(bool a, bool b)
{
    return a||b;
}

bool XOR(bool a, bool b)
{
    return a!=b;
}

bool XNOR(bool a, bool b)
{
    return a==b;
}


// logical operators on matrix
Matrix<bool> NOT(const Matrix<bool>& M)
{
    return apply(M,NOT);
}

Matrix<bool> AND(const Matrix<bool>& A, const Matrix<bool>& B)
{
    return apply(A,B,AND);
}

Matrix<bool> CDT(const Matrix<bool>& A, const Matrix<bool>& B)
{
    return apply(A,B,CDT);
}

Matrix<bool> EQ(const Matrix<bool>& A, const Matrix<bool>& B)
{
    return apply(A,B,EQ);
}

Matrix<bool> NAND(const Matrix<bool>& A, const Matrix<bool>& B)
{
    return apply(A,B,NAND);
}

Matrix<bool> NOR(const Matrix<bool>& A, const Matrix<bool>& B)
{
    return apply(A,B,NOR);
}

Matrix<bool> OR(const Matrix<bool>& A, const Matrix<bool>& B)
{
    return apply(A,B,OR);
}

Matrix<bool> XOR(const Matrix<bool>& A, const Matrix<bool>& B)
{
    return apply(A,B,XOR);
}

Matrix<bool> XNOR(const Matrix<bool>& A, const Matrix<bool>& B)
{
    return apply(A,B,XNOR);
}
