#include "polynomial.hpp"
#include <algorithm>
#include <functional>


Polynomial::Polynomial(const std::vector<float>& v)
{
    coefficients = v;
}

float Polynomial::operator()(unsigned i) const
{
    return coefficients[i];
}

float& Polynomial::operator()(unsigned i)
{
    return coefficients[i];
}

const std::vector<float>& Polynomial::get_coeffs() const
{
    return coefficients;
}

unsigned Polynomial::get_degree() const
{
    return coefficients.size()-1;
}


// operators
void Polynomial::operator+=(const Polynomial& p)
{
    unsigned dp = p.get_degree();
    unsigned dthis = get_degree();
    auto p_coeffs = p.get_coeffs();
    if(dthis<dp)
        coefficients.resize(p.get_degree()+1, 0.0f);
    std::transform(p_coeffs.cbegin(), p_coeffs.cend(), coefficients.cbegin(), coefficients.begin(), std::plus<float>());
}

void Polynomial::operator-=(const Polynomial& p)
{
    unsigned dp = p.get_degree();
    unsigned dthis = get_degree();
    auto p_coeffs = p.get_coeffs();
    if(dthis<dp)
        coefficients.resize(p.get_degree()+1, 0.0f);
    std::transform(coefficients.cbegin(), coefficients.cbegin()+dp+1, p_coeffs.cbegin(), coefficients.begin(), std::minus<float>());
}

void Polynomial::operator*=(const Polynomial& p)
{
    unsigned this_degree = this->get_degree();
    unsigned p_degree = p.get_degree();
    std::vector<float> new_coeffs(this_degree+p_degree+1, 0.0f);
    for(unsigned i=0;i<=this_degree;++i)
    {
        for(unsigned j=0;j<=p_degree;++j)
        {
            new_coeffs[i+j]+=p(j)*coefficients[i];
        }
    }
    coefficients = new_coeffs;
}

void Polynomial::operator/=(const Polynomial& p)
{
}

void Polynomial::operator%=(const Polynomial& p)
{
}

void Polynomial::operator+=(float v)
{
    coefficients[0]+=v;
}

void Polynomial::operator-=(float v)
{
    coefficients[0]-=v;
}

void Polynomial::operator*=(float v)
{
    std::for_each(coefficients.begin(), coefficients.end(), [v](float& c){c*=v;});
}

void Polynomial::operator/=(float v)
{
    std::for_each(coefficients.begin(), coefficients.end(), [v](float& c){c/=v;});
}

Polynomial Polynomial::operator+(const Polynomial& p) const
{
    Polynomial p_new = *this;
    p_new+=p;
    return p_new;
}

Polynomial Polynomial::operator-(const Polynomial& p) const
{
    Polynomial p_new = *this;
    p_new-=p;
    return p_new;
}

Polynomial Polynomial::operator*(const Polynomial& p) const
{
    Polynomial p_new = *this;
    p_new*=p;
    return p_new;
}

Polynomial Polynomial::operator/(const Polynomial& p) const
{
    Polynomial p_new = *this;
    p_new/=p;
    return p_new;
}

Polynomial Polynomial::operator%(const Polynomial& p) const
{
    Polynomial p_new = *this;
    p_new%=p;
    return p_new;
}

Polynomial Polynomial::operator+(float v) const
{
    Polynomial p_new = *this;
    p_new+=v;
    return p_new;
}

Polynomial Polynomial::operator-(float v) const
{
    Polynomial p_new = *this;
    p_new-=v;
    return p_new;
}

Polynomial Polynomial::operator*(float v) const
{
    Polynomial p_new = *this;
    p_new*=v;
    return p_new;
}

Polynomial Polynomial::operator/(float v) const
{
    Polynomial p_new = *this;
    p_new/=v;
    return p_new;
}

Polynomial Polynomial::operator-() const
{

    return (*this)*(-1.0f);
}

bool Polynomial::operator==(const Polynomial& p) const
{
    if(this->get_degree()!=p.get_degree())
        return false;
    else
    {
        auto& coeffs1 = this->get_coeffs();
        auto& coeffs2 = p.get_coeffs();
        //return std::all_of(coeffs1.cbegin(),coeff2.begin)
    }
}

bool Polynomial::operator!=(const Polynomial& p) const
{
    return !(*this==p);
}


// operators
Polynomial operator+(float v, const Polynomial& p)
{
    return p+v;
}

Polynomial operator-(float v, const Polynomial& p)
{
    return -p+v;
}

Polynomial operator*(float v, const Polynomial& p)
{
    return p*v;
}

Polynomial Legendre(unsigned i)
{
    if(i==0)
        return Polynomial({1.0f});
    else if(i==1)
        return Polynomial({0.0f, 1.0f});
    else
        return (Legendre(i-1)*Polynomial({0.0f, 2.0f*(i-1)+1.0f})-Legendre(i-2)*float(i-1))/float(i);
}

Polynomial Laguerre(unsigned i)
{
    if(i==0)
        return Polynomial({1.0f});
    else if(i==1)
        return Polynomial({1.0f, -1.0f});
    else
        return (Laguerre(i-1)*Polynomial({2.0f*(i-1)+1.0f, -1.0f})-Laguerre(i-2)*float(i-1))/float(i);
}

Polynomial Tchebychev(unsigned i)
{
    if(i==0)
        return Polynomial({1.0f});
    else if(i==1)
        return Polynomial({0.0f, 1.0f});
    else
        return Tchebychev(i-1)*Polynomial({0.0f, 2.0f})-Tchebychev(i-2);
}

//print
std::ostream& operator<<(std::ostream& s, const Polynomial& p)
{
    auto& coeffs = p.get_coeffs();
    s<<coeffs[0];
    if(coeffs.size()>=2 && coeffs[1]>0.0f)
        s<<" + "<<coeffs[1]<<"X";
    else if(coeffs.size()>=2)
        s<<" "<<coeffs[1]<<"X";

    for(int i=2;i<coeffs.size();++i)
    {
        if(coeffs[i]==0.0f)
            continue;
        else if(coeffs[i]==1.0f)
            s<<" + "<<"X^"<<i;
        else if(coeffs[i]==-1.0f)
            s<<" - "<<"X^"<<i;
        else if(coeffs[i]>0.0f)
            s<<" + "<<coeffs[i]<<"X^"<<i;
        else
            s<<" "<<coeffs[i]<<"X^"<<i;
    }
    return s;
}
