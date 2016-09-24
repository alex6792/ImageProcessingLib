/*!
 * \author Alexandre Krebs
 * \file polynomial.hpp
 * \brief Polynomial
 */


#pragma once

#ifndef POLYNOMIAL_HPP
#define POLYNOMIAL_HPP


#include <iostream>
#include <vector>


/*!
 * \class Polynomial
 * \brief template class representing a polynomial
 */
class Polynomial
{
    private :
        std::vector<float> coefficients;
    public:
        Polynomial(const std::vector<float>&);

        // getter/setter
        float operator()(unsigned) const;
        float& operator()(unsigned);
        const std::vector<float>& get_coeffs() const;
        unsigned get_degree() const;

        // operators
        void operator+=(const Polynomial&);
        void operator-=(const Polynomial&);
        void operator*=(const Polynomial&);
        void operator/=(const Polynomial&);
        void operator%=(const Polynomial&);

        void operator+=(float);
        void operator-=(float);
        void operator*=(float);
        void operator/=(float);

        Polynomial operator+(const Polynomial&) const;
        Polynomial operator-(const Polynomial&) const;
        Polynomial operator*(const Polynomial&) const;
        Polynomial operator/(const Polynomial&) const;
        Polynomial operator%(const Polynomial&) const;

        Polynomial operator+(float) const;
        Polynomial operator-(float) const;
        Polynomial operator*(float) const;
        Polynomial operator/(float) const;

        Polynomial operator-() const;

        bool operator==(const Polynomial&) const;
        bool operator!=(const Polynomial&) const;
};


// operators
Polynomial operator+(float, const Polynomial&);
Polynomial operator-(float, const Polynomial&);
Polynomial operator*(float, const Polynomial&);

Polynomial Laguerre(unsigned);
Polynomial Legendre(unsigned);
Polynomial Tchebychev(unsigned);

//print
std::ostream& operator<<(std::ostream&, const Polynomial&);


#endif // POLYNOMIAL_HPP
