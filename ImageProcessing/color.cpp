#include <cmath>
#include "color.hpp"


// constructors & accessors
Color::Color() : Color::Color(0, 0, 0, 255)
{
}

Color::Color(unsigned char r, unsigned char g, unsigned char b) : Color::Color(r, g, b, 255)
{
}

Color::Color(unsigned char r, unsigned char g, unsigned char b, unsigned char a)
{
    R = r;
    G = g;
    B = b;
    A = a;
}

unsigned char Color::red() const
{
    return R;
}

unsigned char Color::green() const
{
    return G;
}

unsigned char Color::blue() const
{
    return B;
}

unsigned char Color::alpha() const
{
    return A;
}

unsigned char Color::gray() const
{
    int g = (R+G+B)/3;
    return (unsigned char)g;
}

int Color::hue() const
{
    unsigned char minColor,maxColor;
    if(R>=G && R>=B)
        maxColor = R;
    else if(G>=R && G>=B)
        maxColor = G;
    else
        maxColor = B;
    if(R<=G && R<=B)
        minColor = R;
    else if(G<=R && G<=B)
        minColor = G;
    else
        minColor = B;

    float coeff = 60.0/(maxColor-minColor);
    if(maxColor==minColor)
        return 0;
    else if(maxColor==R)
        return (360+(int)coeff*(G-B))%360;
    else if(maxColor==G)
        return 120+coeff*(B-R);
    else
        return 240+coeff*(R-G);
}

int Color::saturation() const
{
    unsigned char minColor,maxColor;
    if(R>=G && R>=B)
        maxColor = R;
    else if(G>=R && G>=B)
        maxColor = G;
    else
        maxColor = B;
    if(R<=G && R<=B)
        minColor = R;
    else if(G<=R && G<=B)
        minColor = G;
    else
        minColor = B;

    if(minColor==maxColor)
        return 0;
    else
        return (maxColor-minColor)*100.0/maxColor;
}

unsigned char Color::value() const
{
    if(R>=G && R>=B)
        return R;
    else if(G>=R && G>=B)
        return G;
    else
        return B;
}

void Color::red(unsigned char r)
{
    R = r;
}

void Color::green(unsigned char g)
{
    G = g;
}

void Color::blue(unsigned char b)
{
    B = b;
}

void Color::alpha(unsigned char a)
{
    A = a;
}

void Color::HSV2RGB(int H, int S, unsigned char V)
{
    int h = (H/60)%6;
    float f = H/60.0-h;
    float l = float(V)*(1.0-S/100.0);
    float m = float(V)*(1.0-S*f/100.0);
    float n = float(V)*(1.0-1.0*S/100.0+f*S/100.0);
    if(h==0)
    {
        R = V;
        G = std::round(n);
        B = std::round(l);
    }
    else if(h==1)
    {
        R = std::round(m);
        G = V;
        B = std::round(l);
    }
    else if(h==2)
    {
        R = std::round(l);
        G = V;
        B = std::round(n);
    }
    else if(h==3)
    {
        R = std::round(l);
        G = std::round(m);
        B = V;
    }
    else if(h==4)
    {
        R = std::round(n);
        G = std::round(l);
        B = V;
    }
    else if(h==5)
    {
        R = V;
        G = std::round(l);
        B = std::round(m);
    }
}


// operators
void Color::operator+=(const Color& Col)
{
    R+=Col.R;
    G+=Col.G;
    B+=Col.B;
}

void Color::operator-=(const Color& Col)
{
    R-=Col.R;
    G-=Col.G;
    B-=Col.B;
}

void Color::operator*=(const Color& Col)
{
    R*=Col.R;
    G*=Col.G;
    B*=Col.B;
}

void Color::operator+=(const unsigned char& value)
{
    R+=value;
    G+=value;
    B+=value;
}

void Color::operator-=(const unsigned char& value)
{
    R-=value;
    G-=value;
    B-=value;
}

void Color::operator*=(const unsigned char& value)
{
    R*=value;
    G*=value;
    B*=value;
}

void Color::operator/=(const unsigned char& value)
{
    if(value!=0)
    {
        R/=value;
        G/=value;
        B/=value;
    }
}

Color Color::operator+(const Color& Col) const
{
    return Color(R+Col.R, G+Col.G, B+Col.B, A);
}

Color Color::operator-(const Color& Col) const
{
    return Color(R-Col.R, G-Col.G, B-Col.B, A);
}

Color Color::operator*(const Color& Col) const
{
    return Color(R*Col.R, G*Col.G, B*Col.B, A);
}

Color Color::operator+(const unsigned char& value) const
{
    return Color(R+value, G+value, B+value, A);
}

Color Color::operator-(const unsigned char& value) const
{
    return Color(R-value, G-value, B-value, A);
}

Color Color::operator*(const unsigned char& value) const
{
    return Color(R*value, G*value, B*value, A);
}

Color Color::operator/(const unsigned char& value) const
{
    return Color(R/value, G/value, B/value, A);
}

bool Color::operator==(const Color& Col) const
{
    return (R==Col.R && G==Col.G && B==Col.B);
}

bool Color::operator!=(const Color& Col) const
{
    return (R!=Col.R || G!=Col.G || B!=Col.B);
}



// operators
Color operator+(const unsigned char& value, const Color& Col)
{
    return Col+value;
}

Color operator-(const unsigned char& value, const Color& Col)
{
    return Color(value-Col.red(), value-Col.green(), value-Col.blue(), Col.alpha());
}

Color operator*(const unsigned char& value, const Color& Col)
{
    return Col*value;
}

std::ostream& operator<<(std::ostream& s, const Color& my_Color)
{
    s<<"["<<(int)my_Color.red()<<","<<(int)my_Color.green()<<","<<(int)my_Color.blue()<<"]";
    return s;
}
