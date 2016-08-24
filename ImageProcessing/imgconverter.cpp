#include <cfloat>
#include <cmath>
#include "imgconverter.hpp"
#include "morphology.hpp"
#include "../statistics.hpp"


Matrix<bool> gray2bwimage(const Matrix<unsigned char>& img, unsigned char thresh)
{
    return img>=thresh;
}

Matrix<bool> color2bwimage(const Matrix<Color>& img, unsigned char thresh)
{
    return apply(img, &Color::gray)>=thresh;
}

Matrix<unsigned char> bw2grayimage(const Matrix<bool>& img)
{
    return where(img, (unsigned char)255, (unsigned char)0);
}

Matrix<unsigned char> color2grayimage(const Matrix<Color>& img)
{
    return apply(img, &Color::gray);
}

Matrix<Color> bw2colorimage(const Matrix<bool>& img)
{
    return where(img, White, Black);
}

Matrix<Color> gray2colorimage(const Matrix<unsigned char>& img)
{
    return apply<unsigned char, Color>(img, [](unsigned char g){return Color(g, g, g);});
}

Matrix<Color> array2colorimage(const Matrix<std::size_t>& img)
{
    std::size_t minimum = min(img), maximum = max(img);
    Matrix<Color> new_img(img.rowNb(), img.colNb());
    Matrix<std::size_t> normalized_img = 360u*(img-minimum)/(maximum-minimum);
    auto it = normalized_img.cbegin();
    std::for_each(new_img.begin(),
                  new_img.end(),
                  [&it](Color& c){c.HSV2RGB(*it, 255, 127);++it;});
    return new_img;
}

Matrix<std::size_t> histogram(const Matrix<bool>& img)
{
    Matrix<std::size_t> hist(1, 2);
    hist(0, 0) = count(img, false);
    hist(0, 1) = img.size()-hist(0, 0);
    return hist;
}

Matrix<std::size_t> histogram(const Matrix<unsigned char>& img)
{
    Matrix<std::size_t> hist = zeros<std::size_t>(1, 256);
    std::for_each(img.cbegin(), img.cend(), [&hist](unsigned char value){++hist(0, (std::size_t)value);});
    return hist;
}

Matrix<std::size_t> histogram(const Matrix<std::size_t>& img)
{
    std::size_t maximum = max(img);
    Matrix<std::size_t> hist = zeros<std::size_t>(1, maximum+1);
    std::for_each(img.cbegin(), img.cend(), [&hist](unsigned char value){++hist(0, value);});
    return hist;
}
