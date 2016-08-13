#include <cfloat>
#include <cmath>
#include "imgconverter.hpp"
#include "morphology.hpp"
#include "../statistics.hpp"


Matrix<bool> gray2bwimage(const Matrix<unsigned char>& img, unsigned char thresh)
{
    return img>=thresh;
}

Matrix<bool> hysteresis(const Matrix<unsigned char>& img, unsigned char low_thresh, unsigned char high_thresh)
{
    return reconstruct(img>=low_thresh, img>=high_thresh);
}

Matrix<bool> otsu(const Matrix<unsigned char>& img)
{
    Matrix<std::size_t> hist = histogram(img);

    // BG
    Matrix<float> hist_normalized = Matrix<float>(hist)/float(img.size());
    Matrix<float> cum_hist = cumsum(hist_normalized);
    auto XY = meshgrid<float>(hist.rowNb(), hist.colNb());
    Matrix<float>& Y = XY.second;
    Matrix<float> mean_BG = hist_normalized*Y;
    Matrix<float> var_BG = mean_BG*Y;
    mean_BG = cumsum(mean_BG);
    var_BG = cumsum(var_BG);
    std::transform(mean_BG.cbegin(), mean_BG.cend(), cum_hist.cbegin(), mean_BG.begin(), std::divides<float>());
    std::transform(var_BG.cbegin(), var_BG.cend(), cum_hist.cbegin(), var_BG.begin(), std::divides<float>());
    var_BG-=mean_BG*mean_BG;

    // FG
    hist_normalized.fliplr();
    Matrix<float> cum_hist_inv = cumsum(hist_normalized);
    Y.fliplr();
    Matrix<float> mean_FG = hist_normalized*Y;
    Matrix<float> var_FG = mean_FG*Y;
    mean_FG = cumsum(mean_FG);
    var_FG = cumsum(var_FG);
    std::transform(mean_FG.cbegin(), mean_FG.cend(), cum_hist_inv.cbegin(), mean_FG.begin(), std::divides<float>());
    std::transform(var_FG.cbegin(), var_FG.cend(), cum_hist_inv.cbegin(), var_FG.begin(), std::divides<float>());
    var_FG-=mean_FG*mean_FG;
    cum_hist_inv.fliplr();
    mean_FG.fliplr();
    var_FG.fliplr();
    cum_hist.delCol(cum_hist.colNb()-1);
    var_BG.delCol(var_BG.colNb()-1);
    cum_hist_inv.delCol(0);
    var_FG.delCol(0);

    Matrix<float> weigthed_sum = cum_hist*var_BG+cum_hist_inv*var_FG;
    std::replace_if(weigthed_sum.begin(), weigthed_sum.end(), [](float value){return std::isnan(value);}, FLT_MAX);
    Matrix<std::size_t> results = argmin(weigthed_sum);
    return img>results(0, 1);
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
