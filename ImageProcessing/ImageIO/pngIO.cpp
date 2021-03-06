#include "pngIO.hpp"
#include "lodepng/lodepng.h"
#include <fstream>
#include <iostream>
#include <utility>
#include <string>
#include <sstream>


Matrix<Color> read_png(std::string filename)
{
    std::vector<unsigned char> buffer, image;
    lodepng::load_file(buffer, filename);
    unsigned w,h;
    unsigned error = lodepng::decode(image, w, h, buffer);

    if(error)
        return Matrix<Color>();
    else
    {
        Matrix<Color> img(h, w);
        unsigned char* image_ptr = image.data();
        std::for_each(img.begin(), img.end(), [&image_ptr](Color& c){c = Color(*image_ptr, *(image_ptr+1), *(image_ptr+2));image_ptr+=4;});
        return img;
    }
}

void save_png(std::string filename, const Matrix<Color>& img)
{
    std::size_t W = img.colNb();
    std::size_t H = img.rowNb();
    unsigned char* data2 = (unsigned char*) malloc(4*H*W);
    std::size_t cpt = 0;
    for(std::size_t i=0;i<H;++i)
    {
        for(std::size_t j=0;j<W;++j)
        {
            data2[cpt++] = img(i, j).red();
            data2[cpt++] = img(i, j).green();
            data2[cpt++] = img(i, j).blue();
            data2[cpt++] = img(i, j).alpha();
        }
    }
    unsigned error = lodepng_encode32_file(filename.c_str(), data2, W, H);
    if(error)
        std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
}
