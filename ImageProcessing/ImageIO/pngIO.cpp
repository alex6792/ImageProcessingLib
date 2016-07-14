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
    std::cout<<w<< " "<<h<<" "<<image.size()<<std::endl;

    if(error)
        return Matrix<Color>();
    else
    {
        Matrix<Color> img(h, w);
        for(unsigned int i=0;i<h;++i)
        {
            for(unsigned int j=0;j<w;++j)
                img(i, j) = Color(image[(i*w+j)*4], image[(i*w+j)*4+1], image[(i*w+j)*4+2], image[(i*w+j)*4+3]);
        }
        return img;
    }
}

void save_png(std::string filename, const Matrix<Color>& img)
{
    int W = img.colNb();
    int H = img.rowNb();
    unsigned char* data2 = (unsigned char*) malloc(4*H*W);
    int cpt = 0;
    for(int i=0;i<H;++i)
    {
        for(int j=0;j<W;++j)
        {
            data2[cpt++] = img(i, j).red();
            data2[cpt++] = img(i, j).green();
            data2[cpt++] = img(i, j).blue();
            data2[cpt++] = img(i, j).alpha();
        }
    }
    unsigned error = lodepng_encode32_file(filename.c_str(), data2, W, H);
}
