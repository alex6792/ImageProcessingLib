#include "../imgconverter.hpp"
#include "imageIO.hpp"
#include "bmpIO.hpp"
#include "gifIO.hpp"
#include "icoIO.hpp"
#include "jpegIO.hpp"
#include "pngIO.hpp"
#include "pnmIO.hpp"
#include "tgaIO.hpp"
#include "tiffIO.hpp"


Matrix<Color> read_img(std::string filename)
{
    std::size_t n = filename.size();
    std::string extension3 = filename.substr(n-3, n);
    std::string extension4 = filename.substr(n-4, n);
    if(!extension3.compare("bmp"))
    {
        return read_bmp(filename);
    }
    else if(!extension3.compare("gif"))
    {
        return read_gif(filename);
    }
    else if(!extension3.compare("ico"))
    {
        return read_ico(filename);
    }
    else if(!extension3.compare("jpg"))
    {
        return read_jpeg(filename);
    }
    else if(!extension3.compare("png"))
    {
        return read_png(filename);
    }
    else if(!extension3.compare("pbm"))
    {
        return bw2colorimage(read_pbm(filename));
    }
    else if(!extension3.compare("pgm"))
    {
        return gray2colorimage(read_pgm(filename));
    }
    else if(!extension3.compare("ppm"))
    {
        return read_ppm(filename);
    }
    else if(!extension3.compare("tga"))
    {
        return read_tga(filename);
    }
    else if(!extension4.compare("tiff"))
    {
        return read_tiff(filename);
    }
    else
    {
        return Matrix<Color>();
    }
}


void save_img(std::string filename, const Matrix<Color>& img)
{
    std::size_t n = filename.size();
    std::string extension3 = filename.substr(n-3, n);
    std::string extension4 = filename.substr(n-4, n);
    if(!extension3.compare("bmp"))
    {
        save_bmp(filename, img);
    }
    else if(!extension3.compare("gif"))
    {
        save_gif(filename, img);
    }
    else if(!extension3.compare("ico"))
    {
        save_ico(filename, img);
    }
    else if(!extension3.compare("jpg"))
    {
        save_jpeg(filename, img);
    }
    else if(!extension3.compare("png"))
    {
        save_png(filename, img);
    }
    else if(!extension3.compare("pbm"))
    {
        save_pbm(filename, color2bwimage(img));
    }
    else if(!extension3.compare("pgm"))
    {
        save_pgm(filename, color2grayimage(img));
    }
    else if(!extension3.compare("ppm"))
    {
        save_ppm(filename, img);
    }
    else if(!extension3.compare("tga"))
    {
        save_tga(filename, img);
    }
    else if(!extension4.compare("tiff"))
    {
        save_tiff(filename, img);
    }
}
