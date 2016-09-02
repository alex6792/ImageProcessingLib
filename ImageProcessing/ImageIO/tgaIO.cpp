#include "tgaIO.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>


Matrix<Color> read_tga(std::string filename)
{
    std::ifstream myfile(filename, std::ios::in|std::ios::binary);
    if(myfile)
    {
        unsigned char header[18];
        myfile.read(reinterpret_cast<char*>(header), 18);

        unsigned idLength = header[0];
        unsigned colorMapType = header[1];
        unsigned compression = header[2];

        unsigned colorMapIndex = header[3] + header[4] * 256;
        unsigned colorMapLength = header[5] + header[6] * 256;
        unsigned colorMapSize = header[7];

        unsigned xOrigin = header[8] + header[9] * 256;
        unsigned yOrigin = header[10] + header[11] * 256;
        unsigned width = header[12] + header[13] * 256;
        unsigned height = header[14] + header[15] * 256;
        unsigned bpp = header[16];
        unsigned bytes = (bpp + 7) / 8;
        unsigned char alphaBits = header[17] & 0x0f;
        bool flipHoriz = (header[17] & 0x10) ? 1 : 0;
        bool flipVert = (header[17] & 0x20) ? 0 : 1;

        Matrix<Color> img(height, width);
        char r,g,b,a;
        auto it = img.begin();
        std::size_t n = img.size();
        if(compression == 2 && bpp==24)
        {
            for(std::size_t i=0;i<n;++i)
            {
                myfile.get(b);
                myfile.get(g);
                myfile.get(r);
                *it++ = Color(r, g, b);
            }
        }
        else if(compression == 2 && bpp==32)
        {
            for(std::size_t i=0;i<n;++i)
            {
                myfile.get(b);
                myfile.get(g);
                myfile.get(r);
                myfile.get(a);
                *it++ = Color(r, g, b, a);
            }
        }
        else if(compression == 10 && bpp==24)
        {
            auto it_end = img.end();
            while(it!=it_end)
            {
                char chunkheader;
                myfile.get(chunkheader);
                unsigned char chunk2 = chunkheader;
                if(chunk2<128)
                {
                    ++chunk2;
                    for(unsigned i=0;i<chunk2;++i)
                    {
                        myfile.get(b);
                        myfile.get(g);
                        myfile.get(r);
                        *it++ = Color(r,g,b);
                    }
                }
                else
                {
                    myfile.get(b);
                    myfile.get(g);
                    myfile.get(r);
                    chunk2-=127;
                    for(unsigned i=0;i<chunk2;++i)
                    {
                        *it++ = Color(r,g,b);
                    }
                }
            }
        }
        myfile.close();
        if(flipVert)
            img.flipud();
        if(flipHoriz)
            img.fliplr();
        return img;
    }
    std::cout<<"Unable to read the file : "<<filename<<std::endl;
    return Matrix<Color>();
}

void save_tga(std::string filename, const Matrix<Color>& img)
{
    std::ofstream myfile(filename, std::ios::out|std::ios::trunc|std::ios::binary);
    if(myfile)
    {
        unsigned w = img.colNb();
        unsigned h = img.rowNb();
        char seq[] = {0,0,2,0,0,0,0,0,0,0,0,0};
        myfile.write(seq, 12);
        myfile.put(w%256);
        myfile.put(w/256);
        myfile.put(h%256);
        myfile.put(h/256);
        myfile.put(32);
        myfile.put(32);
        std::for_each(img.cbegin(), img.cend(), [&myfile](const Color& c)
                                                {myfile.put(c.blue());
                                                myfile.put(c.green());
                                                myfile.put(c.red());
                                                myfile.put(c.alpha());});
        myfile.close();
    }
    else
        std::cout<<"unable to write in the file : "<<filename<<std::endl;
}
