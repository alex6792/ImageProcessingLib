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
        char header_c[18];
        myfile.read(header_c, 18);
        unsigned char* header = reinterpret_cast<unsigned char*>(header_c);

        int idLength = header[0];
        int colorMapType = header[1];
        int compression = header[2];

        int colorMapIndex = header[3] + header[4] * 256;
        int colorMapLength = header[5] + header[6] * 256;
        int colorMapSize = header[7];

        int xOrigin = header[8] + header[9] * 256;
        int yOrigin = header[10] + header[11] * 256;
        int width = header[12] + header[13] * 256;
        int height = header[14] + header[15] * 256;
        int bpp = header[16];
        int bytes = (bpp + 7) / 8;
        unsigned char alphaBits = header[17] & 0x0f; /* Just the low 4 bits */
        bool flipHoriz = (header[17] & 0x10) ? 1 : 0;
        bool flipVert = (header[17] & 0x20) ? 0 : 1;
        std::cout<<width<<" "<<height<<std::endl;
        std::cout<<compression<<std::endl;
        std::cout<<bpp<<std::endl;
        std::cout<<flipHoriz<<std::endl;
        std::cout<<flipVert<<std::endl;
        Matrix<Color> img(height, width);
        char r,g,b,a;
        if(compression == 2 && bpp==24)
        {
            for(int i=0;i<height;++i)
            {
                for(int j=0;j<width;++j)
                {
                    myfile.get(b);
                    myfile.get(g);
                    myfile.get(r);
                    img(i, j) = Color(r, g, b);
                }
            }
        }
        if(compression == 2 && bpp==32)
        {
            for(int i=0;i<height;++i)
            {
                for(int j=0;j<width;++j)
                {
                    myfile.get(b);
                    myfile.get(g);
                    myfile.get(r);
                    myfile.get(a);
                    img(i, j) = Color(r, g, b, a);
                }
            }
        }
        else if(compression == 10)
        {
            std::size_t cpt = 0;
            while(cpt<img.size())
            {
                char chunkheader;
                myfile.get(chunkheader);
                if((unsigned char)chunkheader<128)
                {
                    ++chunkheader;													// add 1 to get number of following color values
                    for(int counter = 0; counter < (unsigned char)chunkheader; counter++)		// Read RAW color values
                    {
                        myfile.get(b);
                        myfile.get(g);
                        myfile.get(r);
                        img(cpt/width, cpt%width) = Color(r,g,b);
                        ++cpt;														// Return failed
                    }
                }
                else
                {
                    unsigned char chunk2 = chunkheader;
                    myfile.get(b);
                    myfile.get(g);
                    myfile.get(r);
                    for(int i=0;i<chunk2 - 127;++i)
                    {
                        img(cpt/width, cpt%width) = Color(r,g,b);
                        ++cpt;
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
    std::ofstream myfile(filename, std::ios::out | std::ios::trunc |std::ios::binary);
    if(myfile)
    {
        int w = img.colNb();
        int h = img.rowNb();
        char seq[] = {0,0,2,0,0,0,0,0,0,0,0,0};
        myfile.write(seq, 12);
        myfile.put(w%256);
        myfile.put(w/256);
        myfile.put(h%256);
        myfile.put(h/256);
        myfile.put(32);myfile.put(32);
        std::for_each(img.cbegin(), img.cend(), [&myfile](Color c){myfile.put(c.blue());
                                                                    myfile.put(c.green());
                                                                    myfile.put(c.red());
                                                                    myfile.put(c.alpha());});
        myfile.close();
    }
    else
        std::cout<<"unable to write in the file : "<<filename<<std::endl;
}
