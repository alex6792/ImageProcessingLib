#include "pcxIO.hpp"
#include <fstream>
#include <iostream>
#include <utility>
#include <string>
#include <sstream>


Matrix<Color> read_pcx(std::string filename)
{
    std::ifstream myfile(filename, std::ios::in|std::ios::binary);
    if(myfile)
    {
        unsigned char header[128];
        myfile.read(reinterpret_cast<char*>(header), 128);
        unsigned char id = header[0];
        unsigned char version = header[1];
        unsigned char encoding = header[2];
        unsigned char bpp = header[3];
        unsigned xstart = 256*header[5]+header[4];
        unsigned ystart = 256*header[7]+header[6];
        unsigned xend = 256*header[9]+header[8];
        unsigned yend = 256*header[11]+header[10];
        unsigned hres = 256*header[13]+header[12];
        unsigned vres = 256*header[15]+header[14];
        unsigned char* palette = &header[16];
        unsigned char reserved1 = header[64];
        unsigned char numbitplanes = header[65];
        unsigned bytesperline = 256*header[67]+header[66];
        unsigned palettetype = 256*header[69]+header[68];
        unsigned hscreensize = 256*header[71]+header[70];
        unsigned vscreensize = 256*header[73]+header[72];
        unsigned char* reserved2 = &header[74];
        Matrix<Color> img(yend-ystart+1, xend-xstart+1);

        char c, runcount, runvalue;
        std::size_t ScanLineLength = numbitplanes*bytesperline;
        std::size_t LinePaddingSize = ((numbitplanes*bytesperline)*(8/bpp))-((xend - xstart)+1);

        for(std::size_t i=0;i<yend-ystart+1;++i)
        {
            std::size_t idx = 0;
            while(idx < ScanLineLength)
            {
                myfile.get(c);
                if((c & 0xC0) == 0xC0)
                {
                    runcount = c & 0x3F;
                    myfile.get(runvalue);
                }
                else
                {
                    runcount = 1;
                    runvalue = c;
                }
                std::size_t total = 0;
                while(runcount && idx < ScanLineLength)
                {
                    if(idx<ScanLineLength/3)
                        img(i, idx).red(runvalue);
                    else if(idx<2*ScanLineLength/3)
                        img(i, idx-ScanLineLength/3).green(runvalue);
                    else
                        img(i, idx-2*ScanLineLength/3).blue(runvalue);

                    total += runcount;
                    --runcount, ++idx;
                }
            }
        }

        myfile.close();
        return img;
    }
    std::cout<<"Unable to read the file : "<<filename<<std::endl;
    return Matrix<Color>();
}

void save_pcx(std::string filename, const Matrix<Color>& img)
{
}
