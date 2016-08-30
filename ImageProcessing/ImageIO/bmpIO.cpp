#include "bmpIO.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>


Matrix<Color> read_bmp(std::string filename)
{
    std::ifstream myfile(filename, std::ios::in|std::ios::binary);
    if(myfile)
    {
        unsigned char info[54];
        myfile.read(reinterpret_cast<char*>(info), 54);

        std::string type = reinterpret_cast<char*>(info);
        type = type.substr(0,4);
        unsigned offset = *(unsigned*)&info[10];
        unsigned sizeof_header = *(unsigned*)&info[14];
        unsigned w = *(unsigned*)&info[18];
        unsigned h = *(unsigned*)&info[22];
        unsigned nb_plan = info[27]*256+info[26];
        unsigned depth = info[29]*256+info[28];
        unsigned compression = *(unsigned*)&info[30];
        unsigned total_size = *(unsigned*)&info[34];
        unsigned res_h = *(unsigned*)&info[38];
        unsigned res_v = *(unsigned*)&info[42];
        unsigned nb_color = *(unsigned*)&info[46];
        unsigned nb_important_color = *(unsigned*)&info[50];

        Matrix<Color> img(h, w);
        unsigned char palette[offset-54];
        myfile.read(reinterpret_cast<char*>(palette), offset-54);
        char r,g,b;
        for(unsigned i=h;i>0;--i)
        {
            for(unsigned j=0;j<w;++j)
            {
                if(depth == 24)
                {
                    myfile.get(b);
                    myfile.get(g);
                    myfile.get(r);
                    img(i-1, j) = Color(r, g, b);
                }
                else if(depth == 16)
                {
                    char c1,c2;
                    myfile.get(c1);
                    myfile.get(c2);
                    unsigned idx = (unsigned char)c2*256+(unsigned char)c1;
                    if(offset!=54)
                        img(i-1, j) = Color(palette[4*idx+2], palette[4*idx+1], palette[4*idx]);
                }
                else if(depth == 8)
                {
                    myfile.get(g);
                    if(offset==54)
                        img(i-1, j) = Color(g, g, g);
                    else
                        img(i-1, j) = Color(palette[4*(unsigned char)g+2], palette[4*(unsigned char)g+1], palette[4*(unsigned char)g]);
                }
                else if(depth == 4)
                {
                    if(j%2==0)
                    {
                        myfile.get(g);
                        unsigned char idx1 = g/16;
                        unsigned char idx2 = g%16;
                        if(offset!=54)
                        {
                            img(i-1, j) = Color(palette[4*idx1+2], palette[4*idx1+1], palette[4*idx1]);
                            img(i-1, j+1) = Color(palette[4*idx2+2], palette[4*idx2+1], palette[4*idx2]);
                        }
                    }
                }
                else if(depth == 1)
                {
                    if(j%8==0)
                    {
                        myfile.get(g);
                        unsigned char t = 128;
                        for(unsigned k=0;k<8;++k)
                        {
                            bool temp = (unsigned char)g%(2*t)/t;
                            if(j+k<w)
                                img(i-1, j+k) = temp?White:Black;
                            t/=2;
                        }
                    }
                }
            }
            for(unsigned j=0;j<(4-(w*depth/8)%4)%4;++j)
                myfile.get(b);
        }
        myfile.close();
        return img;
    }
    std::cout<<"Unable to read the file : "<<filename<<std::endl;
    return Matrix<Color>();
}

void save_bmp(std::string filename, const Matrix<Color>& img)
{
    std::ofstream myfile(filename, std::ios::out|std::ios::trunc|std::ios::binary);
    if(myfile)
    {
        std::size_t w = img.colNb();
        std::size_t h = img.rowNb();
        std::size_t filesize = 54 + 3*w*h;
        char bmpheader[14] = {'B', 'M',
                                filesize, filesize>>8, filesize>>16, filesize>>24,
                                0,0,0,0,
                                54,0,0,0};
        char bmpinfoheader[40] = {40, 0, 0, 0,
                                    w, w>>8, w>>16, w>>24,
                                    h, h>>8, h>>16, h>>24,
                                    1, 0, 24, 0};
        unsigned  pad_size = (4-(w*3)%4)%4;
        char bmppad[3] = {0, 0, 0};
        myfile.write(bmpheader, 14);
        myfile.write(bmpinfoheader, 40);

        for(unsigned i=h;i>0;--i)
        {
            for(unsigned j=0;j<w;++j)
            {
                const Color& c = img(i-1, j);
                myfile.put(c.blue());
                myfile.put(c.green());
                myfile.put(c.red());
            }
            myfile.write(bmppad, pad_size);
        }
        myfile.close();
    }
    else
        std::cout<<"unable to write in the file : "<<filename<<std::endl;
}
