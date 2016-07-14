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
        char info[54];
        myfile.read(info, 54);

        std::string type = info;
        type = type.substr(0,4);
        int offset = *(int*)&info[10];
        int sizeof_header = *(int*)&info[14];
        int w = *(int*)&info[18];
        int h = *(int*)&info[22];
        int nb_plan = info[27]*256+info[26];
        int depth = info[29]*256+info[28];
        int compression = *(int*)&info[30];
        int total_size = *(int*)&info[34];
        int res_h = *(int*)&info[38];
        int res_v = *(int*)&info[42];
        int nb_color = *(int*)&info[46];
        int nb_important_color = *(int*)&info[50];


        std::cout<<"info"<<std::endl;
        std::cout<<type<<std::endl;
        std::cout<<offset<<std::endl;
        std::cout<<sizeof_header<<std::endl;
        std::cout<<w<<" "<<h<<std::endl;
        std::cout<<nb_plan<<std::endl;
        std::cout<<depth<<std::endl;
        std::cout<<compression<<std::endl;
        std::cout<<total_size<<std::endl;
        std::cout<<res_h<<std::endl;
        std::cout<<res_v<<std::endl;
        std::cout<<nb_color<<std::endl;
        std::cout<<nb_important_color<<std::endl;

        Matrix<Color> img(h, w);
        char palette[offset-54];
        myfile.read(palette, offset-54);
        char r,g,b;
        for(int i=h-1;i>=0;--i)
        {
            for(int j=0;j<w;++j)
            {
                if(depth == 24)
                {
                    myfile.get(b);
                    myfile.get(g);
                    myfile.get(r);
                    img(i, j) = Color(r, g, b);
                }
                else if(depth == 8)
                {
                    myfile.get(g);
                    if(offset==54)
                        img(i, j) = Color(g, g, g);
                    else
                        img(i, j) = Color(palette[4*(unsigned char)g+2], palette[4*(unsigned char)g+1], palette[4*(unsigned char)g+0]);
                }
                else if(depth == 1)
                {
                    if(j%8==0)
                    {
                        myfile.get(g);
                        unsigned char t = 128;
                        for(int k=0;k<8;++k)
                        {
                            bool temp = (unsigned char)g%(2*t)/t;
                            if(j+k<w)
                                img(i, j+k) = temp?White:Black;
                            t/=2;
                        }
                    }
                }
            }
            for(int j=0;j<(4-(w*depth/8)%4)%4;++j)
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
    std::ofstream myfile(filename, std::ios::out | std::ios::trunc |std::ios::binary);
    if(myfile)
    {
        int w = img.colNb();
        int h = img.rowNb();
        int filesize = 54 + 3*w*h;
        char bmpheader[14] = {'B', 'M',
                                filesize, filesize>>8, filesize>>16, filesize>>24,
                                0,0,0,0,
                                54,0,0,0};
        char bmpinfoheader[40] = {40, 0, 0, 0,
                                    w, w>>8, w>>16, w>>24,
                                    h, h>>8, h>>16, h>>24,
                                    1, 0, 24, 0};
        int  pad_size = (4-(w*3)%4)%4;
        char bmppad[3] = {0, 0, 0};
        myfile.write(bmpheader, 14);
        myfile.write(bmpinfoheader, 40);

        for(int i=h-1;i>=0;--i)
        {
            for(int j=0;j<w;++j)
            {
                Color c = img(i, j);
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
