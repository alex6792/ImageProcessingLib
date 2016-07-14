#include "icoIO.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>


Matrix<Color> read_ico(std::string filename)
{
    std::ifstream myfile(filename, std::ios::in|std::ios::binary);
    if(myfile)
    {
        char header[6];
        myfile.read(header, 6);

        int reserved = header[0]+256*header[1];
        int image_type = header[2]+256*header[3];
        int nb_img = header[4]+256*header[5];
        std::cout<<reserved<<std::endl;
        std::cout<<image_type<<std::endl;
        std::cout<<nb_img<<std::endl;

        char image_entry[16];
        myfile.read(image_entry, 16);

        int w = image_entry[0];
        int h = image_entry[1];
        int nb_color = image_entry[2];
        reserved = image_entry[3];
        int color_plane = image_entry[4]+256*image_entry[5];
        int bpp = image_entry[6]+256*image_entry[7];
        int size = *(int*)&image_entry[8];
        int offset = *(int*)&image_entry[12];
        char info[40];
        myfile.read(info, 40);
        int sizeof_header = *(int*)&info[0];
        //w = *(int*)&info[4];
        //h = *(int*)&info[8];
        int nb_plan = info[13]*256+info[12];
        int depth = info[15]*256+info[14];
        int compression = *(int*)&info[16];
        int total_size = *(int*)&info[20];
        int res_h = *(int*)&info[24];
        int res_v = *(int*)&info[28];
        nb_color = *(int*)&info[32];
        int nb_important_color = *(int*)&info[32];
        std::cout<<"info"<<std::endl;
        std::cout<<sizeof_header<<std::endl;
        std::cout<<nb_plan<<std::endl;
        std::cout<<depth<<std::endl;
        std::cout<<compression<<std::endl;
        std::cout<<total_size<<std::endl;
        std::cout<<res_h<<std::endl;
        std::cout<<res_v<<std::endl;
        std::cout<<nb_color<<std::endl;
        std::cout<<nb_important_color<<std::endl;
        Matrix<Color> img(h, w);
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
                    //else
                        //img(i, j) = Color(palette[4*(unsigned char)g+2],palette[4*(unsigned char)g+1],palette[4*(unsigned char)g+0]);
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

                        /*img(i, j) = Color(g&0x04, g&0x04, g&0x04);
                        img(i, j+1) = Color(g&0x04, g&0x04, g&0x04);
                        img(i, j+2) = Color(g&0x04, g&0x04, g&0x04);
                        img(i, j+3) = Color(g&0x04, g&0x04, g&0x04);*/
                    }
                }
            }
            for(int j=0;j<(4-(w*3)%4)%4;++j)
                myfile.get(b);
        }


        myfile.close();
        return img;
    }
    std::cout<<"Unable to read the file : "<<filename<<std::endl;
    return Matrix<Color>();
}

void save_ico(std::string filename, const Matrix<Color>& img)
{

}
