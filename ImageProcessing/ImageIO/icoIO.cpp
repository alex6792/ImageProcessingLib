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
        unsigned char header[6];
        myfile.read(reinterpret_cast<char*>(header), 6);
        unsigned reserved = header[0]+256*header[1];
        unsigned image_type = header[2]+256*header[3];
        unsigned nb_img = header[4]+256*header[5];
        std::cout<<reserved<<std::endl;
        std::cout<<image_type<<std::endl;
        std::cout<<nb_img<<std::endl;

        unsigned char image_entry[16];
        myfile.read(reinterpret_cast<char*>(image_entry), 16);
        unsigned w = image_entry[0];
        unsigned h = image_entry[1];
        unsigned nb_color = image_entry[2];
        reserved = image_entry[3];
        unsigned color_plane = image_entry[4]+256*image_entry[5];
        unsigned bpp = image_entry[6]+256*image_entry[7];
        unsigned size = *(unsigned*)&image_entry[8];
        unsigned offset = *(unsigned*)&image_entry[12];

        std::cout<<"image_entry"<<std::endl;
        for(int i=0;i<4;++i)
        {
            std::cout<<(*(unsigned*)&image_entry[4*i])<<std::endl;
        }

        unsigned char info[40];
        myfile.read(reinterpret_cast<char*>(info), 40);
        unsigned sizeof_header = *(unsigned*)&info[0];
        //w = *(unsigned*)&info[4];
        //h = *(unsigned*)&info[8];
        unsigned nb_plan = info[13]*256+info[12];
        unsigned depth = info[15]*256+info[14];
        unsigned compression = *(unsigned*)&info[16];
        unsigned total_size = *(unsigned*)&info[20];
        unsigned res_h = *(unsigned*)&info[24];
        unsigned res_v = *(unsigned*)&info[28];
        nb_color = *(unsigned*)&info[32];
        unsigned nb_important_color = *(unsigned*)&info[36];
        /*std::cout<<"info"<<std::endl;
        std::cout<<sizeof_header<<std::endl;
        std::cout<<nb_plan<<std::endl;
        std::cout<<depth<<std::endl;
        std::cout<<compression<<std::endl;
        std::cout<<total_size<<std::endl;
        std::cout<<res_h<<std::endl;
        std::cout<<res_v<<std::endl;
        std::cout<<nb_color<<std::endl;
        std::cout<<nb_important_color<<std::endl;*/



////////////////
        std::cout<<"info"<<std::endl;
        for(int i=0;i<10;++i)
        {
            std::cout<<(*(unsigned*)&info[4*i])<<std::endl;
        }







/////////////////////


        Matrix<Color> img(h, w);
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
                else if(depth == 8)
                {
                    myfile.get(g);
                    if(offset==54)
                        img(i-1, j) = Color(g, g, g);
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
                                img(i-1, j+k) = temp?White:Black;
                            t/=2;
                        }

                        /*img(i, j) = Color(g&0x04, g&0x04, g&0x04);
                        img(i, j+1) = Color(g&0x04, g&0x04, g&0x04);
                        img(i, j+2) = Color(g&0x04, g&0x04, g&0x04);
                        img(i, j+3) = Color(g&0x04, g&0x04, g&0x04);*/
                    }
                }
            }
            for(unsigned j=0;j<(4-(w*3)%4)%4;++j)
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
    std::ofstream myfile(filename, std::ios::out|std::ios::trunc|std::ios::binary);
    if(myfile)
    {
        std::size_t w = img.colNb();
        std::size_t h = img.rowNb();
        unsigned pad_size = (4-(w*3)%4)%4;
        std::size_t filesize = 40 + 3*w*h + w*h/8;
        char bmppad[3] = {0, 0, 0};

        char header[6] = {0,0,1,0,1,0};
        myfile.write(header, 6);

///////////
        /*unsigned w = image_entry[0];
        unsigned h = image_entry[1];
        unsigned nb_color = image_entry[2];
        reserved = image_entry[3];
        unsigned color_plane = image_entry[4]+256*image_entry[5];
        unsigned bpp = image_entry[6]+256*image_entry[7];
        unsigned size = *(unsigned*)&image_entry[8];
        unsigned offset = *(unsigned*)&image_entry[12];*/
/////////

        char image_entry[16] = {w, h, 0, 0, 0, 0, 24, 0,
                                filesize, filesize>>8, filesize>>16, filesize>>24,
                                 22, 0, 0, 0};

        std::cout<<"image_entry"<<std::endl;
        for(int i=0;i<4;++i)
        {
            std::cout<<(*(unsigned*)&image_entry[4*i])<<std::endl;
        }

        myfile.write(image_entry, 16);
        filesize = 3*w*h;
        char info[40] = {40, 0, 0, 0,
                        w, w>>8, w>>16, w>>24,
                        h*2, (h*2)>>8, (h*2)>>16, (h*2)>>24,
                        1, 0, 24, 0,
                        0, 0, 0, 0,
                        filesize, filesize>>8, filesize>>16, filesize>>24};
        unsigned char* info_u = (unsigned char*)info;
        std::cout<<"info"<<std::endl;
        for(int i=0;i<10;++i)
        {
            std::cout<<(*(unsigned*)&info_u[4*i])<<std::endl;
        }
        myfile.write(info, 40);

        for(std::size_t i=h;i>0;--i)
        {
            for(std::size_t j=0;j<w;++j)
            {
                const Color& c = img(i-1, j);
                myfile.put(c.blue());
                myfile.put(c.green());
                myfile.put(c.red());
            }
            myfile.write(bmppad, pad_size);
        }
        for(std::size_t i=0;i<w*h/8;++i)
            myfile.put(0);
        myfile.close();
    }
    else
        std::cout<<"unable to write in the file : "<<filename<<std::endl;
}
