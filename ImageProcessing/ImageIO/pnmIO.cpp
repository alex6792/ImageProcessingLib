#include "pnmIO.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>


Matrix<unsigned char> read_pgm(std::string filename)
{
    std::ifstream myfile(filename, std::ios::in|std::ios::binary);
    if(myfile)
    {
        std::string magicNb;
        getline(myfile, magicNb);
        magicNb = magicNb.substr(0,2);
        if(magicNb.compare("P2") && magicNb.compare("P5"))
        {
            std::cout<<"the file : "<<filename<<" is not a pgm file"<<std::endl;
            Matrix<unsigned char> img = zeros<unsigned char>(1);
            return img;
        }
        std::string comments = "#";
        do
        {
            if(!getline(myfile, comments))
            {
                std::cout<<"error reach the end of the file"<<std::endl;
                Matrix<unsigned char> img = zeros<unsigned char>(1);
                return img;
            }
        }while(comments[0]=='#');
        std::istringstream mystream(comments);
        int W, H;
        mystream>>W;
        std::string rem = (mystream.str().substr(mystream.tellg()));
        if(rem.size()>1)
            mystream>>H;
        else
            myfile>>H;
        int maxvalue;
        myfile>>maxvalue;
        char c = maxvalue;
        while('\n'!=c)
        {
            myfile.get(c);
        }
        Matrix<unsigned char> img(H, W);
        if(!magicNb.compare("P2"))
        {
            int data;
            for(int i=0;i<H;++i)
            {
                for(int j=0;j<W;++j)
                {
                    myfile>>data;
                    img(i, j) = data;
                }
            }
        }
        else
        {
            for(int i=0;i<H;++i)
            {
                for(int j=0;j<W;++j)
                {
                    char c;
                    myfile.get(c);
                    img(i, j) = c;
                }
            }
        }
        myfile.close();
        return img;
    }
    std::cout<<"Unable to read the file : "<<filename<<std::endl;
    return zeros<unsigned char>(1);
}

Matrix<bool> read_pbm(std::string filename)
{
    std::ifstream myfile(filename, std::ios::in|std::ios::binary);
    if(myfile)
    {
        std::string magicNb;
        getline(myfile, magicNb);
        magicNb = magicNb.substr(0,2);
        if(magicNb.compare("P1") && magicNb.compare("P4"))
        {
            std::cout<<"the file : "<<filename<<" is not a pbm file"<<std::endl;
            Matrix<bool> img = zeros<bool>(1);
            return img;
        }
        std::string comments = "#";
        do
        {
            if(!getline(myfile, comments))
            {
                std::cout<<"error reach the end of the file"<<std::endl;
                Matrix<bool> img = zeros<bool>(1);
                return img;
            }
        }while(comments[0]=='#');
        std::istringstream mystream(comments);
        int W, H;
        mystream>>W;
        std::string rem = (mystream.str().substr(mystream.tellg()));
        if(rem.size()>1)
            mystream>>H;
        else
            myfile>>H;
        Matrix<bool> img(H, W);
        if(!magicNb.compare("P1"))
        {
            unsigned char data;
            for(int i=0;i<H;++i)
            {
                for(int j=0;j<W;++j)
                {
                    myfile>>data;
                    img(i, j) = !(data-'0');
                }
            }
        }
        else
        {
            char data;
            for(int i=0;i<H;++i)
            {
                for(int j=0;j<W;j+=8)
                {
                    myfile.get(data);
                    unsigned char t = 128;
                    for(int k=0;k<8;++k)
                    {
                        if(j+k<W)
                            img(i, j+k) = !(bool)((unsigned char)data%(2*t)/t);
                        t/=2;
                    }
                }
            }
        }
        myfile.close();
        return img;
    }
    std::cout<<"Unable to read the file : "<<filename<<std::endl;
    Matrix<bool> img = zeros<bool>(1);
    return img;
}

Matrix<Color> read_ppm(std::string filename)
{
    std::ifstream myfile(filename, std::ios::in|std::ios::binary);
    if(myfile)
    {
        std::string magicNb;
        getline(myfile, magicNb);
        magicNb = magicNb.substr(0,2);
        if(magicNb.compare("P3") && magicNb.compare("P6"))
        {
            std::cout<<"the file : "<<filename<<" is not a ppm file"<<std::endl;
            Matrix<Color> img(1);
            return img;
        }
        std::string comments = "#";
        do
        {
            if(!getline(myfile, comments))
            {
                std::cout<<"error reach the end of the file"<<std::endl;
                Matrix<Color> img(1);
                return img;
            }
        }while(comments[0]=='#');
        std::istringstream mystream(comments);
        int W, H;
        mystream>>W;
        std::string rem = (mystream.str().substr(mystream.tellg()));
        if(rem.size()>1)
            mystream>>H;
        else
            myfile>>H;
        int maxvalue;
        myfile>>maxvalue;
        char c = maxvalue;
        while('\n'!=c)
        {
            myfile.get(c);
        }
        Matrix<Color> img(H, W);
        if(!magicNb.compare("P3"))
        {
            int r,g,b;
            for(int i=0;i<H;++i)
            {
                for(int j=0;j<W;++j)
                {
                    myfile>>r;
                    myfile>>g;
                    myfile>>b;
                    img(i, j) = Color(r,g,b);
                }
            }
        }
        else
        {
            char r,g,b;
            for(int i=0;i<H;++i)
            {
                for(int j=0;j<W;++j)
                {
                    myfile.get(r);
                    myfile.get(g);
                    myfile.get(b);
                    img(i, j) = Color(r,g,b);
                }
            }
        }
        myfile.close();
        return img;
    }
    std::cout<<"Unable to read the file : "<<filename<<std::endl;
    return Matrix<Color>(1);
}

void save_pbm(std::string filename, const Matrix<bool>& img)
{
    std::ofstream myfile(filename, std::ios::out | std::ios::trunc);
    if(myfile)
    {
        std::string separator = " ";
        int W = img.colNb();
        int H = img.rowNb();
        myfile<<"P1"<<std::endl;
        myfile<<W<<separator<<H<<std::endl;
        for(int i=0;i<H;++i)
        {
            for(int j=0;j<W-1;++j)
                myfile<<img(i, j)<<separator;
            myfile<<img(i, W-1)<<std::endl;
        }
        myfile.close();
    }
    else
        std::cout<<"unable to write in the file : "<<filename<<std::endl;
}

void save_pgm(std::string filename, const Matrix<unsigned char>& img)
{
    std::ofstream myfile(filename, std::ios::out | std::ios::trunc);
    if(myfile)
    {
        std::string separator = " ";
        int W = img.colNb();
        int H = img.rowNb();
        myfile<<"P2"<<std::endl;
        myfile<<W<<separator<<H<<std::endl;
        myfile<<"255"<<std::endl;
        for(int i=0;i<H;++i)
        {
            for(int j=0;j<W-1;++j)
                myfile<<(int)img(i, j)<<separator;
            myfile<<(int)img(i, W-1)<<std::endl;
        }
        myfile.close();
    }
    else
        std::cout<<"unable to write in the file : "<<filename<<std::endl;
}

void save_ppm(std::string filename, const Matrix<Color>& img)
{
    std::ofstream myfile(filename, std::ios::out | std::ios::trunc);
    if(myfile)
    {
        std::string separator = " ";
        int W = img.colNb();
        int H = img.rowNb();
        myfile<<"P3"<<std::endl;
        myfile<<W<<separator<<H<<std::endl;
        myfile<<"255"<<std::endl;

        Color my_Color;
        for(int i=0;i<H;++i)
        {
            for(int j=0;j<W-1;++j)
            {
                my_Color = img(i, j);
                myfile<<(int)my_Color.red()<<separator<<(int)my_Color.green()<<separator<<(int)my_Color.blue()<<separator<<separator;
            }
            my_Color = img(i, W-1);
            myfile<<(int)my_Color.red()<<separator<<(int)my_Color.green()<<separator<<(int)my_Color.blue()<<std::endl;
        }
        myfile.close();
    }
    else
        std::cout<<"unable to write in the file : "<<filename<<std::endl;
}
