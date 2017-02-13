#include "gifIO.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <queue>

Matrix<Color> read_gif(std::string filename)
{
    std::ifstream myfile(filename, std::ios::in|std::ios::binary);
    if(myfile)
    {
        //HEADER
        char signature[6];
        myfile.read(signature, 6);

        //SCREEN DESCRIPTOR
        unsigned char screen_descriptor[7];
        myfile.read(reinterpret_cast<char*>(screen_descriptor), 7);
        int W = screen_descriptor[0]+256*screen_descriptor[1];
        int H = screen_descriptor[2]+256*screen_descriptor[3];
        bool M = screen_descriptor[4]/128==1;
        bool Zero = (screen_descriptor[4]%16)>8;
        unsigned char Cr = (screen_descriptor[4]/16)%8+1;
        unsigned char Pixel = (screen_descriptor[4]%8)+1;
        std::cout<<signature<<std::endl;
        std::cout<<W<<" "<<H<<std::endl;
        std::cout<<(int)M<<" "<<(int)Zero<<std::endl;
        std::cout<<(int)Cr<<" "<<(int)Pixel<<std::endl;

        // GLOBAL COLOR MAP
        std::size_t n = M?pow(2, Pixel)*3:0;
        std::cout<<n<<std::endl;
        unsigned char global_color_map[n];
        myfile.read(reinterpret_cast<char*>(global_color_map), n);

        while(!myfile.eof())
        {

            unsigned char c1, c2;
            char c;
            myfile.get(c);
            c1=(unsigned char)c;
            myfile.get(reinterpret_cast<char>(c));
            c2=(unsigned char)c;

            if(c1==0x21 && c2==0xFF)// application extension
            {
                std::cout<<"application extension"<<std::endl;
                myfile.get(c);
                while(c!=0x00)
                {
                    char block[(unsigned char) c];
                    myfile.read(block, (unsigned char) c);
                    myfile.get(c);
                }
            }
            else if(c1==0x21 && c2==0xFE)// comment extension
            {
                std::cout<<"comment extension"<<std::endl;
                myfile.get(c);
                while(c!=0)
                {
                    char block[(unsigned char) c];
                    myfile.read(block, (unsigned char) c);
                    myfile.get(c);
                }
            }
            else
            {
                unsigned char extension_introducer[6];
                if(c1==0x21 && c2==0xF9)// graphic control extension
                {
                    std::cout<<"graphic control extension"<<std::endl;
                    myfile.read(reinterpret_cast<char*>(extension_introducer), 6);
                }
                else
                {
                    std::cout<<"no graphic control extension"<<std::endl;
                    myfile.putback(c2);
                    myfile.putback(c1);
                }

                myfile.get(c);
                c1=(unsigned char)c;
                myfile.get(reinterpret_cast<char>(c));
                c2=(unsigned char)c;

                if(c1==0x21 && c2==0x01)// plain text extension
                {
                    std::cout<<"plain text extension"<<std::endl;
                    myfile.get(c);
                    while(c!=0)
                    {
                        char block[(unsigned char) c];
                        myfile.read(block, (unsigned char) c);
                        myfile.get(c);
                    }
                }
                else // image data
                {

                    myfile.putback(c2);
                    myfile.putback(c1);

                    // image descriptor
                    unsigned char image_descriptor[10];
                    myfile.read(reinterpret_cast<char*>(image_descriptor), 10);
                    char separator = image_descriptor[0];
                    int img_left = image_descriptor[1]+256*image_descriptor[2];
                    int img_top = image_descriptor[3]+256*image_descriptor[4];
                    int img_width = image_descriptor[5]+256*image_descriptor[6];
                    int img_height = image_descriptor[7]+256*image_descriptor[8];
                    bool M_img = image_descriptor[9]/128==1;
                    bool I_img = (image_descriptor[9]/64)%2;
                    int Pixel_img = (image_descriptor[9]%8)+1;
                    std::cout<<separator<<" "<<(unsigned char)(0x2C)<<std::endl;
                    std::cout<<img_width<<" "<<img_height<<std::endl;
                    std::cout<<(int)M_img<<" "<<(int)I_img<<std::endl;
                    std::cout<<(int)Cr<<" "<<(int)Pixel_img<<std::endl;

                    Matrix<Color> img(img_height, img_width);
                    auto it = img.begin();

                    // local color table
                    int n_img = M_img?pow(2, Pixel_img)*3:0;
                    std::cout<<n_img<<std::endl;
                    unsigned char local_color_map[n_img];
                    myfile.read(reinterpret_cast<char*>(local_color_map), n_img);

                    // current color map
                    unsigned char* color_map = M_img?local_color_map:global_color_map;
                    int n_color_map = M_img?n_img:n;

                    // LZW_minimum_size
                    unsigned char LZW_minimum_size;
                    char temp;
                    myfile.get(temp);
                    LZW_minimum_size = (unsigned char)temp;
                    std::cout<<"LZW min size "<<(int)LZW_minimum_size<<std::endl;

                    // nb bytes subblock
                    char nb_bytes_following;
                    unsigned char nb_bytes_following_u;
                    myfile.get(nb_bytes_following);
                    nb_bytes_following_u = nb_bytes_following;
                    unsigned code_size = LZW_minimum_size+1;
                    int cpt_subblock = 0;

                    // init dictionnary
                    std::string Dictionnary[4096];
                    bool Existence_in_dic[4096];
                    std::size_t limit_dic = n_color_map/3+2;
                    std::size_t cur_idx_in_Dic = limit_dic;
                    std::cout<<"limit dict "<<limit_dic<<std::endl;
                    std::fill(Dictionnary, Dictionnary+4096, "");
                    std::fill(Existence_in_dic, Existence_in_dic+4096, false);
                    for(std::size_t i=0;i<limit_dic;++i)
                    {
                        Dictionnary[i] = char(i);
                        Existence_in_dic[i] = true;
                    }
                    cur_idx_in_Dic = limit_dic;
                    std::string output = "";

                    while(nb_bytes_following_u!=0)
                    {

                        // read next subblock
                        std::cout<<"new subblock "<<cpt_subblock<<std::endl;
                        ++cpt_subblock;
                        //std::cout<<(int)nb_bytes_following_u<<std::endl;
                        unsigned char subblock[nb_bytes_following_u];
                        myfile.read(reinterpret_cast<char*>(subblock), nb_bytes_following_u);
                        std::queue<bool> cur_bitset;
                        for(unsigned char i=0;i<nb_bytes_following_u;++i)
                        {
                            for(unsigned char j=0;j<8;++j)
                            {
                                cur_bitset.push(subblock[i]%2);
                                subblock[i]/=2;
                            }
                        }
                        //std::cout<<cur_bitset.size()<<std::endl;
                        //std::cout<<code_size<<std::endl;


                        // read first index
                        unsigned code = 0;
                        unsigned previous_code = 0;
                        for(unsigned i=0; i<code_size;++i)
                        {
                            previous_code+=cur_bitset.front()?pow(2, i):0;
                            cur_bitset.pop();
                        }
                        //if(previous_code!=limit_dic-2)
                            //output+=Dictionnary[code];

                        for(unsigned i=0; i<code_size;++i)
                        {
                            code+=cur_bitset.front()?pow(2, i):0;
                            cur_bitset.pop();
                        }
                        previous_code = code;
                        output+=Dictionnary[code];
                        //std::cout<<previous_code<<std::endl;

                        while(cur_bitset.size()>=code_size)
                        {
                            code = 0;
                            for(unsigned i=0; i<code_size;++i)
                            {
                                code+=cur_bitset.front()?pow(2, i):0;
                                cur_bitset.pop();
                            }

                            if(code==limit_dic-2)// clear code
                            {
                                std::cout<<"clear code "<<limit_dic-2<<std::endl;
                                for(std::size_t i=limit_dic;i<4096;++i)
                                {
                                    Existence_in_dic[i] = false;
                                    Dictionnary[i] = "";
                                }
                                cur_idx_in_Dic = limit_dic;
                                code_size = LZW_minimum_size+1;
                                previous_code = code;
                                continue;
                            }
                            else if(code==limit_dic-1)// end code
                            {
                                std::for_each(output.cbegin(), output.cend(),
                                            [&it, &color_map](const unsigned char& c)
                                            {*it = Color(color_map[3*c],
                                                        color_map[3*c+1],
                                                        color_map[3*c+2]);
                                            ++it;});
                                std::cout<<"end code "<<limit_dic-1<<std::endl;
                                std::cout<<"nb de pixels ecrits : ";
                                std::cout<<output.size()<<std::endl;
                                return img;
                            }
                            else if(Existence_in_dic[code])
                            {
                                Dictionnary[cur_idx_in_Dic] = Dictionnary[previous_code]+Dictionnary[code][0];
                                output+=Dictionnary[code];
                                Existence_in_dic[cur_idx_in_Dic] = true;
                                ++cur_idx_in_Dic;

                                if(cur_idx_in_Dic==pow(2, code_size))
                                {
                                    ++code_size;
                                    std::cout<<"new code size "<<code_size<<std::endl;
                                }

                                previous_code = code;
                            }
                            else
                            {
                                Dictionnary[cur_idx_in_Dic] = Dictionnary[previous_code]+Dictionnary[previous_code][0];
                                output+=Dictionnary[cur_idx_in_Dic];
                                Existence_in_dic[cur_idx_in_Dic] = true;
                                ++cur_idx_in_Dic;

                                if(cur_idx_in_Dic==pow(2, code_size))
                                {
                                    ++code_size;
                                    std::cout<<"new code size "<<code_size<<std::endl;
                                }

                                previous_code = code;
                            }
                        }

                        myfile.get(nb_bytes_following);
                        nb_bytes_following_u = nb_bytes_following;
                    }
                    return img;
                }
            }
        }
        myfile.close();
    }
    std::cout<<"Unable to read the file : "<<filename<<std::endl;
    return Matrix<Color>();
}

void save_gif(std::string filename, const Matrix<Color>& img)
{

}
