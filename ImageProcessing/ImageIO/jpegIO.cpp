/*#include <stdio.h>
#include <jpeglib.h>
#include <jerror.h>

#include "jpegIO.hpp"


Matrix<Color> read_jpeg(std::string filename)
{
    struct jpeg_decompress_struct info;
    struct jpeg_error_mgr err;
    FILE* infile;
    int w, h, channels;

    infile = fopen(filename.c_str(), "rb");
    if(infile == NULL)
    {
        std::cout<<"Unable to read the file : "<<filename<<std::endl;
        return Matrix<Color>();
    }

    info.err = jpeg_std_error(& err);
    jpeg_create_decompress(& info);
    jpeg_stdio_src(&info, infile);
    jpeg_read_header(&info, TRUE);
    jpeg_start_decompress(&info);

    w = info.output_width;
    h = info.output_height;
    channels = info.num_components;

    Matrix<Color> img(h, w);
    unsigned char* rowptr[1];
    unsigned char * jdata = (unsigned char *)malloc(channels*w*h);
    while(info.output_scanline < info.output_height)
    {
        rowptr[0] = jdata+channels* info.output_width * info.output_scanline;
        jpeg_read_scanlines(&info, rowptr, 1);
    }
    jpeg_finish_decompress(&info);
    jpeg_destroy_decompress(&info);

    if(channels==1)
        std::transform(jdata, jdata+w*h, img.begin(), [](const unsigned char& c){return Color(c, c, c);});
    else if(channels==3)
    {
        unsigned char* pix_ptr = jdata;
        std::for_each(img.begin(), img.end(), [&pix_ptr](Color& c){c = Color(*pix_ptr, *(pix_ptr+1), *(pix_ptr+2)); pix_ptr+=3;});
    }

    fclose(infile);
    free(jdata);
    return img;
}

void save_jpeg(std::string filename, const Matrix<Color>& img)
{

    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;
    FILE * outfile;
    JSAMPROW row_pointer[1];
    int row_stride;

    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);
    if ((outfile = fopen(filename.c_str(), "wb")) == NULL)
    {
        std::cout<<"unable to write in the file : "<<filename<<std::endl;
        return;
    }

    jpeg_stdio_dest(&cinfo, outfile);
    cinfo.image_width = img.colNb();
    cinfo.image_height = img.rowNb();
    cinfo.input_components = 3;
    cinfo.in_color_space = JCS_RGB;

    int quality = 90;
    jpeg_set_defaults(&cinfo);
    jpeg_set_quality(&cinfo, quality, TRUE);

    jpeg_start_compress(&cinfo, TRUE);


    row_stride = img.colNb() * 3;
    unsigned char* image_buffer = (unsigned char*)malloc(3*img.size());
    int cpt = 0;
    for(std::size_t i=0;i<img.rowNb();++i)
    {
        for(std::size_t j=0;j<img.colNb();++j)
        {
            Color cur_color = img(i, j);
            image_buffer[cpt++] = cur_color.red();
            image_buffer[cpt++] = cur_color.green();
            image_buffer[cpt++] = cur_color.blue();

        }
    }


    while (cinfo.next_scanline < cinfo.image_height)
    {
        row_pointer[0] = & image_buffer[cinfo.next_scanline * row_stride];
        jpeg_write_scanlines(&cinfo, row_pointer, 1);
    }


    jpeg_finish_compress(&cinfo);
    jpeg_destroy_compress(&cinfo);
    fclose(outfile);
}
*/
