#include <fstream>
#include <sstream>
#include <string>
#include "csvIO.hpp"


template <class T> Matrix<T> read_csv(std::string filename, char separator)
{
    // check if the separator is valid
    if(separator!=',' && separator!=';' && separator!='\t' && separator!=' ' && separator!=':')
    {
        std::cout<<separator<<" is not a valid separator"<<std::endl;
        return zeros<T>(1);
    }

    std::ifstream myfile(filename, std::ios::in);
    if(myfile)
    {
        // read the first line
        std::string content;
        if(!getline(myfile, content))
        {
            std::cout<<"the file is empty"<<std::endl;
            return zeros<T>(1);
        }
        std::istringstream temp(content);

        // find the number of column
        std::size_t found = content.find(separator);
        int sizey = 1;
        while(found!=std::string::npos)
        {
            content = content.substr(found+1);
            found = content.find(separator);
            ++sizey;
        }
        int i = 0;
        Matrix<T> M(1,sizey);
        for(int j=0;j<sizey;++j)
        {
            char sep;
            temp.get(sep);
            if(sep!=separator)
            {
                temp.putback(sep);
                temp>>M(i, j);
                temp.get(sep);
            }
        }

        // read the other lines
        while(true)
        {
            if(getline(myfile, content))
            {
                M.newRow();
                ++i;
                std::istringstream temp(content);
                for(int j=0;j<sizey;++j)
                {
                    char sep;
                    temp.get(sep);
                    if(sep!=separator)
                    {
                        temp.putback(sep);
                        temp>>M(i, j);
                        temp.get(sep);
                    }
                }
            }
            else
            {
                myfile.close();
                return M;
            }
        }
    }
    std::cout<<"unable to read the file : "<<filename<<std::endl;
    return zeros<T>(1);
}

template <class T> void save_csv(std::string filename, const Matrix<T>& M, char separator)
{
    // check if the separator is valid
    if(separator!=',' && separator!=';' && separator!='\t' && separator!=' ' && separator!=':')
    {
        std::cout<<separator<<" is not a valid separator"<<std::endl;
        return;
    }

    std::ofstream myfile(filename, std::ios::out | std::ios::trunc);
    if(myfile)
    {
        for(int i=0;i<M.rowNb();++i)
        {
            for(int j=0;j<M.colNb()-1;++j)
                myfile<<M(i, j)<<separator;
            myfile<<M(i, M.colNb()-1)<<'\n';
        }
        myfile.close();
    }
    else
        std::cout<<"unable to write in the file : "<<filename<<std::endl;
}
