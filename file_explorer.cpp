#include <dirent.h>
#include <iostream>
#include "file_explorer.hpp"


bool is_file(std::string elem)
{
    return !is_folder(elem);
}

bool is_folder(std::string elem)
{
    DIR* cur_dir = opendir(elem.c_str());
    return cur_dir!=NULL;
}

std::vector<std::string> get_files(std::string directory)
{
    DIR* cur_dir = opendir(directory.c_str());
    std::vector<std::string > file_list;
    file_list.clear();
    if(cur_dir==NULL)
    {
        std::cout<<"Unable to open the folder "+directory<<std::endl;
        return file_list;
    }
    else
    {
        struct dirent* cur_file = NULL;
        while((cur_file = readdir(cur_dir))!=NULL)
        {
            if(is_file(directory+"/"+cur_file->d_name))
            {
                file_list.push_back(cur_file->d_name);
            }
        }
        closedir(cur_dir);
        return file_list;
    }
}

std::vector<std::string> get_folders(std::string directory)
{
    DIR* cur_dir = opendir(directory.c_str());
    std::vector<std::string > file_list;
    file_list.clear();
    if(cur_dir==NULL)
    {
        std::cout<<"Unable to open the folder "+directory<<std::endl;
        return file_list;
    }
    else
    {
        struct dirent* cur_file = NULL;
        while((cur_file = readdir(cur_dir))!=NULL)
        {
            if(is_folder(directory+"/"+cur_file->d_name))
            {
                file_list.push_back(cur_file->d_name);
            }
        }
        closedir(cur_dir);
        return file_list;
    }
}

std::vector<std::string> get_files_and_folders(std::string directory)
{
    DIR* cur_dir = opendir(directory.c_str());
    std::vector<std::string > file_list;
    file_list.clear();
    if(cur_dir==NULL)
    {
        std::cout<<"Unable to open the folder "+directory<<std::endl;
        return file_list;
    }
    else
    {
        struct dirent* cur_file = NULL;
        while((cur_file = readdir(cur_dir))!=NULL)
        {
            file_list.push_back(cur_file->d_name);
        }
        closedir(cur_dir);
        return file_list;
    }
}

std::vector<std::string> get_files_recursively(std::string directory, std::string extension)
{
    if(directory[directory.size()-1]=='/')
    directory.erase(directory.end()-1);
    DIR* cur_dir = opendir(directory.c_str());
    std::vector<std::string > file_list;
    file_list.clear();
    if(cur_dir==NULL)
    {
        std::cout<<"Unable to open the folder "+directory<<std::endl;
        return file_list;
    }
    else
    {
        struct dirent* cur_file = NULL;
        while((cur_file = readdir(cur_dir))!=NULL)
        {
            std::string cur_filname = cur_file->d_name;
            if(is_file(directory+"/"+cur_filname) && cur_filname.size()>extension.size())
            {

                if(!extension.compare(cur_filname.substr(cur_filname.size()-extension.size(),extension.size())))
                    file_list.push_back(directory+"/"+cur_file->d_name);
            }
            else if(is_folder(directory+"/"+cur_filname) && cur_filname.compare("..") && cur_filname.compare("."))
            {
                std::vector<std::string > subfolders_list = get_files_recursively(directory+"/"+cur_filname,extension);
                for(std::size_t i=0;i<subfolders_list.size();++i)
                {
                    file_list.push_back(subfolders_list[i]);
                }
            }
        }
        closedir(cur_dir);
        return file_list;
    }
}
