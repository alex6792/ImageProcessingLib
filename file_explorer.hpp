/*!
 * \author Alexandre Krebs
 * \file file_explorer.hpp
 * \brief file manager
 */


#pragma once

#ifndef FILE_EXPLORER_HPP
#define FILE_EXPLORER_HPP


#include <string>
#include <vector>


bool is_file(std::string);
bool is_folder(std::string);

std::vector<std::string> get_files(std::string);
std::vector<std::string> get_folders(std::string);
std::vector<std::string> get_files_and_folders(std::string);

std::vector<std::string> get_files_recursively(std::string, std::string);


#endif // FILE_EXPLORER_HPP
