/*!
 * \author Alexandre Krebs
 * \file SDL_interface.hpp
 * \brief functions to show images
 */


#pragma once

#ifndef SDL_INTERFACE_HPP
#define SDL_INTERFACE_HPP


#define SDL_MAIN_HANDLED
#include <SDL2/SDL.h>
#include "matrix.hpp"
#include "ImageProcessing/color.hpp"


SDL_Window* init_SDL();

SDL_Renderer* get_SDL_renderer(SDL_Window* screen);

void pause();

void show_matrix(SDL_Renderer* renderer, const Matrix<Color>& img);
void show_matrix(SDL_Renderer* renderer, const Matrix<unsigned char>& img);
void show_matrix(SDL_Renderer* renderer, const Matrix<bool>& img);
void show_matrix(SDL_Renderer* renderer, const Matrix<std::size_t>& img);


#include "SDL_interface.hpp"


#endif // SDL_INTERFACE_HPP
