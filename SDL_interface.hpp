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

