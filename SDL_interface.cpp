#include "SDL_interface.hpp"
#include "ImageProcessing/imgconverter.hpp"

SDL_Window* init_SDL()
{
    SDL_Init(SDL_INIT_VIDEO);
    return SDL_CreateWindow("hello", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 800, 600, SDL_WINDOW_SHOWN);
}

SDL_Renderer* get_SDL_renderer(SDL_Window* screen)
{
    return SDL_CreateRenderer(screen, -1, SDL_RENDERER_PRESENTVSYNC);
}

void pause()
{
    bool done = false;
    SDL_Event event;
    while(!done)
    {
        while(SDL_PollEvent(&event))
        {
            if(event.type==SDL_QUIT)
                done = true;
        }
    }
}

void show_matrix(SDL_Renderer* renderer, const Matrix<Color>& img)
{
    SDL_RenderClear(renderer);
    SDL_Texture* tex = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ABGR8888, SDL_TEXTUREACCESS_STREAMING, img.colNb(), img.rowNb());
    Color* pColor = NULL;

    pColor = (Color*) malloc(sizeof(Color)*img.size());
    std::copy(img.cbegin(), img.cend(), pColor);

    SDL_UpdateTexture(tex, NULL, pColor, img.colNb()*sizeof(Color));
    SDL_RenderCopy(renderer, tex, NULL, NULL);
    SDL_RenderPresent(renderer);
    SDL_DestroyTexture(tex);
    free(pColor);
    pause();
}



void show_matrix(SDL_Renderer* renderer, const Matrix<unsigned char>& img)
{
    show_matrix(renderer, gray2colorimage(img));
}

void show_matrix(SDL_Renderer* renderer, const Matrix<bool>& img)
{
    show_matrix(renderer, bw2colorimage(img));
}

void show_matrix(SDL_Renderer* renderer, const Matrix<std::size_t>& img)
{
    show_matrix(renderer, array2colorimage(img));
}
