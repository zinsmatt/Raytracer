#include <iostream>
#include <SDL/SDL.h>
#include <ctime>
#include <cstdlib>
#include <random>
#include <cmath>
#define WIDTH 640
#define HEIGHT 480

using namespace std;

int main()
{
	std::srand(std::time(0));

	cout << "Hello World!" << endl;
	SDL_Surface *screen;
	if( SDL_Init(SDL_INIT_VIDEO) == -1)
	{
		return EXIT_FAILURE;
	}

	atexit(SDL_Quit);
	screen = SDL_SetVideoMode(WIDTH,HEIGHT,32,SDL_HWSURFACE);

	if( screen == NULL)
	{
		return EXIT_FAILURE;
	}


	std::default_random_engine generator;
	std::normal_distribution<double> distribution(100.0,20.0);


	char *p = (char*)screen->pixels;
	int bpp = screen->format->BytesPerPixel;
	std::cout << "bpp = " << bpp << std::endl;
	for(int i=0;i<HEIGHT*WIDTH*bpp;i++)
	{
		//p[i] = std::rand()%256;
		p[i] = (int)floor(distribution(generator));
	}

	SDL_Flip(screen);

	SDL_Delay(3000);

	SDL_FreeSurface(screen);
	return 0;
}

