#include <iostream>
#include <SDL/SDL.h>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <random>
#include <cmath>
#include "raytracer.h"
#define WIDTH 640
#define HEIGHT 480
#define MAX_SPHERES 8
#define MAX_SOURCES 1
#define PI 3.14159265359
using namespace std;


Sphere tabSpheres[MAX_SPHERES];
int nbSpheres = 0;

int nbSources = 0;
Source tabSources[MAX_SOURCES];
Source *source1 = &tabSources[0];

Point obs = {0,0,0};
double znear = 1;

void init()
{
	tabSpheres[0].xc = 0;
	tabSpheres[0].yc = 0;
	tabSpheres[0].zc = -2;
	tabSpheres[0].r = 1;
	tabSpheres[0].ka = 0.7;
	tabSpheres[0].kd = 0.7;
	nbSpheres++;

	tabSources[0].Ia = {1.0,1.0,0.0};
	tabSources[0].Ip = {1.0,1.0,0.0};
	tabSources[0].pos = {0,5,-2};
	nbSources++;
}

/*
void raytrace(Point start, Vector dir, int nb, Intensity it)
{
	if(nb == 0)
	{
		it = {0,0,0};
	}else if(! intersect(start, dir, pointI, face))
	{
		it = {0,0,0};
	}else
	{
		if(face.ks == 0)
		{
			Is = 0;
		}else
		{
			computeReflected(dir,face.normal,r);
			raytrace(pointI,r,nb-1,Is);
		}

		if(face.kt == 0)
		{
			It = 0;
		}else
		{
			computeTransmitted(dir,face.normal,face.kro,p);
			raytrace(pointI,p,nb-1,It);
		}

		I = Ia*face.ka;

		for(int iter=0; iter<nbSources; ++iter)
		{
			if(!shadow(pointI,face.normal,dir,j)){
				//I = I + Idj*(face.kd) ....
			}
		}
		I += face.ks * Is + face.kt * It;
		if(start != observer)
		{
			computeDistance(start,pointI,d);
			I /= d;
		}
		intensite = I;

	}
}






void draw()
{
	Intensity it;
	for(int i=0;i<HEIGHT; ++i)
	{
		for(int j=0;j<WIDTH; ++j)
		{
			Vector direction;
			Point start;
			raytrace(start,direction,2,it);
		}
	}
}*/



bool findNextIntersection(Point start, Vector dir, Point& pt, Face& face)
{
	bool intersection = false;
	double tmin = HUGE_VALF;
	for(int sphereIter = 0; sphereIter<nbSpheres; ++sphereIter)
	{
		Sphere *sphere = &tabSpheres[sphereIter];
		double A = dir.x*dir.x + dir.y*dir.y + dir.z*dir.z;
		double B = 2*(dir.x*(start.x-sphere->xc)+dir.y*(start.y-sphere->yc)+dir.z*(start.z-sphere->zc));
		double C = (start.x-sphere->xc)*(start.x-sphere->xc)+
		           (start.y-sphere->yc)*(start.y-sphere->yc)+
		           (start.z-sphere->zc)*(start.z-sphere->zc)-
		           sphere->r*sphere->r;
		double delta = B*B-4*A*C;
		if(delta == 0)
		{
			double t = -B / (2*A);
			if(t>0 && t<tmin)
			{
				intersection = true;
				tmin = t;
				pt = {pt.x+dir.x*t,pt.y+dir.y*t,pt.z+dir.z*t};
				Point center = {sphere->xc,sphere->yc,sphere->zc};
				Vector normal = vectNormFromPoints(center,pt);
				face.ka = sphere->ka;
				face.kd = sphere->kd;
				face.normal = normal;
			}
		}else if(delta>0)
		{
			double t1 = (-B-sqrt(delta))/(2*A);
			double t2 = (-B+sqrt(delta))/(2*A);
			double t;
			if(t1>0 || t2>0)
			{
				if(t1<0)	t = t2;
				else if(t2<0)	t = t1;
				else	t = (t1<t2)?t1:t2;

				if(t<tmin)
				{
					intersection = true;
					tmin = t;
					pt = {pt.x+dir.x*t,pt.y+dir.y*t,pt.z+dir.z*t};
					Point center = {sphere->xc,sphere->yc,sphere->zc};
					Vector normal = vectNormFromPoints(center,pt);
					face.ka = sphere->ka;
					face.kd = sphere->kd;
					face.normal = normal;
				}
			}
		}
	}
	return intersection;
}






void raytrace2(Point start, Vector dir, int nb, Intensity& it)
{
	if(nb == 0)
	{
		it = {0,0,0};
	}else
	{
		Point pointI;
		Face face;
		if(!findNextIntersection(start,dir,pointI,face))
		{
			it = {0,0,0};
		}else
		{
			Intensity temp;
			temp = mult(source1->Ia,face.ka);	//ambient only with the first source
			// diffuse with all the sources
			for(int srcIter = 0; srcIter<nbSources; ++srcIter)
			{
				Vector Lj = vectNormFromPoints(pointI,tabSources[srcIter].pos);
				double fAtt = 1/dist(tabSources[srcIter].pos,pointI);
				double coef = face.kd*prodScalPos(face.normal,Lj);
				add(temp,mult(tabSources[srcIter].Ip,coef));
			}
			it = temp;
		}
	}
}

void setPixel(SDL_Surface *screen, int i, int j, uint8_t r, uint8_t g, uint8_t b)
{
	int index = (i*WIDTH+j);
	uint32_t *pixel = (uint32_t*) screen->pixels;
	pixel += index;
	uint32_t mask = g;
	mask <<= 8;
	uint32_t value = r;
	value <<= 16;
	value |= mask;
	mask = b;
	value |= mask;
	*pixel = value;
}


void draw2(SDL_Surface *screen)
{
	char *p = (char*)screen->pixels;


	double fovy = 90;
	double ratio = (double)HEIGHT/WIDTH;
	double right = znear * tan((fovy*PI/180)/2);
	double top = ratio * right;
	double left = -right;
	double bottom = -top;
	double stepH = (right-left) / WIDTH;
	double stepV = (top-bottom) / HEIGHT;


	Intensity it;
	Vector direction;
	for(int i=0;i<HEIGHT; ++i)
	{
		for(int j=0;j<WIDTH; ++j)
		{
			direction.x = left + j * stepH;
			direction.y = top - i * stepV;
			direction.z = -znear;

			normalize(direction);
			raytrace2(obs,direction,2,it);

			uint8_t r = (uint8_t)(floor(it.r*255));
			uint8_t g = (uint8_t)(floor(it.g*255));
			uint8_t b = (uint8_t)(floor(it.b*255));
			setPixel(screen,i,j,r,g,b);
		}
	}
}






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
	//SDL_SetAlpha(screen, SDL_SRCALPHA, 255);

	SDL_PixelFormat *f;
	f = screen->format;
	cout << "r mask " << f->Rmask << "\n";

	if( screen == NULL)
	{
		return EXIT_FAILURE;
	}


	std::default_random_engine generator;
	std::normal_distribution<double> distribution(100.0,20.0);


	/*
	char *p = (char*)screen->pixels;
	int bpp = screen->format->BytesPerPixel;
	std::cout << "bpp = " << bpp << std::endl;
	for(int i=0;i<HEIGHT*WIDTH*bpp;i++)
	{
		//p[i] = std::rand()%256;
		p[i] = (int)floor(distribution(generator));
	}*/

	init();
	draw2(screen);
	std::cout << "Drawing end\n";

	SDL_Flip(screen);

	bool quit = false;
	SDL_Event event;
	while(SDL_WaitEvent(&event) >=0 && !quit)
	{
		switch(event.type)
		{
			case SDL_QUIT:
				quit = true;
				break;
		}
		if(quit)
			break;
	}

	//SDL_Delay(3000);


	SDL_FreeSurface(screen);
	return 0;
}
