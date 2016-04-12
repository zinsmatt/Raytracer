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
	tabSpheres[0].mat.ka = { 0.3,0.3,0.0};
	tabSpheres[0].mat.kd = { 0.7,0.7,0.0};
	tabSpheres[0].mat.ks = { 0.7,0.7, 0.7};
	tabSpheres[0].mat.shininess = 38;
	nbSpheres++;

	tabSpheres[1].xc = -0.2;
	tabSpheres[1].yc = 0.3;
	tabSpheres[1].zc = -1;
	tabSpheres[1].r = 0.3;
	tabSpheres[1].mat.ka = { 0.0,0.3,0.0};
	tabSpheres[1].mat.kd = { 0.0,0.7,0.0};
	tabSpheres[1].mat.ks = { 1.0,1.0,1.0};
	tabSpheres[1].mat.shininess = 52;
	nbSpheres++;

	tabSpheres[2].xc = 0.4;
	tabSpheres[2].yc = 0.3;
	tabSpheres[2].zc = -1.1;
	tabSpheres[2].r = 0.05;
	tabSpheres[2].mat.ka = { 0.0,0.0,0.4};
	tabSpheres[2].mat.kd = { 0.0,0.0,0.7};
	tabSpheres[2].mat.ks = { 1.0,1.0,1.0};
	tabSpheres[2].mat.shininess = 128;
	nbSpheres++;

	tabSources[0].Ia = {0.1,0.1,0.1};
	tabSources[0].Ip = {0.95,0.95,0.95};
	tabSources[0].pos = {3,3,2};
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
				pt = {start.x+dir.x*t,start.y+dir.y*t,start.z+dir.z*t};
				Point center = {sphere->xc,sphere->yc,sphere->zc};
				Vector normal;
				normal.fromPoints(center,pt).normalize();
				face.mat = sphere->mat;
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
					pt = {start.x+dir.x*t,start.y+dir.y*t,start.z+dir.z*t};
					Point center = {sphere->xc,sphere->yc,sphere->zc};
					Vector normal;
					normal.fromPoints(center,pt).normalize();
					face.mat = sphere->mat;
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
			Vector Lj, reflected;
			//ambient only with the first source
			temp = source1->Ia * face.mat.ka;
			// diffuse with all the sources
			Vector vision;
			vision.fromPoints(pointI,obs).normalize();
			for(int srcIter = 0; srcIter<nbSources; ++srcIter)
			{
				Lj.fromPoints(pointI,tabSources[srcIter].pos).normalize();
				double d = dist(tabSources[srcIter].pos,pointI);
				double fAtt = 1/(0.01*d*d+0.004*d+0.003);
				if(fAtt>1) fAtt = 1;
				temp +=  (tabSources[srcIter].Ip * fAtt * face.mat.kd) * (face.normal * Lj);
				if(face.normal * Lj >0)
				{
					reflected = (face.normal * 2 *(face.normal * Lj)) - Lj;
					temp += (tabSources[srcIter].Ip * fAtt * face.mat.ks) * pow(reflected*vision,(double)face.mat.shininess);
				}
			}
			it = temp;
			if(it.r>1.0) it.r = 1.0;
			if(it.g>1.0) it.g = 1.0;
			if(it.b>1.0) it.b = 1.0;
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
	SDL_LockSurface(screen);
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

			direction.normalize();
			raytrace2(obs,direction,2,it);

			uint8_t r = (uint8_t)(floor(it.r*255));
			uint8_t g = (uint8_t)(floor(it.g*255));
			uint8_t b = (uint8_t)(floor(it.b*255));
			setPixel(screen,i,j,r,g,b);
		}
	}
	SDL_UnlockSurface(screen);
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

