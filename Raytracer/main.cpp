#include <iostream>
#include <SDL/SDL.h>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <random>
#include <cmath>
#include <unistd.h>
#include "raytracer.h"
#define WIDTH 640
#define HEIGHT 480
#define MAX_SPHERES 8
#define MAX_SOURCES 4
#define PI 3.14159265359
using namespace std;


Sphere tabSpheres[MAX_SPHERES];
int nbSpheres = 0;

int nbSources = 0;
Source tabSources[MAX_SOURCES];
Source *source1 = &tabSources[0];

Point obs = {0,0,0};
double znear = 1;

bool antialiasing = false;

void init()
{
	tabSpheres[0].xc = 0;
	tabSpheres[0].yc = 0;
	tabSpheres[0].zc = -2;
	tabSpheres[0].r = 1;
	tabSpheres[0].mat.ka = { 0.3,0.3,0.0};
	tabSpheres[0].mat.kd = { 0.7,0.7,0.0};
	tabSpheres[0].mat.ks = { 1.0,1.0,1.0};
	tabSpheres[0].mat.kt = { 0.0,0.0,0.0};
	tabSpheres[0].mat.ro = 1.0;
	tabSpheres[0].mat.shininess = 78;
	nbSpheres++;

	tabSpheres[1].xc = -0.4;
	tabSpheres[1].yc = 0.3;
	tabSpheres[1].zc = -0.8;
	tabSpheres[1].r = 0.1;
	tabSpheres[1].mat.ka = { 0.0,0.3,0.0};
	tabSpheres[1].mat.kd = { 0.0,0.7,0.0};
	tabSpheres[1].mat.ks = { 1.0,1.0,1.0};
	tabSpheres[1].mat.kt = { 0.0,0.0,0.0};
	tabSpheres[1].mat.ro = 1.0;
	tabSpheres[1].mat.shininess = 52;
	nbSpheres++;

	tabSpheres[2].xc = 0.4;
	tabSpheres[2].yc = -0.3;
	tabSpheres[2].zc = -1.1;
	tabSpheres[2].r = 0.05;
	tabSpheres[2].mat.ka = { 0.0,0.0,0.4};
	tabSpheres[2].mat.kd = { 0.0,0.0,0.7};
	tabSpheres[2].mat.ks = { 1.0,1.0,1.0};
	tabSpheres[2].mat.kt = { 0.0,0.0,0.0};
	tabSpheres[2].mat.ro = 1.0;
	tabSpheres[2].mat.shininess = 128;
	nbSpheres++;

	tabSources[0].Ia = {0.1,0.1,0.1};
	tabSources[0].Ip = {0.95,0.95,0.95};
	tabSources[0].pos = {5,5,1};
	nbSources++;

	//tabSources[1].Ia = {0.1,0.1,0.1};
	//tabSources[1].Ip = {0.7,0.7,0.7};
	//tabSources[1].pos = {-0,10,2};
	//nbSources++;
}


bool findNextIntersection(Point start, Vector dir, const Sphere* current_sphere, Point& pt, Face& face)
{
	bool intersection = false;
	double tmin = HUGE_VALF;
	for(int sphereIter = 0; sphereIter<nbSpheres; ++sphereIter)
	{
		Sphere *sphere = &tabSpheres[sphereIter];
		if(sphere == current_sphere)
			continue;
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
				face.sphere = sphere;
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
					face.sphere = sphere;
					face.normal = normal;
				}
			}
		}
	}
	return intersection;
}


bool ombre(const Point& start,const Vector& dir,const Face& face, const Source& src)
{
	double d_src = dist(start,src.pos);
	for(int sphereIter = 0; sphereIter<nbSpheres; ++sphereIter)
	{
		Sphere *sphere = &tabSpheres[sphereIter];
		if(face.sphere == sphere)
			continue;

		double A = dir.x*dir.x + dir.y*dir.y + dir.z*dir.z;
		double B = 2*(dir.x*(start.x-sphere->xc)+dir.y*(start.y-sphere->yc)+dir.z*(start.z-sphere->zc));
		double C = (start.x-sphere->xc)*(start.x-sphere->xc)+
		           (start.y-sphere->yc)*(start.y-sphere->yc)+
		           (start.z-sphere->zc)*(start.z-sphere->zc)-
		           sphere->r*sphere->r;
		double delta = B*B-4*A*C;
		if(delta<0)
			continue;
		else if(delta == 0)
		{
			double t = -B / (2*A);
			if(t>0 && t<d_src)
				return true;
		}else
		{
			double t1 = (-B-sqrt(delta))/(2*A);
			double t2 = (-B+sqrt(delta))/(2*A);
			if((t1>0 && t1<d_src) || (t2>0 && t2<d_src))
				return true;
		}
	}
	return false;
}



void raytrace2(Point start, Vector dir, const Sphere* current_sphere, double current_ro, int nb, Intensity& I)
{
	if(nb == 0)
	{
		I = {0,0,0};
	}else
	{
		Point pointI;
		Face face;
		if(!findNextIntersection(start,dir,current_sphere,pointI,face))
		{
			I = {0,0,0};
		}else
		{

			Intensity Is = {0,0,0}, It = {0,0,0};
			Vector Lj, vision, reflected, transmitted;
			vision.fromPoints(pointI,obs).normalize();

			if(face.sphere->mat.ks.r == 0 && face.sphere->mat.ks.g == 0 && face.sphere->mat.ks.b == 0)
			{
				Is = {0, 0, 0};
			}else
			{
				if(face.normal*Lj>0)
				{
					// compute reflected ray
					reflected = (face.normal * 2 *(face.normal * Lj)) - Lj;
					//reflected.normalize();
					raytrace2(pointI,reflected,face.sphere,current_ro,nb-1,Is);
				}
			}


			//if(face.sphere->mat.kt.r != 0 || face.sphere->mat.kt.g != 0 || face.sphere->mat.kt.b != 0 )
			//{
			//	Vector L = dir*-1;
			//	if(face.normal * L > 0)
			//	{
			//		double ro = current_ro / face.sphere->mat.ro;
			//
			//		double tmp = ro*(face.normal*L) - sqrt(1-ro*ro*(1-pow(face.normal*L,2)));
			//		transmitted = (face.normal*tmp)-(L*ro);
			//		//transmitted.normalize();
			//		//raytrace2(pointI,transmitted,NULL,face.sphere->mat.ro,nb-1,It);
			//		std::cout << "############\n################\n##############\n#############\n";
			//	}
			//
			//}else
			//{
			//	It = {0,0,0};
			//}
	int xxx;double aaa=0,bbb=0,ccc=0,rt=0,are=0,fr,s,bdf,o;
	//Vector aaaare = dir;
			//Vector Lll = {-dir.x,-dir.y,-dir.z};
			//double ro_it = current_ro/face.sphere->mat.ro;
			//double tmp1 = (face.normal * L) * ro_it;
			//double tmp2 = sqrt(1-(ro_it*ro_it*(1-pow(face.normal * L,2))));


			//ambient only with the first source
			I = source1->Ia * face.sphere->mat.ka;
			// diffuse with all the sources
			for(int srcIter = 0; srcIter<nbSources; ++srcIter)
			{
				Lj.fromPoints(pointI,tabSources[srcIter].pos).normalize();
				if(face.normal * Lj >0)	//eclaire les faces bien orientees par rapport a la source
				{
					if(ombre(pointI,Lj,face,tabSources[srcIter]))
						continue;

					double d = dist(tabSources[srcIter].pos,pointI);
					double fAtt = 1/(0.01*d*d+0.04*d+0.03);
					if(fAtt>1) fAtt = 1;

					// diffuse
					I +=  (tabSources[srcIter].Ip * fAtt * face.sphere->mat.kd) * vmax((face.normal * Lj),0);

					// specular
					if(face.sphere->mat.ks.r != 0 || face.sphere->mat.ks.g != 0 || face.sphere->mat.ks.b != 0)
						I += (tabSources[srcIter].Ip * fAtt * face.sphere->mat.ks) * pow(vmax(reflected*vision,0), (double)face.sphere->mat.shininess);
				}
			}

			I += Is * face.sphere->mat.ks;
			//I += It * face.sphere->mat.kt;

			if(start != obs)
			{
				double d = dist(start,pointI);
				I *= (1/(1+d));
			}
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
			if(!antialiasing)
			{
				direction = {left + j * stepH, top - i * stepV, -znear};
				direction.normalize();
				raytrace2(obs,direction,NULL,1.0,2,it);

				setPixel(screen,i,j,floor(it.r*255), floor(it.g*255),floor(it.b*255));
			}else
			{	//anti-aliasing
				double rm=0,gm=0,bm=0;
				double offset = 0.002;
				direction = {left + j * stepH, top - i * stepV, -znear};
				raytrace2(obs,direction,NULL,1.0,2,it);
				rm += it.r; gm +=it.g; bm += it.b;

				direction = {left + j * stepH -offset, top - i * stepV, -znear};
				raytrace2(obs,direction,NULL,1.0,2,it);
				rm += it.r; gm +=it.g; bm += it.b;

				direction = {left + j * stepH + offset, top - i * stepV, -znear};
				raytrace2(obs,direction,NULL,1.0,2,it);
				rm += it.r; gm +=it.g; bm += it.b;

				direction = {left + j * stepH, top - i * stepV - offset, -znear};
				raytrace2(obs,direction,NULL,1.0,2,it);
				rm += it.r; gm +=it.g; bm += it.b;

				direction = {left + j * stepH, top - i * stepV + offset, -znear};
				raytrace2(obs,direction,NULL,1.0,2,it);
				rm += it.r; gm +=it.g; bm += it.b;

				rm /= 5; gm /= 5; bm /= 5;
				setPixel(screen,i,j,floor(rm*255), floor(gm*255),floor(bm*255));

			}
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

