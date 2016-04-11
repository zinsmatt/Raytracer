#ifndef RAYTRACER_H
#define RAYTRACER_H

#include <cmath>

struct Sphere{
	double xc;
	double yc;
	double zc;
	double r;
	double g;
	double b;
	double ka;
	double kd;
	double ks;
};

struct Intensity{
	double r;
	double g;
	double b;
//	double a;
};

struct Vector{
	double x;
	double y;
	double z;
};

struct Point{
	double x;
	double y;
	double z;
};

struct Face{	//represents an intersection between an object and a ray
	Vector normal;
	double ka;
	double kd;
	double kro;
};

struct Intersection{
	Point pt;
	Face face;
};

struct Source{
	Point pos;
	Intensity Ia;
	Intensity Ip;
};


Intensity mult(const Intensity &it, double k)
{
	Intensity temp = it;
	temp.r *= k;
	temp.g *= k;
	temp.b *= k;
	return temp;
}

void add(Intensity &it, const Intensity &it2)
{
	it.r += it2.r;
	it.g += it2.g;
	it.b += it2.b;
}

double dist(const Point &a, const Point &b)
{
	return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)+(a.z-b.z)*(a.z-b.z));
}

double prodScalPos(const Vector& a, const Vector& b)
{
	double res = a.x*b.x + a.y*b.y + a.z*b.z;
	return (res>0)?res:0;
}

void normalize(Vector &v)
{
	double n = sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
	v.x /= n;
	v.y /= n;
	v.z /= n;
}

Vector vectNormFromPoints(const Point &a, const Point &b)
{
// return the normalized vector from a to b
	Vector v = {b.x-a.x,b.y-a.y,b.z-a.z};
	normalize(v);
	return v;
}


#endif // RAYTRACER_H

