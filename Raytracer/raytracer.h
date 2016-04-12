#ifndef RAYTRACER_H
#define RAYTRACER_H

#include <cmath>

struct Coef{
	double r;
	double g;
	double b;
};

struct Intensity{
	double r;
	double g;
	double b;

	Intensity operator*(double k);
	Intensity operator*(const Coef& k);
	Intensity operator+(const Intensity& it);
	void operator+=(const Intensity& it);
};

struct Point{
	double x;
	double y;
	double z;
};

struct Vector{
	double x;
	double y;
	double z;

	double operator*(const Vector& v2);
	void normalize();
	Vector& fromPoints(const Point& a, const Point& b);
};

struct Sphere{
	double xc;
	double yc;
	double zc;
	double r;
	double g;
	double b;
	Coef ka;
	Coef kd;
	Coef ks;
};

struct Face{	//represents an intersection between an object and a ray
	Vector normal;
	Coef ka;
	Coef kd;
	Coef ks;
	Coef kro;
};

struct Source{
	Point pos;
	Intensity Ia;
	Intensity Ip;
};



inline double dist(const Point &a, const Point &b)
{
	return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)+(a.z-b.z)*(a.z-b.z));
}



#endif // RAYTRACER_H

