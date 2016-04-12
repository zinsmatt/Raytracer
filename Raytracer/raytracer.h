#ifndef RAYTRACER_H
#define RAYTRACER_H

#include <cmath>

struct Coef{
	double r;
	double g;
	double b;
};

struct Material{
	Coef ka;
	Coef kd;
	Coef ks;
	unsigned int shininess;
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
	Vector operator*(double k);
	Vector operator-(const Vector& v2);
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
	Material mat;

};

struct Face{	//represents an intersection between an object and a ray
	Vector normal;
	Sphere* sphere;
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

inline double power(double x, int n)
{
	double res = 1.0;
	for(int i=0;i<n;++i)
		res *= x;
	return res;
}


#endif // RAYTRACER_H

