#include "raytracer.h"




Intensity Intensity::operator*(double k)
{
	Intensity temp = *this;
	temp.r *= k;
	if(temp.r>1) temp.r = 1.0;
	temp.g *= k;
	if(temp.g>1) temp.g = 1.0;
	temp.b *= k;
	if(temp.b>1) temp.b = 1.0;
	return temp;
}

void Intensity::operator*=(double k)
{
	r *= k;
	if(r>1) r=1;
	g *= k;
	if(g>1) g=1;
	b *= k;
	if(b>1) b=1;
}

Intensity Intensity::operator*(const Coef& k)
{
	Intensity temp = *this;
	temp.r *= k.r;
	if(temp.r>1) temp.r = 1.0;
	temp.g *= k.g;
	if(temp.g>1) temp.g = 1.0;
	temp.b *= k.b;
	if(temp.b>1) temp.b = 1.0;
	return temp;
}

Intensity Intensity::operator+(const Intensity& it)
{
	Intensity temp;
	temp.r = this->r + it.r;
	if(temp.r>1) temp.r = 1.0;
	temp.g = this->g + it.g;
	if(temp.g>1) temp.g = 1.0;
	temp.b = this->b + it.b;
	if(temp.b>1) temp.b = 1.0;
	return temp;
}

void Intensity::operator+=(const Intensity& it)
{
	this->r = this->r + it.r;
	if(this->r>1) this->r = 1.0;
	this->g = this->g + it.g;
	if(this->g>1) this->g = 1.0;
	this->b = this->b + it.b;
	if(this->b>1) this->b = 1.0;
}

double Vector::operator*(const Vector& v2)
{
	double res = this->x*v2.x + this->y*v2.y + this->z*v2.z;
	return res;
}

Vector Vector::operator*(double k)
{
	Vector temp = *this;
	temp.x *= k;
	temp.y *= k;
	temp.z *= k;
	return temp;
}

Vector Vector::operator-(const Vector& v2)
{
	Vector temp;
	temp.x = x - v2.x;
	temp.y = y - v2.y;
	temp.z = z - v2.z;
	return temp;
}

void Vector::normalize()
{
	double n = sqrt(x*x+y*y+z*z);
	x /= n;
	y /= n;
	z /= n;
}

Vector& Vector::fromPoints(const Point &a, const Point &b)
{
	x = b.x-a.x;
	y = b.y-a.y;
	z = b.z - a.z;
	return *this;
}
