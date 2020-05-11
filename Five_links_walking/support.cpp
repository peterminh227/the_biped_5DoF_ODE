#include <ode/ode.h>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "libxl.h"
#include <iomanip>
#include <locale>
#include <armadillo> // algebra library
#include <cmath>
#include <drawstuff/drawstuff.h>
#include "support.h"
int sign(dReal value)
{
	int sign;
	if (value >= 0) sign = 1;
	if (value < 0) sign = -1;
	return sign;
}
double min(dReal x, dReal y)
{
	if (x <= y)
	{
		return x;
	}
	return y;
}
double max(dReal x, dReal y)
{
	if (x >= y)
	{
		return x;
	}
	return y;
}
double getRotateAngle(dVector3 vec1, dVector3 vec2)
{
	const double epsilon = 1.0e-6;
	double angle = 0;
	// normalize
	dVector3 norVec1, norVec2;
	norVec1[0] = vec1[0] / sqrt(pow(vec1[0],2)+ pow(vec1[1],2) + pow(vec1[2],2));
	norVec1[1] = vec1[1] / sqrt(pow(vec1[0],2)+ pow(vec1[1],2) + pow(vec1[2],2));
	norVec1[2] = vec1[2] / sqrt(pow(vec1[0],2)+ pow(vec1[1],2) + pow(vec1[2],2));
	norVec2[0] = vec2[0] / sqrt(pow(vec2[0],2)+ pow(vec2[1],2) + pow(vec2[2],2));
	norVec2[1] = vec2[1] / sqrt(pow(vec2[0],2)+ pow(vec2[1],2) + pow(vec2[2],2));
	norVec2[2] = vec2[2] / sqrt(pow(vec2[0],2)+ pow(vec2[1],2) + pow(vec2[2],2));
	// dot product
	double dotProd = (norVec1[0]*norVec2[0] + norVec1[1]*norVec2[1] + norVec1[2]*norVec2[2]);
	if ( abs(dotProd - 1.0) <= epsilon )
	angle = 0;
	else if ( abs(dotProd + 1.0) <= epsilon )
	angle = PI;
	else {
		double cross_x = 0;
		angle = acos(dotProd);
		cross_x = (norVec1[0]*norVec2[2]-norVec1[2]*norVec2[0]);
		if (cross_x < 0) // vec1 rotate clockwise to vec2
		angle = 2 * PI - angle; 

	}
	return angle;
}
double deg2rad(dReal angle)
{
	dReal radian;
	radian = angle/180*PI;
	return radian;
}

//// 
