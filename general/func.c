#include "func.h"

#define PI 3.14159265358979323846264

double arctan(double x, double y) { /*use x,y coords to find angle from 0 to 2*pi instead of the usual -pi/2 to pi/2*/
    double answer=0.0;
    if (x>0 && y>0) answer=atan(y/x);
    else if (x==0 && y>0) answer=PI/2.0;
    else if (x==0 && y<0) answer=PI*3.0/2.0;
    else if (x==0 && y==0) answer=0.0;
    else if (x<0) answer=PI+atan(y/x);
    else if (y<0) answer=atan(y/x)+2.0*PI;
    return answer;
}

int toggle(int* A) {
    *A=1-*A;
    return *A;
}

/*float fmod(float A, float B) {
	return A*float(int(A/B)%1);
}*/

int rrDArrayMinIndex(double* a, int n) {
    double themin = a[0];
    int minind=0;
    int i;
    for(i=1; i<n; i++) {
	if(themin>a[i]) {
	    themin = a[i];
	    minind = i;
	}
    }
    return minind;
}
int rrDArrayMaxIndex(double* a, int n) {
    double themax = a[0];
    int maxind=0;
    int i;
    for(i=1; i<n; i++) {
	if(themax<a[i]) {
	    themax = a[i];
	    maxind = i;
	}
    }
    return maxind;
}
double rrDArrayMax(double* a, int n) {
    double themax = a[0];
    int i;
    for(i=1; i<n; i++) {
	themax = dmax(themax,a[i]);
    }
    return themax;
}
double rrDArrayMin(double* a, int n) {
    double themin = a[0];
    int i;
    for(i=1; i<n; i++) {
	themin = dmin(themin,a[i]);
    }
    return themin;
}
int imax(int A, int B) {
	return (A>B)?A:B;
}
int imin(int A, int B) {
	return (A<B)?A:B;
}
double dmin(double A, double B) {
	return (A<B)?A:B;
}
double dmax(double A, double B) {
	return (A>B)?A:B;
}
float min(float A, float B) {
	return (A<B)?A:B;
}
float max(float A, float B) {
	return (A>B)?A:B;
}

double dmax3(double A, double B, double C) {
    double maxBC = (B>C)?(B):(C);
    return (A>maxBC)?(A):(maxBC);
}
double dmin3(double A, double B, double C) {
    double minBC = (B<C)?(B):(C);
    return (A<minBC)?(A):(minBC);
}
double dabs(double A){ return (A<0.0)?(-A):(A); }

int dfloor(double A) {
    return (A<0)?(-(int)(-A+1.0)):((int)A);
}
