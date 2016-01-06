#ifndef __FUNC_H__
#define __FUNC_H__

#include <math.h>

#define SQR(A) ((A)*(A))
#define MAX(A,B) ((A)>(B)?(A):(B))
#define MIN(A,B) ((A)<(B)?(A):(B))
#define ABS(A) ((A)<0?(-(A)):(A))
#define SGN(A) ((A)==0?0:((A)>0?1:-1))


double arctan(double x, double y);
/*float fmod(float A, float B);*/
int imax(int A, int B);
int imin(int A, int B);
double dmin(double A, double B);
double dmax(double A, double B);
float min(float A, float B);
float max(float A, float B);
/*double abs(double A);*/
void swap(int *A, int *B);

double rrDArrayMin(double* a, int n);
double rrDArrayMax(double* a, int n);
int rrDArrayMinIndex(double* a, int n);
int rrDArrayMaxIndex(double* a, int n);
double dmax3(double A, double B, double C);
double dmin3(double A, double B, double C);
double dabs(double A);
int dfloor(double A);
int toggle(int* A);
#endif
