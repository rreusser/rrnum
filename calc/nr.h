#ifndef __NR_H__
#define __NR_H__
#include <stdio.h>

double rr_newton_raphson(double (*func)(double), double guess, double tol, double dx, int nmax, int verbose);

double rr_newton_raphson_param(double (*func)(double,double*), double guess, double tol, double dx, int nmax, int verbose,double* param);

double rr_newton_raphson_param_bounded(double (*func)(double,double*), double guess, double lower, double upper, double tol, double dx, int nmax, int verbose,double* param);

#endif
