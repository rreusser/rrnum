#ifndef __INTEGRATE_H__
#define __INTEGRATE_H__
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <rrdata.h>

#ifndef __EULER_H__
#include "euler.h"
#include "rk2.h"
#include "rk4.h"
#include "ab2.h"
#endif

#define EULER 1
#define RK2   2
#define RK4   4
#define AB2   12

int rr_integrator_type(const char* type);

double** rr_integrate_malloc(const rr_ibase* tbase, const int n);
void rr_integrate_free(double **y, const int n);

void rr_integrate(const char* cintegrator, double* y, void (*deriv_func)(double,double*,int,double*),const rr_ibase* tbase, const int n);

void rr_integrate_store(const char* cintegrator, double** y, void (*deriv_func)(double,double*,int,double*),const rr_ibase* tbase, const int n);

#endif
