#ifndef __TECIO_H__
#define __TECIO_H__
#include <stdio.h>
#include <rrgrid.h>
#include <rrdata.h>

void rr_tecplotConformalMap(char *fname, double (*xfuncPtr)(double,double), double (*yfuncPtr)(double,double), const rr_ibase* xbase, const rr_ibase* ybase, char *mname);

void rr_tecplotODE(char *fname, double **y, char *mname, const rr_ibase* ibase, int n);

void rr_tecplotXYArray(char *fname, double *x, double *y, int nsteps);

void rr_tecplotArray(char *fname, double *y, char *mname, const rr_ibase* ibase);

void rr_tecplotXYFunc(char *fname, double (*xfuncPtr)(double), double(*yfuncPtr)(double), char *mname, const rr_ibase* ibase);

void rr_tecplotFunc(char *fname, double (*funcPtr)(double), char *mname, const rr_ibase* ibase);

void rr_tecplotStructuredGrid(char *fname, rr_structured_grid* grid, double* T);

void rr_tecplotPolar(char *fname, double (*rFunc)(double), double (*thetaFunc)(double), const rr_ibase* tbase);

#endif
