#ifndef __GNUPLOTIO_H__
#define __GNUPLOTIO_H__
#include <stdio.h>
#include <rrgrid.h>
#include <rrdata.h>

void rr_gnuplotODE(char *fname, double **y, char *mname, const rr_ibase* ibase, int n);

void rr_gnuplotXYArray(char *fname, double *x, double *y, int nsteps);

void rr_gnuplotArray(char *fname, double *y, char *mname, const rr_ibase* ibase);

void rr_gnuplotXYFunc(char *fname, double (*xfuncPtr)(double), double(*yfuncPtr)(double), const rr_ibase* ibase);

void rr_gnuplotFunc(char *fname, double (*funcPtr)(double), const rr_ibase* ibase);


#endif
