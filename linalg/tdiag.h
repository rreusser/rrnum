#ifndef __TDIAG_H__
#define __TDIAG_H__

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void rr_tridiag(double** A, double* b, const int n);
void rr_pentadiag(double** A, double* b, const int n);
void rr_tridiag_periodic(double** A, double* b, const int n, double* work1, double* work2);

#endif
