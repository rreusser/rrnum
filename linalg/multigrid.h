#ifndef __MGRID_H__
#define __MGRID_H__
#include <rrgeneral.h>
#include <stdlib.h>

typedef struct {
    int lambda; //G-S iterations per level
    int maxlevels;
    int* isize;
    int* jsize;
    int* inc;
    double** du;
    double** pdu;
    double** r;
    double** Ldu;
    double** A;
    double** M;
    rr_structured* grid;
} rrmgSolver;


void rrmgAlloc(int isize, int jsize, int maxlevels);

void mgFree(double** mgdu, double** mgpdu, double** mgr, double** Ldu, int n);

double mgGaussSeidelZero(double* x, double h, int n, int iters, double omega);
double mgGaussSeidelNeg(double* x, double* b, double h, int n, int iters, double omega);
double mgGaussSeidel(double* x, double* b, double h, int n, int iters, double omega);
void mgResid(double* x, double* b, double h, int n);


void mgRestrict(double* rr, double* r, int n);
void mgProlongate(double* pdu, double* du, int n);

double mgLevel(double** du, double** pdu, double** r, double** Ldu, double h0, int n, double omega, int minlevel, int lambda);

double mgIter(double* Field, double** pdu, double** du, double** r, double** Ldu, double h0, int n, double tol, int itermax, double omega, int minlevel, int lambda);

#endif
