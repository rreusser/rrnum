#ifndef __NDER_H__
#define __NDER_H__

double rr_nder(double (*funcPtr)(double), const double x0, const double tol);

double rr_nder_param(double (*funcPtr)(double, double*), const double x0, const double tol, double* param);

double dummy(double a);

#endif
