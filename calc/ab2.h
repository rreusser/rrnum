#ifndef __AB2_H__
#define __AB2_H__

void rr_ab2(double* y, void (*deriv_func)(double, double*, int, double*), double x, const int n, const double h, double* ypold);

#endif
