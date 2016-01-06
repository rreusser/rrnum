#ifndef __RK4_H__
#define __RK4_H__
#define RK4 4

void rr_rk4(double* y, void (*deriv_func)(double, double*, int, double*), double x, const int n, const double h);

#endif
