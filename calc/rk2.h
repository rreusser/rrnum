#ifndef __RK2_H__
#define __RK2_H__

void rr_rk2(double* y, void (*deriv_func)(double, double*, int, double*), double x, const int n, const double h);

#endif
