#ifndef __EULER_H__
#define __EULER_H__
#include <stdlib.h>

void rr_euler(double* y, void (*deriv_func)(double, double*, int, double*), double x, const int n, const double h);

#endif
