#include "euler.h"
void rr_euler(double* y, void (*deriv_func)(double, double*, int, double*), double x, const int n, const double h) {
    int i;
    double dydx[n];
    deriv_func(x,y,n,dydx);
    for(i=0;i<n;i++) {
	y[i]+=dydx[i]*h;
    }
}
