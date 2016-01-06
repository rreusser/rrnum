#include "ab2.h"

void rr_ab2(double* y, void (*deriv_func)(double, double*, int, double*), double x, const int n, const double h, double* ypold) {
    int i;
    double dydx[n];
    deriv_func(x,y,n,dydx);
    for(i=0;i<n;i++) {
	y[i]+=(1.5*dydx[i]-0.5*ypold[i])*h;
	ypold[i]=dydx[i];
    }
}
