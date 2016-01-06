#include "rk2.h"

void rr_rk2(double* y, void (*deriv_func)(double, double*, int, double*), double x, const int n, const double h) {
    int i;
    double k1[n], k2[n], ym[n];
    double xm;
    deriv_func(x,y,n,k1);
    for(i=0;i<n;i++) {
	ym[i]=y[i]+k1[i]*0.5*h;
    }
    xm=x+0.5*h;
    deriv_func(xm,ym,n,k2);
    for(i=0;i<n;i++) {
	y[i]=y[i]+k2[i]*h;
    }
}
