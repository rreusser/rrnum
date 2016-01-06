
#include "rk4.h"
void rr_rk4(double* y, void (*deriv_func)(double, double*, int, double*), double x, const int n, const double h) {
    int i;
    double k1[n], k2[n], k3[n], k4[n], ym[n];
    double xm, slope;
    deriv_func(x,y,n,k1);
    for(i=0;i<n;i++) {
	ym[i]=y[i]+k1[i]*0.5*h;
    }
    xm=x+0.5*h;
    deriv_func(xm,ym,n,k2);
    for(i=0;i<n;i++) {
	ym[i]=y[i]+k2[i]*0.5*h;
    }
    deriv_func(xm,ym,n,k3);
    for(i=0;i<n;i++) {
	ym[i]=y[i]+k3[i]*h;
    }
    xm=x+h;
    deriv_func(xm,ym,n,k4);
    for(i=0;i<n;i++) {
	slope=(k1[i]+2.0*(k2[i]+k3[i])+k4[i])/6.0;
	y[i]=y[i]+slope*h;
    }
}
