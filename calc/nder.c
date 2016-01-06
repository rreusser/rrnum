#include "nder.h"

double rr_nder(double (*funcPtr)(double), const double x0, const double tol) {
    double ans = (funcPtr(x0+tol)-funcPtr(x0-tol))/(2.0*tol);
    return ans;
}

double rr_nder_param(double (*funcPtr)(double,double*), const double x0, const double tol, double* param) {
    double ans = (funcPtr(x0+tol,param)-funcPtr(x0-tol,param))/(2.0*tol);
    return ans;
}

double dummy(double a) {
    return -a;
}
