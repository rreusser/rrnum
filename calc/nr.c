#include "nr.h"
#include "nder.h"

double rr_newton_raphson(double (*func)(double), double guess, double tol, double dx, int nmax, int verbose) {
    int n=0;
    double err,fp,pguess=0.0,f,xt;

    if(verbose) {
	printf("Newton-Raphson iterations:\n");
    }
    err=tol+1.0;
    while(err > tol && n<nmax) {
	fp = rr_nder(func,guess,dx);
	f=func(guess);
	xt=guess;
	guess -= f/fp;
	if(n>nmax) break;
	if(n>0) {
	    err=guess-pguess;
	    err=err<0?-err:err;
	    if(verbose) {
		printf("iter=%i,  x=%f,  f(x)=%f,  err=%f\n",n+1,xt,f,err);
	    }
	} else {
	    if(verbose) {
		printf("iter=%i,  x=%f,  f(x)=%f\n",n+1,xt,f);
	    }
	}
	pguess=guess;
	n++;
    }

    return guess;
}

double rr_newton_raphson_param(double (*func)(double,double*), double guess, double tol, double dx, int nmax, int verbose, double* param) {
    int n=0;
    double err,fp,pguess=0.0,f,xt;

    if(verbose) {
	printf("Newton-Raphson iterations:\n");
    }
    err=tol+1.0;
    while(err > tol && n<nmax) {
	fp = rr_nder_param(func,guess,dx,param);
	f=func(guess,param);
	xt=guess;
	guess -= f/fp;
	if(n>nmax) break;
	if(n>0) {
	    err=guess-pguess;
	    err=err<0?-err:err;
	    if(verbose) {
		printf("iter=%i,  x=%e,  f(x)=%e,  err=%e\n",n+1,xt,f,err);
	    }
	} else {
	    if(verbose) {
		printf("iter=%i,  x=%e,  f(x)=%e\n",n+1,xt,f);
	    }
	}
	pguess=guess;
	n++;
    }

    return guess;
}

double rr_newton_raphson_param_bounded(double (*func)(double,double*), double guess, double lower, double upper, double tol, double dx, int nmax, int verbose,double* param) {
    double bound_tol = 1.0e-5;
    int n=0;
    double err,fp,pguess=0.0,f,xt;
    if(verbose) {
	printf("Newton-Raphson iterations:\n");
    }
    err=tol+1.0;
    while(err > tol && n<nmax) {
	fp = rr_nder_param(func,guess,dx,param);
	f=func(guess,param);
	xt=guess;
	guess -= f/fp;
	if(n>nmax) break;
	if(n>0) {
	    err=guess-pguess;
	    err=err<0?-err:err;
	    if(verbose) {
		printf("iter=%i,  x=%e,  f(x)=%e,  err=%e\n",n+1,xt,f,err);
	    }
	} else {
	    if(verbose) {
		printf("iter=%i,  x=%e,  f(x)=%e\n",n+1,xt,f);
	    }
	}
	pguess=guess;
	n++;
    }
    if(guess>upper+bound_tol) {
	fprintf(stderr,"Newton-Raphson: %e is out of bounds. Defaulting to upper bound, %e\n", guess, upper);
	guess = upper;
    }
    if(guess<lower-bound_tol) {
	fprintf(stderr,"Newton-Raphson: %e is out of bounds. Defaulting to lower bound, %e.\n", guess, lower);
	guess = lower;
    }

    return guess;
}
