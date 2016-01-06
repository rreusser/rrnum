#include "stencil.h"

//Use this routine to automatically work around boundary conditions,
//at the expense of precise stencil patterns.  It just chooses the most
//symmetric pattern available, or asymmetric if the wall interferes.
void rr_choose_points(int npoints,rr_ibase* ibase, int i, int* ip, double* dp) {
    int j;
    double center;
    center=ibase->t[i];
    switch(npoints) {
    default:
	printf("choose points: error: unexpected number of points\n");
    break;
    case (2):
	if(i<0) {
	    printf("error: i<0\n");
	    ip[0]=ip[1]=0;
	} else if(i<1) {
	    ip[0]=0;
	    ip[1]=1;
	} else if(i<ibase->n-1) {
	    ip[0]=i-1;
	    ip[1]=i+1;
	} else if(i<ibase->n) {
	    ip[0]=ibase->n-2;
	    ip[1]=ibase->n-1;
	} else {
	    printf("Error, i out of bounds.\n");
	    ip[0]=ip[1]=0;
	}
	for(j=0;j<2;j++) {
	    dp[j]=ibase->t[ip[j]]-center;
	}
    break;
    case(3):
	if(i<0) {
	    printf("error: i<0\n");
	    ip[0]=ip[1]=ip[2]=0;
	} else if(i<1) {
	    ip[0]=0;
	    ip[1]=1;
	    ip[2]=2;
	} else if(i<ibase->n-1) {
	    ip[0]=i-1;
	    ip[1]=i;
	    ip[2]=i+1;
	} else if(i<ibase->n) {
	    ip[0]=ibase->n-3;
	    ip[1]=ibase->n-2;
	    ip[2]=ibase->n-1;
	} else {
	    printf("Error, i out of bounds.\n");
	    ip[0]=ip[1]=ip[2]=0;
	}
	for(j=0;j<3;j++) {
	    dp[j]=ibase->t[ip[j]]-center;
	}
    break;
    case (4):
	if(i<0) {
	    printf("error: i<0\n");
	    ip[0]=ip[1]=ip[2]=ip[3]=0;
	} else if(i<2) {
	    ip[0]=0;
	    ip[1]=1;
	    ip[2]=2;
	    ip[3]=3;
	} else if(i<ibase->n-2) {
	    ip[0]=i-2;
	    ip[1]=i-1;
	    ip[2]=i+1;
	    ip[3]=i+2;
	} else if(i<ibase->n) {
	    ip[0]=ibase->n-4;
	    ip[1]=ibase->n-3;
	    ip[2]=ibase->n-2;
	    ip[3]=ibase->n-1;
	} else {
	    printf("Error, i out of bounds.\n");
	    ip[0]=ip[1]=ip[2]=ip[3]=0;
	}
	for(j=0;j<4;j++) {
	    dp[j]=ibase->t[ip[j]]-center;
	}
    break;
    }
}

//General solutions for dy/dx so they can be used for non-uniform
//grids - obtained by differentiating polynomial fits.
void rr_dy_dx(int npoints, double *xp, double *an) {
    double x0,x1,x2;
    switch(npoints) {
    default:
	printf("stencil: error: unexpected number of points\n");
    break;
    case(2):
	an[0]=-1.0/(xp[1]-xp[0]);
	an[1]= 1.0/(xp[1]-xp[0]);
    break;
    case(3):
	x0=xp[0];
	x1=xp[1];
	x2=xp[2];
	an[0]=-(x2+x1)/(x1-x0)/(x2-x0);
	an[1]= (x2+x0)/(x1-x0)/(x2-x1);
	an[2]=-(x1+x0)/(x2-x0)/(x2-x1);
    break;
    }
}

//General solutions for d2y/dx2 so they can be used for non-uniform
//grids - obtained by differentiating polynomial fits.
void rr_d2y_dx2(int npoints, double* xp, double *an) {
    double x0,x1,x2,x3;
    switch(npoints) {
    default:
	printf("stencil: error: unexpected number of points\n");
    break;
    case(3):
	x0=xp[0];
	x1=xp[1];
	x2=xp[2];
	an[0]= 2.0/(x1-x0)/(x2-x0);
	an[1]=-2.0/(x1-x0)/(x2-x1);
	an[2]= 2.0/(x2-x0)/(x2-x1);
    break;
    case(4):
	x0=xp[0];
	x1=xp[1];
	x2=xp[2];
	x3=xp[3]; //Thank you Maxima...
	an[0]= 2.0*(x3+x2+x1)/(x1-x0)/(x2-x0)/(x3-x0);
	an[1]=-2.0*(x3+x2+x0)/(x1-x0)/(x2-x1)/(x3-x1);
	an[2]= 2.0*(x3+x1+x0)/(x2-x0)/(x2-x1)/(x3-x2);
	an[3]=-2.0*(x2+x1+x0)/(x3-x0)/(x3-x1)/(x3-x2);
    break;
    }
}

