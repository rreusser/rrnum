#include "integrate.h"

int rr_integrator_type(const char* type) {
    int itype=1;
    if(strcasecmp(type,"EULER")==0) {
	itype=1;
    } else if(strcasecmp(type,"RK2")==0) {
	itype=2;
    } else if(strcasecmp(type,"RK4")==0) {
	itype=4;
    } else if(strcasecmp(type,"AB2")==0) {
	itype=12;
    } else {
	printf("rr_integrator_type: Unknown type '%s': defaulting to EULER.\n",type);
    }
    return itype;
}

    
//Automatically malloc the right sized array for use with rr_integrate_store
double** rr_integrate_malloc(const rr_ibase* tbase, const int n) {
    const int nsteps = tbase->n;
    double** temp;
    int i;
    temp = malloc(n*sizeof(double*));
    for(i=0; i<n; i++) {
	temp[i]=malloc(nsteps*sizeof(double));
    }
    return temp;
}

//Automatically free the array
void rr_integrate_free(double **y, const int n) {
    int i;
    for(i=0; i<n; i++) {
	free(y[i]);
    }
    free(y);
}

//March an ODE forward in time from t0 to t1 in 'nsteps' steps.
//Pass the initial conditions of y_i in y[i][0] and store the nth
//integration step in y[i][n].  Dimensions of y must be y[i][nsteps+1]
void rr_integrate_store(const char* cintegrator, double** y, void (*deriv_func)(double,double*,int,double*),const rr_ibase* tbase, const int n) {

    const double dt = tbase->dt;
    const double t0 = tbase->t0;
    const int nsteps = tbase->n-1;
    double t;
    double ym[n];
    double *ypold;
    int i,j;
    t = tbase->t0;
    
    //transfer the initial y vector into ym
    //in order to pass to the actual integrator
    for(i=0; i<n; i++) {
	ym[i] = y[i][0];
    }

    switch(rr_integrator_type(cintegrator)) {
    default:
	printf("rr_integrate: Defaulting to Euler integration\n");

    case EULER:
	for(i=0; i<nsteps; i++) {
	    rr_euler( ym, deriv_func, t, n, dt );
	    for(j=0; j<n; j++) {
		y[j][i+1]=ym[j];
	    }
	    t=t0+dt*(i+1);
	    if(tbase->verbosity) printf("%i \t%f \t%f\n",i,t,y[0][i+1]);
	}
    break;
    
    case RK2:
	for(i=0; i<nsteps; i++) {
	    rr_rk2( ym, deriv_func, t, n, dt );
	    for(j=0; j<n; j++) {
		y[j][i+1]=ym[j];
	    }
	    t=t0+dt*(i+1);
	    if(tbase->verbosity) printf("%i \t%f \t%f\n",i,t,y[0][i+1]);
	}
    break;

    case RK4:
	for(i=0; i<nsteps; i++) {
	    rr_rk4( ym, deriv_func, t, n, dt );
	    for(j=0; j<n; j++) {
		y[j][i+1]=ym[j];
	    }
	    t=t0+dt*(i+1);
	    if(tbase->verbosity) printf("%i \t%f \t%f\n",i,t,y[0][i+1]);
	}
    break;

    case AB2:

	ypold=malloc(sizeof(double)*n);

	//initialize with an euler step
	deriv_func(t,ym,n,ypold);
	rr_euler( ym, deriv_func, t, n, dt );
	for(j=0; j<n; j++) {
	    y[j][1]=ym[j];
	}
	t=t0+dt;
	//now jump into adams-bashforth integration
	for(i=1; i<nsteps; i++) {
	    rr_ab2( ym, deriv_func, t, n, dt, ypold );
	    for(j=0; j<n; j++) {
		y[j][i+1]=ym[j];
	    }
	    t=t0+dt*(i+1);
	}
	free(ypold);
    break;
    }
}

//March an ODE forward in time from t0 to t1 in 'nsteps' steps.
//Pass the initial conditions of y_i in y[i] and only store the
//final integration step.
void rr_integrate(const char* cintegrator, double* y, void (*deriv_func)(double,double*,int,double*),const rr_ibase* tbase, const int n) {
    const double dt = tbase->dt;
    const double t0 = tbase->t0;
    const int nsteps = tbase->n-1;
    double t, *ypold;
    double ym[n];
    int i,j;
    t=t0;
    //printf("Integrator: \n\tdt = %f\n\tt0 = %f\n\tt1 = %f\n",dt,t0,t1);
    if(tbase->verbosity) printf("%f \t%f\n",t,y[0]);

    switch(rr_integrator_type(cintegrator)) {
    default:
	printf("rr_integrate: Defaulting to Euler integration\n");
    case EULER:
	for(i=0; i<nsteps; i++) {
	    rr_euler( y, deriv_func, t, n, dt );
	    if(tbase->verbosity) printf("%f \t%f\n",t,y[0]);
	    t=t0+dt*(i+1);
	}
    break;
    case RK2:
	for(i=0; i<nsteps; i++) {
	    rr_rk2( y, deriv_func, t, n, dt );
	    if(tbase->verbosity) printf("%f \t%f\n",t,y[0]);
	    t=t0+dt*(i+1);
	}

    break;
    case RK4:
	for(i=0; i<nsteps; i++) {
	    rr_rk4( y, deriv_func, t, n, dt );
	    if(tbase->verbosity) printf("%f \t%f\n",t,y[0]);
	    t=t0+dt*(i+1);
	}

    break;
    case AB2:
	ypold=malloc(sizeof(double)*n);

	//initialize with an euler step
	deriv_func(t,ym,n,ypold);
	rr_euler( ym, deriv_func, t, n, dt );
	for(j=0; j<n; j++) {
	    y[1]=ym[j];
	}
	t=t0+dt;
	//now jump into adams-bashforth integration
	for(i=1; i<nsteps; i++) {
	    rr_ab2( ym, deriv_func, t, n, dt, ypold );
	    for(j=0; j<n; j++) {
		y[i+1]=ym[j];
	    }
	    t=t0+dt*(i+1);
	}
	free(ypold);
    break;
    }
}
