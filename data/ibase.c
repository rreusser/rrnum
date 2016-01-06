#include "ibase.h"


rr_ibase* rr_ibase_alloc_t(double pt0, double pt1, int pn) {
    int i;
    rr_ibase* ibase;
    ibase = malloc(sizeof(rr_ibase));
    ibase->t0 = pt0;
    ibase->t1 = pt1;
    ibase->n = pn;
    ibase->dt = (pt1-pt0)/(pn-1);
    double *temp;
    temp = malloc(sizeof(double)*pn);
    ibase->t = temp;
    for(i=0; i<pn; i++) {
	ibase->t[i]=pt0 + (pt1-pt0)*i/(pn-1);
    }
    ibase->verbosity=0;
    ibase->name = malloc(sizeof(char)*129);
    return ibase;
}

rr_ibase* rr_ibase_alloc_dt(double pt0, double pdt, double pn) {
    int i;
    rr_ibase* ibase;
    ibase = malloc(sizeof(rr_ibase));
    ibase->t0 = pt0;
    ibase->n = pn;
    ibase->dt = pdt;
    ibase->t1 = pdt*(pn-1);
    double *temp;
    temp = malloc(sizeof(double)*pn);
    ibase->t = temp;
    for(i=0; i<pn; i++) {
	ibase->t[i]=pt0 + (ibase->t1-pt0)*i/(pn-1);
    }
    ibase->verbosity=0;
    ibase->name = malloc(sizeof(char)*128);
    return ibase;
}

void rr_ibase_label(rr_ibase* ibase, const char* name) {
    strncpy(ibase->name,name,128);
}

void rr_ibase_set_verbosity(rr_ibase* ibase,int verbose) {
    ibase->verbosity=verbose;
}

void rr_ibase_populate(rr_ibase* ibase) {
    int i;
    for(i=0; i<ibase->n; i++) {
	ibase->t[i]=ibase->t0 + (ibase->t1-ibase->t0)*i/(ibase->n-1);
    }
}

void rr_ibase_populate_with_ratio(rr_ibase* ibase, double ratio) {
    int i;
    double max;
    double inc = 1.0;
    ibase->t[0]=0.0;
    for(i=1; i<ibase->n; i++) {
	ibase->t[i]=ibase->t[i-1]+inc;
	inc*=ratio;
    }
    max=ibase->t[ibase->n-1];
    for(i=0; i<ibase->n; i++) {
	ibase->t[i]=ibase->t[i]/max*(ibase->t1-ibase->t0)+ibase->t0;
    }
}

void rr_ibase_populate_with_lin_ratio(rr_ibase* ibase, double ratio1, double ratio2) {
    int i;
    double max, interp;
    double inc = 1.0;
    ibase->t[0]=0.0;
    for(i=1; i<ibase->n; i++) {
	interp = (double)i/(ibase->n-1);
	ibase->t[i]=ibase->t[i-1]+inc;
	inc*=(1.0-interp)*ratio1 + interp*ratio2;
    }
    max=ibase->t[ibase->n-1];
    for(i=0; i<ibase->n; i++) {
	ibase->t[i]=ibase->t[i]/max*(ibase->t1-ibase->t0)+ibase->t0;
    }
}


void rr_ibase_recite(rr_ibase* ibase) {
    if(ibase->verbosity) {
	printf("dt = \t%f\n",ibase->dt);
	printf("n  = \t%i\n",ibase->n);
	printf("t0 = \t%f\n",ibase->t0);
	printf("t1 = \t%f\n",ibase->t1);
	int i;
	for(i=0; i<ibase->n; i++) {
	    printf("%i\t%f\n",i,ibase->t[i]);
	    fflush(stdout);
	}
    }
}
void rr_ibase_free(rr_ibase* ibase) {
    free(ibase->name);
    free(ibase->t);
    free(ibase);
}
