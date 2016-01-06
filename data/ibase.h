#ifndef __IBASE_H__
#define __IBASE_H__
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

typedef struct {
    char *name;
    unsigned int n;
    double t0, t1, dt;
    double *t;
    unsigned int verbosity;
} rr_ibase;

/*Allocate t by providing t0, t1, and n*/
rr_ibase* rr_ibase_alloc_t(double pt0, double pt1, int pn);

/*Allocate t by providing t0, dt, and n*/
rr_ibase* rr_ibase_alloc_dt(double pt0, double pdt, double pn);

/*Specify a label*/
void rr_ibase_label(rr_ibase* ibase, const char* name);

/*List the time base, element by element*/
void rr_ibase_recite(rr_ibase* ibase);

/*Free allocated time base*/
void rr_ibase_free(rr_ibase* pt);

/*Set the verbosity threshold*/
void rr_ibase_set_verbosity(rr_ibase* ibase, int verbose);
void rr_ibase_populate(rr_ibase* ibase);
void rr_ibase_populate_with_ratio(rr_ibase* ibase, double ratio);
void rr_ibase_populate_with_lin_ratio(rr_ibase* ibase, double ratio1, double ratio2);

#endif
