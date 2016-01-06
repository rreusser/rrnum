#ifndef __STENCIL_H__
#define __STENCIL_H__

#include <rrdata.h>

void rr_choose_points(int npoints, rr_ibase* ibase, int i, int* ip, double* dp);

void rr_d2y_dx2(int npoints, double *xp, double* an);
void rr_dy_dx(int npoints, double* xp, double *an);

#endif
