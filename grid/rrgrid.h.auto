#ifndef __STRUCTURED_H__
#define __STRUCTURED_H__

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <rrdata.h>

#define NORMAL 0
#define H_TYPE 0
#define O_TYPE 1
#define C_TYPE 2

typedef struct {
    unsigned int type;
    unsigned int isize,jsize;
    double *x, *y;
} rr_structured_grid;

/*Allocate a structured grid*/
rr_structured_grid* rr_alloc_structured_grid(unsigned int ni, unsigned int nj, unsigned int type);

/*Build grid with trans-finite interpolation*/
void rr_transfinite_interpolation(rr_structured_grid* grid, rr_ibase* xbase, rr_ibase* ybase);

/*Apply elliptic smoothing to the grid*/
double rr_elliptic_smoother(rr_structured_grid* grid, unsigned int iterations, double underrelaxation);

/*March the grid out from j=0*/
void rr_hyperbolic_grid(rr_structured_grid* grid, double speed);

/*Free the allocated grid*/
void rr_structured_grid_free(rr_structured_grid* grid);

/*i-row grid definition for x-coords via function*/
void rr_i_row_gridx(rr_structured_grid* grid, double (*funcPtr)(double), rr_ibase* ibase, unsigned int row);

/*j-row grid definition for x-coords via function*/
void rr_j_row_gridx(rr_structured_grid* grid, double (*funcPtr)(double), rr_ibase* ibase, unsigned int row);

/*i-row grid definition for y-coords via function*/
void rr_i_row_gridy(rr_structured_grid* grid, double (*funcPtr)(double), rr_ibase* ibase, unsigned int row);

/*j-row grid definition for y-coords via function*/
void rr_j_row_gridy(rr_structured_grid* grid, double (*funcPtr)(double), rr_ibase* ibase, unsigned int row);

#endif
