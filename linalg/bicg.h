#ifndef __BICG_H__
#define __BICG_H__

#include <stdlib.h>
#include <stdio.h>
#include <rrgeneral.h>
#include <rrdata.h>
#include <rrfileio.h>
#ifndef __SPARSE_H__
#include "sparse.h"
#endif

double inner_product(double* a, double* b, int n);
rr_sparse_mat* create_ilu_factorization( rr_sparse_mat* A );
void ilu_factorization(rr_sparse_mat* A);
void ilu_preconditioner(rr_sparse_mat* LU, double* x, double* b);
void bicg(rr_sparse_mat* A, double* x, double* b, rr_sparse_mat* LU, double tol, int itermax);
void bicgstab(rr_sparse_mat* A, double* x, double* b, rr_sparse_mat* LU, double tol, int itermax);
#endif
