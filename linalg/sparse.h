#ifndef __SPARSE_H__
#define __SPARSE_H__

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include <rrgeneral.h>


struct node;

typedef struct node{
    struct node* next;
    double data;
    int index;
} node;

typedef struct rr_sparse_mat {
    int size;
    int dim;
    int* ind;
    int* length;
    int* permutation;
    double* dat;
    double* diagonal;
    double* scale;
    node** row;
    int* cind;
    int* cdat;
    int col_indexed;
} rr_sparse_mat;

rr_sparse_mat* rr_sparse_mat_alloc(int n);
void rr_sparse_mat_dealloc(rr_sparse_mat* A);
int rr_sparse_mat_is_col_indexed(rr_sparse_mat* A);
void rr_sparse_mat_copy(rr_sparse_mat* B, rr_sparse_mat *A);
void rr_sparse_mat_copyT(rr_sparse_mat* A, rr_sparse_mat *B);
void rr_sparse_mat_calculate_scale(rr_sparse_mat* A);
void rr_sparse_mat_partial_pivot( rr_sparse_mat* A );
void rr_sparse_mat_cleanup( rr_sparse_mat* A);
int rr_sparse_mat_dimension(rr_sparse_mat* A);
int rr_sparse_mat_dyn_size(rr_sparse_mat* A);
int rr_sparse_mat_alloc_size(rr_sparse_mat* A);
int rr_sparse_mat_nonzeros(rr_sparse_mat* A);
void rr_sparse_mat_unlink(rr_sparse_mat* A);
void rr_sparse_mat_output(rr_sparse_mat* A);
void rr_sparse_mat_insert(rr_sparse_mat* A, int i, int j, double data);
void rr_sparse_mat_build(rr_sparse_mat* A);
void vector_mult(rr_sparse_mat* A, double* x, double* b);
void vector_multT(rr_sparse_mat* A, double* x, double* b);
void rr_sparse_mat_col_output(rr_sparse_mat* A);
void rr_sparse_mat_ind_output(rr_sparse_mat* A);
void rr_sparse_mat_col_index(rr_sparse_mat* A);
void rr_sparse_mat_ppm_output(rr_sparse_mat* A, char fname[], float scale);
void stats(rr_sparse_mat* A);

double inner_product(double* a, double* b, int n);



int* get_cind(rr_sparse_mat* A, int i);
int* get_cdat(rr_sparse_mat* A, int i);
int* get_ind(rr_sparse_mat* A, int i);
double* get_dat(rr_sparse_mat* A, int i);
double* diag(rr_sparse_mat* A, int i);
void rr_sparse_mat_stats(rr_sparse_mat* A);

#endif
