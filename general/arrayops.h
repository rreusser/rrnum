#ifndef __ARRAYOPS_H__
#define __ARRAYOPS_H__
void rr_fill(double* arr, double c, int size);
void rr_ones(double* arr, int size);
void rr_zeros(double* arr, int size);
void rr_subtract(double* A, double* B, double* C, int size);
void rr_add(double* A, double* B, double* C, int size);
void rr_minusEquals(double* A, double* B, int size);
void rr_plusEquals(double* A, double* B, int size);
#endif
