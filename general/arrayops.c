#include "arrayops.h"

void rr_fill(double* arr, double c, int size) {
    int i=0;
    for(i=0; i<size; i++) {
	arr[i] = c;
    }
}
void rr_ones(double* arr, int size) {
    int i=0;
    for(i=0; i<size; i++) {
	arr[i] = 1.0;
    }
}
void rr_zeros(double* arr, int size) {
    int i=0;
    for(i=0; i<size; i++) {
	arr[i] = 0.0;
    }
}
void rr_subtract(double* A, double* B, double* C, int size) {
    int i=0;
    for(i=0; i<size; i++) {
	C[i]=A[i]-B[i];
    }
}
void rr_add(double* A, double* B, double* C, int size) {
    int i=0;
    for(i=0; i<size; i++) {
	C[i]=A[i]+B[i];
    }
}
void rr_minusEquals(double* A, double* B, int size) {
    int i=0; 
    for(i=0; i<size; i++) {
	A[i]-=B[i];
    }
}
void rr_plusEquals(double* A, double* B, int size) {
    int i=0; 
    for(i=0; i<size; i++) {
	A[i]+=B[i];
    }
}
