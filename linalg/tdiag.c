#include "tdiag.h"

// Tridiagonal solver: expects A in the form of three n-length vectors that define
// [ b1 c1          ]
// [ a2 b2 c2       ]
// [    a3 b3 c3    ]
// [      .  .  .   ]
// [        .  .  . ]
// [          an bn ]
void rr_tridiag(double** A, double* b, const int n) {
    int i=0;
    double fac;
    //Eliminate
    for(i=n-2; i>=0; i--) {
        fac = A[0][i]/A[1][i+1];
        A[1][i] = A[1][i]-fac*A[0][i+1];
        b[i] = b[i]-fac*b[i+1];
    }
    //Back-substitute
    b[0]=b[0]/A[1][0];
    for(i=1;i<n;i++) {
        b[i] = (b[i]-A[0][i]*b[i-1])/A[1][i];
    }
}

void rr_pentadiag(double** A, double* b, const int n) {
    int i;
    double fac;
    //Eliminate 1
    for(i=2;i<n;i++) {
        fac=A[0][i]/A[1][i-1];
        if(A[1][i-1]==0) { fprintf(stderr,"Pentadiag failed due to zero pivot element.  May I suggest you get out and walk?\n"); exit(1); }
        A[1][i] -= fac*A[2][i-1];
        A[2][i] -= fac*A[3][i-1];
        A[3][i] -= fac*A[4][i-1];
        b[i] -= fac*b[i-1];
    }
    //Eliminate 2
    for(i=1;i<n;i++) {
        fac=A[1][i]/A[2][i-1];
        if(A[2][i-1]==0) { fprintf(stderr,"Pentadiag failed on elimination #2.  It's gonna be a long night.\n"); exit(1); }
        A[2][i] -= fac*A[3][i-1];
        A[3][i] -= fac*A[4][i-1];
        b[i] -= fac*b[i-1];
    }

    //Back-substitute
    b[n-1]/=A[2][n-1];
    b[n-2]=(b[n-2]-A[3][n-1]*b[n-1])/A[2][n-1];
    for(i=n-3;i>=0;i--) {
        b[i] = (b[i]-A[3][i]*b[i+1]-A[4][i]*b[i+2])/A[2][i];
    }
}

// INEFFICIENT, but it works...
// Expects input as three n-length vectors in A with:
//     [ b1 c1       a1 ]
//     [ a2 b2 c2       ]
// A = [    a3 b3 c3    ]
//     [      .  .  .   ]
//     [        .  .  . ]
//     [ cn       an bn ]
void rr_tridiag_periodic(double** A, double* b, const int n, double* work1, double* work2) {
    int i;
    double fac,A10;
    A10 = A[1][0];
    memcpy(work1, b, sizeof(double)*n);
    memset(work2, 0, sizeof(double)*n);
    work2[0] = A[1][0];
    work2[n-1] = A[2][n-1];
    A[1][n-1] -= A[0][n-1]*A[0][0]/A[1][0];
    A[1][0] = 0.0;
    //Eliminate
    for(i=n-2; i>=0; i--) {
        fac = A[0][i]/A[1][i+1];
        A[1][i] = A[1][i]-fac*A[0][i+1];
        work1[i] = work1[i]-fac*work1[i+1];
        work2[i] = work2[i]-fac*work2[i+1];
    }
    //Back-substitute
    work1[0]=work1[0]/A[1][0];
    work2[0]=work2[0]/A[1][0];
    for(i=1;i<n;i++) {
        work1[i] = (work1[i]-A[0][i]*work1[i-1])/A[1][i];
        work2[i] = (work2[i]-A[0][i]*work2[i-1])/A[1][i];
    }
    memcpy(b, work1, sizeof(double)*n);
    fac = (A10*work1[0]+A[0][0]*work1[n-1]) / (A10*(1.0+work2[0])+A[0][0]*work2[n-1]);
    for(i=0;i<n;i++) {
        b[i] = work1[i] - fac*work2[i];
    }
}
    
