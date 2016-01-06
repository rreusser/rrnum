#include "multigrid.h"
#include "rrgeneral.h"

#define poisson 4.0

void rrmgAlloc(int isize, int jsize, int maxlevels) {
    mg->isize = malloc(maxlevels*sizeof(int));
    mg->jsize = malloc(maxlevels*sizeof(int));
    double itest, jtest;
    mg->isize[0] = isize;
    mg->jsize[0] = jsize;
    mg->inc[0] = 1;
    for(i=1; i<maxlevels; i++) {
	mg->inc[i] = mg->inc[i-1]*2;
	mg->isize[i] = (mg->isize[i-1]-1)/2+1;
	itest = (double)(mg->isize[i-1]-1)/2.0+1;
	if((int)itest != mg->isize[i]) {
	    printf("Error: rrmgAlloc: isize not evenly divisible.\n");
	    continue;
	}
	mg->jsize[i] = (mg->jsize[i-1]-1)/2+1;
	jtest = (double)(mg->isize[i-1]-1)/2.0+1;
	if((int)jtest != mg->jsize[i]) {
	    printf("Error: rrmgAlloc: jsize not evenly divisible.\n");
	    continue;
	}
    }
    rrmgSolver* mg = malloc(sizeof(rrmgSolver));
    mg->A = malloc(maxlevels*sizeof(double*));
    mg->M = malloc(maxlevels*sizeof(double*));
    mg->r = malloc(maxlevels*sizeof(double*));
    mg->du = malloc(maxlevels*sizeof(double*));
    mg->pdu = malloc(maxlevels*sizeof(double*));
    mg->Ldu = malloc(maxlevels*sizeof(double*));
    int i;
    for(i=0; i<maxlevels; i++) {
	A[i] = malloc(isize[i]*jsize[i]*sizeof(double)*9);
	M[i] = malloc(isize[i]*jsize[i]*sizeof(double)*4);
	r[i] = malloc(isize[i]*jsize[i]*sizeof(double));
	du[i] = malloc(isize[i]*jsize[i]*sizeof(double));
	pdu[i] = malloc(isize[i]*jsize[i]*sizeof(double));
	Ldu[i] = malloc(isize[i]*jsize[i]*sizeof(double));
    }
    return mg;
}
void rrmgFree(rrmgSolver* mg) {
    int i;
    for(i=0; i<maxlevels; i++) {
	free(A[i]);
	free(M[i]);
	free(r[i]);
	free(du[i]);
	free(pdu[i]);
	free(Ldu[i]);
    }
    free(isize);
    free(jsize);
    free(A);
    free(M);
    free(r);
    free(du);
    free(pdu);
    free(Ldu);
}

void rrmgSetGrid(rrmgSolver* mg, rr_structured* grid) {
    mg->grid = grid;
}

void rrmgMetric() {
    double* y = mg->grid->y;
    double* x = mg->grid->x;
    int i,j,nx=mg->grid->isize, ny=mg->grid->jsize;
    double dydeta, dxdeta, dxdxi, dydxi, J;
    int i0,im,ip, j0,jm,jp, nxl, nyl;
    double dxi, deta;
    int otype = mg->grid->type==O_TYPE;
    int htype = mg->grid->type==H_TYPE;
    for(lev=0; lev<mg->maxlevels; lev++) {
	nxl = mg->isize[lev];
	nyl = mg->jsize[lev];
	for(i=0, il=0; i<nx; i+=inc[lev], il++) {
	    for(j=0, jl=0; j<ny; j+=inc[lev], jl++) {
		dxi = deta = 2.0;
		i0 = i
		im = i-inc[lev];
		ip = i+inc[lev];
		if(i==0) { im=0; dxi=1.0; }
		if(i==nx-inc[lev]) { ip=nx-inc[lev]; dxi=1.0; }
		j0 = nx*j;
		jp = nx*(j+inc[lev]);
		jm = nx*(j-inc[lev]);
		if(j==0 && htype) { jm=0; deta=1.0; }
		if(j==0 && otype) { jm=nx*(ny-2*inc[lev]); }
		if(j==ny-inc[lev] && htype) { jp=nx*(j+inc[lev]); deta=1.0; }
		if(j==ny-inc[lev] && otype) { jp=nx*inc[lev]; }
		dxdxi = (x[ip+j0] - x[im+j0])/dxi;
		dydxi = (y[ip+j0] - y[im+j0])/dxi;
		dxdeta = (x[i0+jp] - x[i0+jm])/deta;
		dydeta = (y[i0+jp] - y[i0+jm])/deta;
		J = dxidxi*dydeta - dxdeta*dydxi;
		M[(il+jl*nxl)*4] = dydeta/J;
		M[(il+jl*nxl)*4+1] = -dydxi/J;
		M[(il+jl*nxl)*4+2] = -dxdeta/J;
		M[(il+jl*nxl)*4+3] = dxdxi/J;
	    }
	}
    }
}
void rrmgLaplacianAMatrix() {
    int i,j;
    int nxl,nyl;
    for(lev=0; lev<mg->maxlevels; lev++) {
	nxl = mg->isize[lev];
	nyl = mg->jsize[lev];
	for(i=0; i<nxl; i++) {
	    for(j=0; j<nyl; j++) {
		i0=i;
		im=i-1;
		ip=i+1;
		j0=nxl*j;
		jp=nxl*(j+1);
		jm=nxl*(j-1);
		if(j==0 && otype==1) { jm=nxl*(nyl-2); }
		i0j04 = (i0+j0)*4;
		i0j09 = (i0+j0)*9;
		A[i0j09  ]  = (M[i0j04  ]*M[(im+j0)*4+1]+M[i0j04+1]*M[(i0+jm)*4  ])*0.25
			    + (M[i0j04+2]*M[(im+j0)*4+3]+M[i0j04+3]*M[(i0+jm)*4+2])*0.25;
		A[i0j09+1]  = M[i0j04+1]*(M[i0j04+1]+M[(i0+jm)*4+1])*0.5
			    + M[i0j04+3]*(M[i0j04+3]+M[(i0+jm)*4+3])*0.5;
		A[i0j09+2]  =-(M[i0j04  ]*M[(ip+j0)*4+1]+M[i0j04+1]*M[(i0+jm)*4  ])*0.25
			    - (M[i0j04+2]*M[(ip+j0)*4+3]+M[i0j04+3]*M[(i0+jm)*4+2])*0.25;
		A[i0j09+3]  = M[i0j04  ]*(M[i0j04  ]+M[(im+j0)*4  ])*0.5
			    + M[i0j04+2]*(M[i0j04+2]+M[(im+j0)*4+2])*0.5;
		A[i0j09+4]  = -M[i0j04  ]*(M[i0j04  ]+0.5*(M[(ip+j0)*4  ]+M[(im+j0)*4  ]))
			    -  M[i0j04+2]*(M[i0j04+2]+0.5*(M[(ip+j0)*4+2]+M[(im+j0)*4+2]))
			    -  M[i0j04+1]*(M[i0j04+1]+0.5*(M[(i0+jp)*4+1]+M[(i0+jm)*4+1]))
			    -  M[i0j04+3]*(M[i0j04+3]+0.5*(M[(i0+jp)*4+3]+M[(i0+jm)*4+3]));
		A[i0j09+5]  = M[i0j04  ]*(M[i0j04  ]+M[(ip+j0)*4  ])*0.5
			    + M[i0j04+2]*(M[i0j04+2]+M[(ip+j0)*4+2])*0.5;
		A[i0j09+6]  =-(M[i0j04  ]*M[(im+j0)*4+1]+M[i0j04+1]*M[(i0+jp)*4  ])*0.25
			    - (M[i0j04+2]*M[(im+j0)*4+3]+M[i0j04+3]*M[(i0+jp)*4+2])*0.25;
		A[i0j09+7]  = M[i0j04+1]*(M[i0j04+1]+M[(i0+jp)*4+1])*0.5
			    + M[i0j04+3]*(M[i0j04+3]+M[(i0+jp)*4+3])*0.5;
		A[i0j09+8]  = (M[i0j04  ]*M[(ip+j0)*4+1]+M[i0j04+1]*M[(i0+jp)*4  ])*0.25
			    + (M[i0j04+2]*M[(ip+j0)*4+3]+M[i0j04+3]*M[(i0+jp)*4+2])*0.25;
	    }
	}
    }
}

double mgGaussSeidelConst(double* x, double h, int n, int iters, double omega, double c) {
    int iter;
    int i,j;
    double h2=h*h;
    int size=(2<<n)+1;
    double temp,temp2,resid=10.0;
    double ph2 = poisson*h2;
    for(iter=0; iter<iters; iter++) {
	for(i=1; i<size-1; i++) {
	    for(j=1; j<size-1; j++) {
	       x[i+size*j]=0.25*(x[i+1+size*j]+x[i-1+size*j]+x[i+size*(j+1)]+x[i+size*(j-1)] - ph2);
	    }
	}
    }
    return resid;
}

double mgGaussSeidelZero(mgSolver* mg, double* x, int lev) {
   double *A=mg->A[lev];
   int i,j,i0,j0,im,jm,ip,jp,i0j09;
   int nx=mg->isize[lev];
   int ny=mg->jsize[lev];
   int otype = (mg->grid->type==O_TYPE);
   for(iter=0; iter<mg->lambda; iter++) {
       for(i=1; i<nx-1; i++) {
	   for(j=1-otype; j<ny-1; j++) {
		i0=i;
		ip=i+1;
		im=i-1;
		j0=nx*j;
		jp=nx*(j+1);
		jm=nx*(j-1);
		if(j==0){ jm=nx*(ny-2); }
		if(j==ny-2 && otype){ jp=0; }
		i0j09=(i0+j0)*9;
		x[i+nx*j]=(
		    x[im+jm]*A[i0j09  ] +
		    x[i0+jm]*A[i0j09+1] +
		    x[ip+jm]*A[i0j09+2] +
		    x[im+j0]*A[i0j09+3] +
		    x[ip+j0]*A[i0j09+5] +
		    x[im+jp]*A[i0j09+6] +
		    x[i0+jp]*A[i0j09+7] +
		    x[ip+jp]*A[i0j09+8]) /
		    A[i0j09+5];
	    }
	}
    }
}

/*double mgGaussSeidelZero(double* x, double h, int n, int iters, double omega) {
    int iter;
    int i,j;
    int size=(2<<n)+1;
    double h2 = h*h;
    double temp,temp2,resid=0.0;
    for(iter=0; iter<iters; iter++) {
	for(i=1; i<size-1; i++) {
	    for(j=1; j<size-1; j++) {
	       x[i+size*j]=0.25*(x[i+1+size*j]+x[i-1+size*j]+x[i+size*(j+1)]+x[i+size*(j-1)]);
	       //temp=0.25*(x[i+1+size*j]+x[i-1+size*j]+x[i+size*(j+1)]+x[i+size*(j-1)]);
	       //temp2=temp-x[i+size*j];
	       //resid+=temp2*temp2;
	       //x[i+size*j]+=temp2*omega;
	    }
	}
    }
    return resid;
}*/

double mgGaussSeidelNeg(mgSolver* mg, double* x, double* b, int lev) {
   double *A=mg->A[lev];
   int i,j,i0,j0,im,jm,ip,jp,i0j09;
   int nx=mg->isize[lev];
   int ny=mg->jsize[lev];
   int otype = (mg->grid->type==O_TYPE);
   for(iter=0; iter<mg->lambda; iter++) {
       for(i=1; i<nx-1; i++) {
	   for(j=1-otype; j<ny-1; j++) {
		i0=i;
		ip=i+1;
		im=i-1;
		j0=nx*j;
		jp=nx*(j+1);
		jm=nx*(j-1);
		if(j==0){ jm=nx*(ny-2); }
		if(j==ny-2 && otype){ jp=0; }
		i0j09=(i0+j0)*9;
		x[i+nx*j]=(
		    x[im+jm]*A[i0j09  ] +
		    x[i0+jm]*A[i0j09+1] +
		    x[ip+jm]*A[i0j09+2] +
		    x[im+j0]*A[i0j09+3] +
		    x[ip+j0]*A[i0j09+5] +
		    x[im+jp]*A[i0j09+6] +
		    x[i0+jp]*A[i0j09+7] +
		    x[ip+jp]*A[i0j09+8]) /
		    A[i0j09+5];
	    }
	}
    }
}

double mgGaussSeidelNeg(double* x, double* b, double h, int n, int iters, double omega) {
    int iter;
    int i,j;
    int size=(2<<n)+1;
    double h2=h*h;
    double temp,temp2,resid=0.0;
    for(iter=0; iter<iters; iter++) {
	for(i=1; i<size-1; i++) {
	    for(j=1; j<size-1; j++) {
	       //temp=0.25*(x[i+1+size*j]+x[i-1+size*j]+x[i+size*(j+1)]+x[i+size*(j-1)]+h2*b[i+size*j]);
	       //temp2=temp-x[i+size*j];
	       //resid+=temp2*temp2;
	       //x[i+size*j]+=temp2*omega;
	       x[i+size*j]=0.25*(x[i+1+size*j]+x[i-1+size*j]+x[i+size*(j+1)]+x[i+size*(j-1)]+h2*b[i+size*j]);
	    }
	}
    }
    return resid;
}

void mgResidp(double* x, double* b, double h, int n) {
    int i,j;
    int size=(2<<n)+1;
    double h2=h*h;
    int im,jp,jm,ip;
    for(i=1; i<size-1; i++) {
	for(j=1; j<size-1; j++) {
	    //im = (i==0)?(0):(i-1);
	    //jp = (j==size-1)?(j*size):(size*(j+1));
	    //jm = (j==0)?(0):(j-1)*size;
	    //ip = (i==size-1)?(i):(i+1);
	    //b[i+size*j]=0.25*(x[ip+size*j]+x[im+size*j]+x[i+jp]+x[i+jm]-4.0*x[i+size*j])/h2 - poisson;
	    b[i+size*j]=(x[i+1+size*j]+x[i-1+size*j]+x[i+size*(j+1)]+x[i+size*(j-1)]-4.0*x[i+size*j])/h2 - poisson;
	}
    }
}
void mgResid(double* x, double* b, double h, int n) {
    int i,j;
    int size=(2<<n)+1;
    double h2=h*h;
    int im,jp,jm,ip;
    for(i=1; i<size-1; i++) {
	for(j=1; j<size-1; j++) {
	    im = (i==0)?(0):(i-1);
	    jp = (j==size-1)?(j*size):(size*(j+1));
	    jm = (j==0)?(0):(j-1)*size;
	    ip = (i==size-1)?(i):(i+1);
	    b[i+size*j]=0.25*(x[ip+size*j]+x[im+size*j]+x[i+jp]+x[i+jm]-4.0*x[i+size*j])/h2;
	}
    }
}

double mgGaussSeidel(double* x, double* b, double h, int n, int iters, double omega) {
    int i,j;
    int size=(2<<n)+1;
    double h2=h*h;
    double temp,temp2,resid=0.0;
    for(i=1; i<size-1; i++) {
	for(j=1; j<size-1; j++) {
	   x[i+size*j]=0.25*(x[i+1+size*j]+x[i-1+size*j]+x[i+size*(j+1)]+x[i+size*(j-1)]-h2*b[i+size*j]);
	   //temp=0.25*(x[i+1+size*j]+x[i-1+size*j]+x[i+size*(j+1)]+x[i+size*(j-1)]-h2*b[i+size*j]);
	   //temp2=temp-x[i+size*j];
	   //resid+=temp2*temp2;
	   //x[i+size*j]+=temp2*omega;
	}
    }
    return resid;
}

void mgRestrict(double* rr, double* r, int n) {
    int i,j,ip,jp,s,sp;
    s=(2<<n)+1;
    sp=(2<<(n+1))+1;
    for(i=0,ip=0; i<s; i++,ip+=2) {
	for(j=0,jp=0; j<s; j++,jp+=2) {
	    rr[i+s*j]=r[ip+sp*jp];
	}
    }
}
void mgProlongate(double* pdu, double* du, int n) {
    int i,j,im,jm,s,sm;
    s=(2<<n)+1;
    sm=(2<<(n-1))+1;
    for(i=0,im=0; i<s; i+=2,im++) {
	for(j=0,jm=0; j<s; j+=2,jm++) {
	    pdu[i+s*j]=du[im+sm*jm];
	}
    }
    for(i=1,im=0; i<s; i+=2,im++) {
	for(j=0,jm=0; j<s; j+=2,jm++) {
	    pdu[i+s*j]=0.5*(du[im+sm*jm]+du[im+1+sm*jm]);
	}
    }
    for(i=0,im=0; i<s; i+=2,im++) {
	for(j=1,jm=0; j<s; j+=2,jm++) {
	    pdu[i+s*j]=0.5*(du[im+sm*jm]+du[im+sm*(jm+1)]);
	}
    }
    for(i=1; i<s; i+=2) {
	for(j=1; j<s; j+=2) {
	    pdu[i+s*j]=0.5*(pdu[i-1+s*j]+pdu[i+1+s*j]);
	}
    }
}

double mgLevel(double** du, double** pdu, double** r, double** Ldu, double h0, int n, double omega, int minlevel, int lambda) {
    int size=(2<<n)+1;
    double resid=1e10;
    mgRestrict(r[n],r[n+1],n);
    rr_zeros(du[n],size*size);
    mgGaussSeidelNeg(du[n],r[n],h0/size,n,lambda,omega);
    if(n>minlevel) {
	mgResid(du[n],Ldu[n],h0/size,n);
	rr_add(r[n],du[n],r[n],size*size);
	mgLevel(du,pdu,r,Ldu,h0,n-1,omega,minlevel,lambda);
	mgProlongate(pdu[n],du[n-1],n);
	rr_plusEquals(du[n],pdu[n],size*size);
	resid = mgGaussSeidelNeg(du[n],r[n],h0/size,n,lambda,omega);
    }
    return resid;
}

double mgL0Norm(double* r, int size) {
    int i;
    double lmax=0.0;
    for(i=0; i<size; i++) {
	lmax=max(dabs(r[i]),lmax);
    }
    return lmax;
}
double mgL2Norm(double* r, int size) {
    int i;
    double sum=0.0;
    for(i=0; i<size; i++) {
	sum+=r[i]*r[i];
    }
    return sqrt(sum);
}
double mgIter(double* Field, double** pdu, double** du, double** r, double** Ldu, double h0, int n, double tol, int itermax, double omega, int minlevel, int lambda) {
    double resid = tol+1e10;
    int l0 = n-1;
    int size=(2<<l0)+1;
    resid = mgGaussSeidelConst(Field, h0/size, l0, lambda/2, omega, poisson);
    mgResidp(Field,r[l0],h0/size,l0);
    resid = mgLevel(du,pdu,r,Ldu,h0,l0-1, omega,minlevel, lambda/2);
    mgProlongate(pdu[l0],du[l0-1],l0);
    rr_plusEquals(Field,pdu[l0],size*size);
    rr_zeros(r[l0],size*size);
    mgResidp(Field,r[l0],h0/size,l0);
    return mgL2Norm(r[l0],size*size);
}


