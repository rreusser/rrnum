#include "bicg.h"

int assert(int i, int n) {
    if(i<0 || i > n-1) {
	fprintf(stderr,"ERROR i=%i, n=%i\n",i,n);
    }
    return i;
}
double inner_product(double* a, double* b, int n) {
    int i;
    double sum=0;
    for(i=0; i<n; i++) {
	sum+=a[i]*b[i];
    }
    return sum;
}


rr_sparse_mat* create_ilu_factorization( rr_sparse_mat* A ) {
    rr_sparse_mat* LU;
    LU = rr_sparse_mat_alloc( A->dim );
    rr_sparse_mat_copy( LU, A );
    ilu_factorization( LU );
    return LU;
}

void ilu_factorization(rr_sparse_mat* A) {
    int r;
    if(!(rr_sparse_mat_is_col_indexed(A))) {
	rr_sparse_mat_col_index(A);
    }
    int n=A->dim;
    double d,e;
    int rind,cind, row, col,index, i, j, k;
    for(r=0;r<n;r++) {
	if(A->dat[r]==0.0) {
	    fprintf(stderr,"ILU: failed due to zero pivot element.\n");
	    break;
	} else {
	    d=1.0/A->dat[r];
	    rind=A->cind[r];
	    for(i=rind;i<A->cind[r+1];i++) {
		row=A->cind[i];
		if(row > r) {
		    e = d*A->dat[A->cdat[i]];
		    //assert(A->cdat[i],A->size);
		    A->dat[A->cdat[i]] = e;
		    cind=A->cind[row];
		    for(j=cind;j<A->ind[row+1];j++) {
			col=A->ind[j];
			if(col>r) {
			    index=-1;
			    for(k=A->ind[r];k<A->ind[r+1];k++) {
				if(A->ind[k]==col) {
				    index = k;
				    break;
				}
			    }
			    if(index!=-1) {
				//assert(j,A->size);
				A->dat[j] -= e*A->dat[index];
			    }
			}
		    }
		    for(k=A->ind[r];k<A->ind[r+1];k++) {
			if(A->ind[k]==row) {
			    //assert(row,A->size);
			    A->dat[row] -= e*A->dat[k];
			}
		    }
		}
	    }
	}
    }
}

void ilu_preconditioner(rr_sparse_mat* LU, double* x, double* b) {
    register int i,j,k,j1,j2;
    double* dat = LU->dat;
    int* ind = LU->ind;
    int n = LU->dim;
    double* y;
    y = malloc(sizeof(double)*n);
    //forward substitution
    y[0] = b[0];
    memcpy(y,b,sizeof(double)*n);
    for(i=1;i<n;i++) {
	j1 = ind[ind[i]];
	j2 = ind[ind[i+1]];
	for(k=ind[i];k<ind[i+1];k++) {
	    j=ind[k];
	    if(j>i-1) {
		break;
	    } else {
		y[i] -= dat[k]*y[j];
	    }
	}
    }
    //back substitution
    memcpy(x,y,sizeof(double)*n);
    x[n-1] /= dat[n-1];
    for(i=n-2;i>=0;i--) {
	for(k=ind[i];k<ind[i+1];k++) {
	    j=ind[k];
	    if(j>i) {
		x[i] -= dat[k]*x[j];
	    }
	}
	x[i] /= dat[i];
    }
    free(y);
}



void bicg(rr_sparse_mat* A, double* x, double* b, rr_sparse_mat* LU, double tol, int itermax) {
    fprintf(stderr,"beginning allocation.\n");
    int i;
    int n = A->dim;
    double *p, *pt, *z, *zt, *q, *qt, *r, *rt;
    double rho1=0.0, rho2=0.0, beta=0.0, err=0.0, alpha=0.0;
    p  = malloc(sizeof(double)*n);
    pt = malloc(sizeof(double)*n);
    z  = malloc(sizeof(double)*n);
    zt = malloc(sizeof(double)*n);
    q  = malloc(sizeof(double)*n);
    qt = malloc(sizeof(double)*n);
    r  = malloc(sizeof(double)*n);
    rt = malloc(sizeof(double)*n);

    if(A->permutation!=NULL) {
	for(i=0;i<n;i++) {
	    p[i]=b[A->permutation[i]];
	}
	memcpy(b,p,sizeof(double)*n);
    }

    ilu_preconditioner(LU,x,b);

    //r=b-A*x
    vector_mult(A,x,r);
    for(i=0;i<n;i++) {
	r[i]=b[i]-r[i];
	rt[i]=r[i];
    }
    

    fprintf(stderr,"beginning iterations.\n");

    int converged=0;
    int iter=0;
    while(iter++<=itermax && converged==0) {
	//ilu_preconditioner(LU,z,r);
	//ilu_preconditioner(LU_T,zt,rt);
	//rho1=inner_product(z,rt,n);
	rho1=inner_product(r,rt,n);
	if(rho1==0) {
	    converged = 2;
	}
	if(iter==1) {
	    for(i=0;i<n;i++) {
		p[i]=r[i];
		pt[i]=rt[i];
		//p[i]=z[i];
		//pt[i]=zt[i];
	    }
	} else {
	    beta=rho1/rho2;
	    for(i=0;i<n;i++) {
		p[i]=r[i]+beta*p[i];
		pt[i]=rt[i]+beta*pt[i];
		//p[i]=z[i]+beta*p[i];
		//pt[i]=zt[i]+beta*pt[i];
	    }
	}
	vector_mult(A,p,q);
	vector_multT(A,pt,qt);
	alpha = rho1/inner_product(pt,q,n);
	err = 0.0;
	for(i=0; i<n; i++) {
	    x[i] = x[i] + alpha*p[i];
	    r[i] = r[i] - alpha*q[i];
	    rt[i] = rt[i] - alpha*qt[i];
	    err += r[i]*r[i];
	}
	rho2=rho1;
	err=sqrt(err);
	//if(iter%10==0) {
	    fprintf(stderr,"%cIteration %i: Residual = %e",13,iter,err);
	//}
	if(err*0 != 0) {
	    converged = 3;
	    break;
	}
	if(err<tol) {
	    converged = 1;
	    break;
	}
    }
    free(p);
    free(pt);
    free(z);
    free(zt);
    free(q);
    free(qt);
    free(r);
    free(rt);

    if(converged==1) {
	fprintf(stderr,"%cbi-cg: converged after %i iterations.\n",13,iter);
    } else if(converged == 0) {
	fprintf(stderr,"%cbi-cg: max iterations reached after %i iterations.\n",13,itermax);
    } else if(converged == 2) {
	fprintf(stderr,"%cbi-cg: failed to converge after %i iterations.\n",13,iter);
    } else if(converged == 3) {
	fprintf(stderr,"%cbi-cg: error: nan encountered after %i iterations.\n",13,iter);
    }
}



void bicgstab(rr_sparse_mat* A, double* x, double* b, rr_sparse_mat* LU, double tol, int itermax) {
    //double resid[itermax];
    double* resid = malloc(sizeof(double)*itermax);
    int i;
    int n = A->dim;
    double *p, *s, *t, *v, *r, *shat, *phat, *rtilde;
    double rho1=0.0, rho2=0.0, alpha=0.0, beta=0.0, omega=0.0, err=0.0, err0=0.0;
    p = malloc(sizeof(double)*n);
    s = malloc(sizeof(double)*n);
    t = malloc(sizeof(double)*n);
    v = malloc(sizeof(double)*n);
    r = malloc(sizeof(double)*n);
    shat = malloc(sizeof(double)*n);
    phat = malloc(sizeof(double)*n);
    rtilde= malloc(sizeof(double)*n);

    if(A->permutation!=NULL) {
	for(i=0;i<n;i++) {
	    p[i]=b[A->permutation[i]];
	}
	memcpy(b,p,sizeof(double)*n);
    }

    ilu_preconditioner(LU,x,b);

    //r = b - A*x
    vector_mult(A,x,r);
    for(i=0; i<n; i++) {
	r[i]=b[i]-r[i];
    }
    memcpy(rtilde,r,sizeof(double)*n);

    int converged = 0;
    int iter=0;
    while(iter++<=itermax && converged==0) {
	rho1 = inner_product(rtilde,r,n);
	if(iter==1) {
	    memcpy(p,r,sizeof(double)*n);
	} else {
	    beta=(rho1/rho2)*(alpha/omega);
	    for(i=0;i<n;i++) p[i]=r[i]+beta*(p[i]-omega*v[i]);
	}
	ilu_preconditioner(LU,phat,p);
	vector_mult(A,phat,v); //vector_mult(A,p,v);
	alpha=rho1/inner_product(rtilde,v,n);
	for(i=0;i<n;i++) s[i]=r[i]-alpha*v[i];
	ilu_preconditioner(LU,shat,s);
	vector_mult(A,shat,t); //vector_mult(A,s,t);
	omega=inner_product(t,s,n)/inner_product(t,t,n);
	for(i=0;i<n;i++) x[i]+=alpha*phat[i]+omega*shat[i];
	for(i=0;i<n;i++) r[i]=s[i]-omega*t[i]; //x[i]+=alpha*p[i]+omega*s[i];
	for(err=0.0,i=0;i<n;i++) err+=r[i]*r[i];
	err=sqrt(err);
	if(iter==1) err0=err;
	resid[iter]=err/err0;
	rho2=rho1;
	//if(iter%10==0 || iter==1) {
	    fprintf(stderr,"%cIteration %i: Residual = %e",13,iter,err/err0);
	//}
	if(isnan(err)) {
	    converged = 3;
	    break;
	}
	if(err/err0<tol) {
	    converged = 1;
	    break;
	}
    }
    if(converged==1) {
	//fprintf(stderr,"iter = %i,  err = %e\n",iter-1,err);
	fprintf(stderr,"%cbi-cgstab: converged after %i iterations to residual = %e\n",13,iter,err/err0);
	rr_ibase* iters;
	iters = rr_ibase_alloc_t(1.0,(double)iter,iter);
	rr_tecplotArray("residual.plt",resid,"residual",iters);
	rr_ibase_free(iters);
    } else if(converged == 0) {
	fprintf(stderr,"%cbi-cgstab: max iterations reached after %i iterations.\n",13,itermax);
    } else if(converged == 2) {
	fprintf(stderr,"%cbi-cgstab: failed to converge after %i iterations.\n",13,iter);
    } else if(converged == 3) {
	fprintf(stderr,"%cbi-cgstab: error: nan encountered after %i iterations.\n",13,iter);
    }
    free(resid);
    free(p);
    free(s);
    free(t);
    free(v);
    free(r);
    free(shat);
    free(phat);
    free(rtilde);
}
