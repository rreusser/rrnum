#include "structured.h"

/*Allocate a structured grid*/
rr_structured_grid* rr_alloc_structured_grid(unsigned int ni, unsigned int nj, unsigned int type) {
    rr_structured_grid* grid;
    grid = malloc(sizeof(rr_structured_grid));
    grid->isize = ni;
    grid->jsize = nj;
    grid->type = type;
    grid->x = malloc(sizeof(double)*ni*nj);
    grid->y = malloc(sizeof(double)*ni*nj);
    return grid;
}

/*Build grid with trans-finite interpolation*/
void rr_transfinite_interpolation(rr_structured_grid* grid, rr_ibase* xbase, rr_ibase* ybase) {
    int i,j,ni,nj;
    double xi,eta,*x,*y;
    x=grid->x;
    y=grid->y;
    ni=grid->isize;
    nj=grid->jsize;
    for(i=1; i<ni-1; i++) {
        //xi=(double)i/(ni-1);
        xi = xbase->t[i];
        for(j=1; j<nj-1; j++) {
            eta = ybase->t[j];
            //eta=(double)j/(nj-1);
            x[i+ni*j] =
                (1.0-xi)*x[ni*j] +
                xi*x[ni-1+ni*j] +
                (1.0-eta)*x[i] +
                eta*x[i+ni*(nj-1)] -
                (1.0-xi)*(1.0-eta)*x[0] -
                xi*(1.0-eta)*x[ni-1] -
                (1.0-xi)*eta*x[(nj-1)*ni] -
                xi*eta*x[ni-1+(nj-1)*ni];
            y[i+ni*j] =
                (1.0-xi)*y[ni*j] +
                xi*y[ni-1+ni*j] +
                (1.0-eta)*y[i] +
                eta*y[i+ni*(nj-1)] -
                (1.0-xi)*(1.0-eta)*y[0] -
                xi*(1.0-eta)*y[ni-1] -
                (1.0-xi)*eta*y[(nj-1)*ni] -
                xi*eta*y[ni-1+(nj-1)*ni];
        }
    }
}

/*Apply elliptic smoothing to the grid*/
double rr_elliptic_smoother(rr_structured_grid* grid, unsigned int iterations, double underrelaxation) {
    int i,j,ni,nj,iter, i0,i1,j0,j1,is,js, in,ip,im, jp,jm,jn,ot;
    double *x,*y,ax,ay,bx,by,g11,g12,g22,errx,erry,err=0.0,xtemp,ytemp;
    x=grid->x;
    y=grid->y;
    ni=grid->isize;
    nj=grid->jsize;
    ot=grid->type==O_TYPE;
    printf("Elliptic smoother with Gauss-Seidel iterations:\n");
    i0=0; i1=1; is=1;
    j0=0; j1=1; js=1;
    double underone = 1.0-underrelaxation;
    for(iter=0; iter<iterations; iter++) {
        err=0.0;
        //Alternate directions for gauss-seidel iteration
        for(i=1+(ni-3)*i0; i>0 && i<ni-1; i+=is) {
            for(j=1-ot+(nj-3+ot)*j0; j>0-ot && j<nj-1; j+=js) {
                //printf("i,j = %i,%i\n",i,j);
                if(i==0) { im=ni-2;      } else { im=i-1;      }
                if(j==0) { jm=(nj-2)*ni; } else { jm=(j-1)*ni; }
                in=i;
                ip=i+1;
                jn=j*ni;
                jp=(j+1)*ni;

                ax = x[ip+jn]-x[im+jn];
                ay = y[ip+jn]-y[im+jn];
                bx = x[i+jp]-x[i+jm];
                by = y[i+jp]-y[i+jm];

                g11=(ax*ax+ay*ay)*0.25;
                g22=(bx*bx+by*by)*0.25;
                g12=((x[ip+jn]-x[im+jn])*(x[i+jp]-x[i+jm])+
                    (y[ip+jn]-y[im+jn])*(y[i+jp]-y[i+jm]))*0.25;

                xtemp = 0.5/(g11+g22)*(
                    g22*x[ip+jn] + g11*x[in+jp] + g11*x[in+jm] + g22*x[im+jn]
                    + 0.5*(g12*x[ip+jm]+g12*x[im+jp] - g12*x[ip+jp] - g12*x[im+jm]));

                ytemp = 0.5/(g11+g22)*(
                    g22*y[ip+jn] + g11*y[in+jp] + g11*y[in+jm] + g22*y[im+jn]
                    + 0.5*(g12*y[ip+jm]+g12*y[im+jp] - g12*y[ip+jp] - g12*y[im+jm]));

                errx = x[i+jn]-xtemp;
                erry = y[i+jn]-ytemp;
                err = err + errx*errx + erry*erry;

                x[i+jn] = x[i+jn]*underone + underrelaxation*xtemp;
                y[i+jn] = y[i+jn]*underone + underrelaxation*ytemp;
                if(j==0) {
                    x[in+ni*(nj-1)] = x[in+jn];
                    y[in+ni*(nj-1)] = y[in+jn];
                }
                if(i==0) {
                    x[ni-1+jn] = x[in+jn];
                    y[ni-1+jn] = y[in+jn];
                }
            }
        }
        if((iter+1)%(iterations/10)==0 || iter==1) {
            printf("iter=%i,  err=%f\n",iter+1,sqrt( err/(double)((nj-2)*(ni-2)) ));
        }
        i0=!i0; i1=!i1; is*=-1;
        j0=!j0; j1=!j1; js*=-1;
    }
    return sqrt(err);
}

/*March the grid out from j=0*/
void rr_hyperbolic_grid(rr_structured_grid* grid, double speed) {

}

/*Free the allocated grid*/
void rr_structured_grid_free(rr_structured_grid* grid) {
    free(grid->x);
    free(grid->y);
    free(grid);
    grid=NULL;
}

/*i-row grid definition for x-coords via function*/
void rr_i_row_gridx(rr_structured_grid* grid, double (*funcPtr)(double), rr_ibase* ibase, unsigned int row) {
    int i;
    int ni=grid->isize;
    double x;
    for(i=0; i<ni; i++) {
        x = ibase->t[i];
        //x = (double)i/(ni-1);
        grid->x[i+ni*row] = funcPtr(x);
    }
}

/*j-row grid definition for x-coords via function*/
void rr_j_row_gridx(rr_structured_grid* grid, double (*funcPtr)(double), rr_ibase* ibase, unsigned int row) {
    int j;
    int ni=grid->isize;
    double y;
    //printf("size of ibase = %i\n",ibase->n);
    for(j=0; j<grid->jsize; j++) {
        y = ibase->t[j];
        //printf("j=%i,\ty=%f\n",j,y);
        //y = (double)j/(grid->jsize-1);
        grid->x[row+ni*j] = funcPtr(y);
    }
}

/*i-row grid definition for y-coords via function*/
void rr_i_row_gridy(rr_structured_grid* grid, double (*funcPtr)(double), rr_ibase* ibase, unsigned int row) {
    int i;
    int ni=grid->isize;
    double x;
    for(i=0; i<ni; i++) {
        x = ibase->t[i];
        //x = (double)i/(ni-1);
        grid->y[i+ni*row] = funcPtr(x);
    }
}

/*j-row grid definition for y-coords via function*/
void rr_j_row_gridy(rr_structured_grid* grid, double (*funcPtr)(double), rr_ibase* ibase, unsigned int row) {
    int j;
    int ni=grid->isize;
    double y;
    for(j=0; j<grid->jsize; j++) {
        y = ibase->t[j];
        //y = (double)j/(grid->jsize-1);
        grid->y[row+ni*j] = funcPtr(y);
    }
}
