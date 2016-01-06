#include "tecio.h"

void rr_tecplotXYFunc(char *fname, double (*xfuncPtr)(double), double (*yfuncPtr)(double), char *mname, const rr_ibase* ibase) {
    int i;
    FILE *fp = NULL;
    double t;
    fp = fopen(fname,"w");
    
    fprintf(fp,"VARIABLES =\n\"t\"\n\"y\"\n");
    fprintf(fp,"Zone T=\"%s\", I=%i, F=point\n",mname,ibase->n);
    for(i=0; i<ibase->n+1; i++) {
	t=ibase->t[i];
	fprintf(fp,"%e\t%e\t\n",xfuncPtr(t),yfuncPtr(t));
    }

    fclose(fp);
}
void rr_tecplotConformalMap(char *fname, double (*xfuncPtr)(double,double), double (*yfuncPtr)(double,double), const rr_ibase* xbase, const rr_ibase* ybase, char *mname) {
    int i;
    FILE *fp = NULL;
    double x,y;
    fp = fopen(fname,"w");
    fprintf(fp,"VARIABLES =\n\"x\"\n\"y\"\n");
    fprintf(fp,"Zone T=\"%s\", I=%i, F=point\n",mname,xbase->n);
    for(i=0; i<xbase->n+1; i++) {
	x=xbase->t[i];
	y=ybase->t[i];
	fprintf(fp,"%e\t%e\t\n",xfuncPtr(x,y),yfuncPtr(x,y));
    }
    fclose(fp);
}

void rr_tecplotFunc(char *fname, double (*funcPtr)(double), char *mname, const rr_ibase* ibase) {
    int i;
    FILE *fp = NULL;
    double t;
    fp = fopen(fname,"w");
    
    fprintf(fp,"VARIABLES =\n\"t\"\n\"y\"\n");
    fprintf(fp,"Zone T=\"%s\", I=%i, F=point\n",mname,ibase->n);
    for(i=0; i<ibase->n+1; i++) {
	t=ibase->t[i];
	fprintf(fp,"%e\t%e\t\n",t,funcPtr(t));
    }

    fclose(fp);
}

void rr_tecplotXYArray(char *fname, double *x, double *y, int nsteps) {
    int i;
    FILE *fp = NULL;
    fp = fopen(fname,"w");
    
    fprintf(fp,"VARIABLES =\n\"x\"\n\"y\"\n");
    fprintf(fp,"Zone T=\"Array\", I=%i, F=point\n",nsteps);
    for(i=0; i<nsteps; i++) {
	fprintf(fp,"%e\t%e\t\n",x[i],y[i]);
    }
    fclose(fp);
}

void rr_tecplotArray(char *fname, double *y, char *mname, const rr_ibase* ibase) {
    int i;
    FILE *fp = NULL;
    fp = fopen(fname,"w");
    
    fprintf(fp,"VARIABLES =\n\"t\"\n\"y\"\n");
    fprintf(fp,"Zone T=\"%s\", I=%i, F=point\n",mname,ibase->n);
    for(i=0; i<ibase->n; i++) {
	fprintf(fp,"%e\t%e\t\n",ibase->t[i],y[i]);
    }
    fclose(fp);
}



void rr_tecplotODE(char *fname, double **y, char *mname, const rr_ibase* ibase, int n) {
    int i,j;
    FILE *fp = NULL;
    fp = fopen(fname,"w");
    fprintf(fp,"VARIABLES =\n\"t\"\n");
    for(i=0; i<n; i++) {
	fprintf(fp,"\"y%i\"\n",i);
    }
    /*stable up to y9... oh well... for now...*/
    //for(i=0; i<n; i++) {
	//fprintf(fp,"Zone T=\"%s%c\", I=%i, F=point\n",mname,(char)(48+i),ibase->n);
	//for(j=0; j<ibase->n; j++) {
	    //fprintf(fp,"%e\t%e\t\n",ibase->t[j],y[i][j]);
	//}
    //}

    fprintf(fp,"Zone T=\"%s%c\", I=%i, F=point\n",mname,(char)(48+i),ibase->n);
    for(j=0; j<ibase->n; j++) {
	fprintf(fp,"%e\t",ibase->t[j]);
	for(i=0; i<n; i++) {
	    fprintf(fp,"%e\t",y[i][j]);
	}
	fprintf(fp,"\n");
    }

    fclose(fp);
}


void rr_tecplotStructuredGrid(char *fname, rr_structured_grid* grid, double* T) {
    int i,j;
    FILE *fp = NULL;
    fp = fopen(fname,"w");
    fprintf(fp,"VARIABLES =\n\"x\"\n\"y\"\n\"f\"\n");
    /*stable up to y9... oh well... for now...*/
    fprintf(fp,"Zone T=\"Grid\", I=%i, J=%i, F=point\n", grid->jsize, grid->isize);
    double data;
    for(i=0; i<grid->isize; i++) {
	for(j=0; j<grid->jsize; j++) {
	    data = T[i+grid->isize*j];
	    if(isnan(data)==1) {data=0.0;}
	    fprintf(fp,"%e %e %e\n",grid->x[i+grid->isize*j],grid->y[i+grid->isize*j],data);
	}
    }
    fclose(fp);
}

void rr_tecplotPolar(char *fname, double (*rFunc)(double), double (*thetaFunc)(double), const rr_ibase* tbase) {
    int i;
    FILE *fp = NULL;
    double t, r, theta, x, y;
    fp = fopen(fname,"w");
    
    fprintf(fp,"VARIABLES =\n\"x\"\n\"y\"\n");
    fprintf(fp,"Zone T=\"PolarPlot\", I=%i, F=point\n",tbase->n);
    for(i=0; i<tbase->n+1; i++) {
	t=tbase->t[i];
	r = rFunc(t);
	theta = thetaFunc(t);
	x = r*cos(theta);
	y = r*sin(theta);
	fprintf(fp,"%e\t%e\t\n",x,y);
    }

    fclose(fp);
}
