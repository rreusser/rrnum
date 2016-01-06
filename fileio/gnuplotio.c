#include "gnuplotio.h"


void rr_gnuplotXYFunc(char *fname, double (*xfuncPtr)(double), double (*yfuncPtr)(double), const rr_ibase* ibase) {
    int i;
    FILE *fp = NULL;
    double t;
    fp = fopen(fname,"w");
    
    for(i=0; i<ibase->n+1; i++) {
	t=ibase->t[i];
	fprintf(fp,"%e\t%e\t\n",xfuncPtr(t),yfuncPtr(t));
    }

    fclose(fp);
}

void rr_gnuplotFunc(char *fname, double (*funcPtr)(double), const rr_ibase* ibase) {
    int i;
    FILE *fp = NULL;
    double t;
    fp = fopen(fname,"w");
    
    //fprintf(fp,"VARIABLES =\n\"t\"\n\"y\"\n");
    //fprintf(fp,"Zone T=\"%s\", I=%i, F=point\n",mname,ibase->n);
    for(i=0; i<ibase->n; i++) {
	t=ibase->t[i];
	fprintf(fp,"%e\t%e\t\n",t,funcPtr(t));
    }

    fclose(fp);
}

void rr_gnuplotXYArray(char *fname, double *x, double *y, int nsteps) {
    int i;
    FILE *fp = NULL;
    fp = fopen(fname,"w");
    
    //fprintf(fp,"VARIABLES =\n\"x\"\n\"y\"\n");
    //fprintf(fp,"Zone T=\"Array\", I=%i, F=point\n",nsteps);
    for(i=0; i<nsteps; i++) {
	fprintf(fp,"%e\t%e\t\n",x[i],y[i]);
    }
    fclose(fp);
}

void rr_gnuplotODE(char *fname, double **y, char *mname, const rr_ibase* ibase, int n) {
    int i,j;
    FILE *fp = NULL;
    fp = fopen(fname,"w");
    //fprintf(fp,"VARIABLES =\n\"t\"\n");
    //for(i=0; i<n; i++) {
	//fprintf(fp,"\"y%i\"\n",i);
    //}
    /*stable up to y9... oh well... for now...*/
    //for(i=0; i<n; i++) {
	//fprintf(fp,"Zone T=\"%s%c\", I=%i, F=point\n",mname,(char)(48+i),ibase->n);
	//for(j=0; j<ibase->n; j++) {
	    //fprintf(fp,"%e\t%e\t\n",ibase->t[j],y[i][j]);
	//}
    //}

    //fprintf(fp,"Zone T=\"%s%c\", I=%i, F=point\n",mname,(char)(48+i),ibase->n);
    for(j=0; j<ibase->n; j++) {
	fprintf(fp,"%e\t",ibase->t[j]);
	for(i=0; i<n; i++) {
	    fprintf(fp,"%e\t",y[i][j]);
	}
	fprintf(fp,"\n");
    }

    fclose(fp);
}

