#include "drawcontours.h"

void rrDrawContours(double* T, rr_structured_grid* grid, double cmin, double cmax, int num, int skip) {
    int nx = grid->isize;
    int ny = grid->jsize;

    int i,j, hit1, hit2, hit3, c;
    double x00,x10,x11,x01,y00,y10,y11,y01,T00,T10,T11,T01;
    double min1,max1,interp1,interp2,cstart,cend,f1=0.0,f2=0.0,f3=0.0;
    double start=0.0, cont;
    glBegin( GL_LINES );
    int skipx=skip,skipy=skip;
    int repeati=1, repeatj=1;
    i=0; j=0;
    //for(i=0; i<nx-skipx; i+=skipx) {
	//for(j=0; j<ny-skipy; j+=skipy) {
    while(repeati) {
	skipy=skip;
	repeatj=1;
	j=0;
	while(repeatj) {
	    x00 = grid->x[i+nx*j];
	    x10 = grid->x[i+skipx+nx*j];
	    x11 = grid->x[i+skipx+nx*(j+skipy)];
	    x01 = grid->x[i+nx*(j+skipy)];
	    y00 = grid->y[i+nx*j];
	    y10 = grid->y[i+skipx+nx*j];
	    y11 = grid->y[i+skipx+nx*(j+skipy)];
	    y01 = grid->y[i+nx*(j+skipy)];
	    T00 = T[i+nx*j];
	    T10 = T[i+skipx+nx*j];
	    T11 = T[i+skipx+nx*(j+skipy)];
	    T01 = T[i+nx*(j+skipy)];
	    min1 = dmin3(T00,T10,T01);
	    max1 = dmax3(T00,T10,T01);
	    interp1 = (max1-cmin)/(cmax-cmin);
	    interp2 = (max1-cmin)/(cmax-cmin);
	    cstart = dfloor(interp1*(double)(num-1));
	    cstart = (cstart<0)?(0):(cstart);
	    cend = dfloor(interp2*(double)(num-1));
	    cend = (cend>num-1)?(num-1):(cend);
	    for(c=start; c<cend+1; c++) {
		hit1=hit2=hit3=0;
		cont=cmin+(cmax-cmin)*c/(num-1);
		if(T10-T00!=0) {
		    f1=(cont-T00)/(T10-T00);
		    if(f1>0.0&&f1<1.0) { hit1++; }
		}
		if(T10-T01!=0) {
		    f2=(cont-T10)/(T01-T10);
		    if(f2>0.0&&f2<1.0) { hit2++; }
		}
		if(T01-T00!=0) {
		    f3=(cont-T01)/(T00-T01);
		    if(f3>0.0&&f3<1.0) { hit3++; }
		}
		if(hit1+hit2+hit3==2) {
		    if(hit1) glVertex2f(x00*(1.0-f1)+x10*f1,y00*(1.0-f1)+y10*f1);
		    if(hit2) glVertex2f(x10*(1.0-f2)+x01*f2,y10*(1.0-f2)+y01*f2);
		    if(hit3) glVertex2f(x01*(1.0-f3)+x00*f3,y01*(1.0-f3)+y00*f3);
		}
	    }
		//upper left triangle
	    min1 = dmin3(T10,T11,T01);
	    max1 = dmax3(T10,T11,T01);
	    interp1 = (max1-cmin)/(cmax-cmin);
	    interp2 = (max1-cmin)/(cmax-cmin);
	    cstart = dfloor(interp1*(double)(num-1));
	    cstart = (cstart<0)?(0):(cstart);
	    cend = dfloor(interp2*(double)(num-1));
	    cend = (cend>num-1)?(num-1):(cend);
	    for(c=start; c<cend+1; c++) {
		hit1=hit2=hit3=0;
		cont=cmin+(cmax-cmin)*c/(num-1);
		if(T10-T11!=0) {
		    f1=(cont-T10)/(T11-T10);
		    if(f1>0.0&&f1<1.0) { hit1++; }
		}
		if(T11-T01!=0) {
		    f2=(cont-T11)/(T01-T11);
		    if(f2>0.0&&f2<1.0) { hit2++; }
		}
		if(T01-T10!=0) {
		    f3=(cont-T01)/(T10-T01);
		    if(f3>0.0&&f3<1.0) { hit3++; }
		}
		if(hit1+hit2+hit3==2) {
		    if(hit1) glVertex2f(x10*(1.0-f1)+x11*f1,y10*(1.0-f1)+y11*f1);
		    if(hit2) glVertex2f(x11*(1.0-f2)+x01*f2,y11*(1.0-f2)+y01*f2);
		    if(hit3) glVertex2f(x01*(1.0-f3)+x10*f3,y01*(1.0-f3)+y10*f3);
		}
	    }
	    j+=skipy;
	    if(j>ny-2) repeatj=0;
	    if(j+skipy>ny-1) {
		skipy = (ny-1) - j;
	    }
	}
	i+=skipx;
	if(i>nx-2) repeati=0;
	if(i+skipx>nx-1) {
	    skipx = (nx-1) - i;
	}
    }
    glEnd();
}

