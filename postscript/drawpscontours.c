#include "drawpscontours.h"

void rrOpenPSFile(rrPSWindow* win, char* fname) {
    win->psfile = fopen(fname, "wt");
    if(win->psfile==NULL) {
	printf("Error opening file %s for output.\n",fname);
    } else {
	fprintf(win->psfile, "%%!PS\n");
	win->curtype=-1;
	win->state=-1;
    }
    win->xlist = malloc(4*sizeof(double));
    win->ylist = malloc(4*sizeof(double));
} 
void rrPSSetLineWidth(rrPSWindow *pswin, double width) {
    fprintf(pswin->psfile, "%f setlinewidth\n",width*pswin->glwin->yran2);
}
void rrPSCreateWindow(rrPSWindow *pswin, rr_window2d* glwin) {
    pswin->glwin = glwin;
}
void rrPSSetOffset(rrPSWindow* pswin, double x, double y) {
    pswin->pageoffset_x = x;
    pswin->pageoffset_y = y;
}
void rrPSBeginSection(rrPSWindow* pswin) {
    double width = 3.5;
    double height = width * pswin->glwin->yran2/pswin->glwin->xran2;
    fprintf(pswin->psfile, "gsave\n");
    fprintf(pswin->psfile, "72 72 scale\n");
    fprintf(pswin->psfile, "%f %f translate\n",4.25+pswin->pageoffset_x, 5.5+pswin->pageoffset_y);
    fprintf(pswin->psfile, "%f %f %f %f rectclip\n",-width,-height,2.0*width,2.0*height);
    fprintf(pswin->psfile, "%f %f scale\n",1.0/pswin->glwin->xran2*3.5,1.0/pswin->glwin->xran2*3.5);
    fprintf(pswin->psfile, "%f %f translate\n",-pswin->glwin->xcen,-pswin->glwin->ycen);
}
void rrPSEndSection(rrPSWindow* pswin) {
    fprintf(pswin->psfile,"grestore\n");
}

void rrClosePSFile(rrPSWindow* win) {
    fprintf(win->psfile,"showpage\n");
    fclose(win->psfile);
}
void rrPSBegin(rrPSWindow* win, int type) {
    win->curtype = type;
    if(win->curtype==RRPS_LINES) {
	fprintf(win->psfile,"newpath\n");
    }
    win->state=0;
    win->voidcur=0;
}
void rrPSEnd(rrPSWindow* win) {
    if(win->curtype==RRPS_LINES) {
	fprintf(win->psfile,"stroke\n");
    } else if (win->curtype==RRPS_TRIANGLES) {
	//fprintf(win->psfile,"fill\n");
    } else if (win->curtype==RRPS_QUADS) {
	//fprintf(win->psfile,"fill\n");
    }
}
int rrIsOutsideWindow(rr_window2d* win, double x, double y,double fudge) {
    double fx = fudge*win->xran2*2.0;
    double fy = fudge*win->yran2*2.0;
    return (x<win->xmin-fx || x>win->xmax+fx || y<win->ymin-fy || y>win->ymax+fy);
}

int rrIsInWindow(rr_window2d* win, double x, double y,double fudge) {
    double fx = fudge*win->xran2*2.0;
    double fy = fudge*win->yran2*2.0;
    return (x>win->xmin-fx && x<win->xmax+fx && y>win->ymin-fy && y<win->ymax+fy);
}
double rrFitInWindowY(rr_window2d* win, double y) {
    return (y-win->ymin)/(2.0*win->yran2);
}
double rrFitInWindowX(rr_window2d* win, double x) {
    return (x-win->xmin)/(2.0*win->xran2);
}

void rrPSSetColor(rrPSWindow *pswin, double r, double g, double b) {
    fprintf(pswin->psfile,"%f %f %f setrgbcolor\n",r,g,b);
}

void rrPSSetRotation(rrPSWindow* win, double rotdeg) {
    win->rot=rotdeg*M_PI/180.0;
}
void rrPSVertex(rrPSWindow* win, double x, double y) {
    double rot=win->rot;
    int i;
    if(win->psfile!=NULL) {
	if(win->state==0) {
	    win->xlist[win->state] = x;
	    win->ylist[win->state] = y;
	    win->voidcur += rrIsInWindow(win->glwin,x,y,0.1);
	} else if (win->state < win->curtype-1) {
	    win->xlist[win->state] = x;
	    win->ylist[win->state] = y;
	    win->voidcur += rrIsInWindow(win->glwin,x,y,0.1);
	} else {
	    if(win->voidcur>0) {
		if(win->curtype>2) fprintf(win->psfile,"newpath\n");
		fprintf(win->psfile, "%f %f moveto\n",win->xlist[0]*cos(rot)+win->ylist[0]*sin(rot),-win->xlist[0]*sin(rot)+win->ylist[0]*cos(rot));
		for(i=1; i<win->curtype-1; i++) {
		    fprintf(win->psfile, "%f %f lineto\n",win->xlist[i]*cos(rot)+win->ylist[i]*sin(rot),-win->xlist[i]*sin(rot)+win->ylist[i]*cos(rot));
		}
		fprintf(win->psfile, "%f %f lineto\n",x*cos(rot)+y*sin(rot),-x*sin(rot)+y*cos(rot));
		if(win->curtype>2) fprintf(win->psfile,"fill\n");
	    }
	    win->voidcur=0;
	}
	win->state = (win->state+1)%(win->curtype);
    } else {
	printf("Error writing to file.\n");
    }
}



void rrPSContours(rrPSWindow* win, double* T, rr_structured_grid* grid, double cmin, double cmax, int num, int skip) {
    int nx = grid->isize;
    int ny = grid->jsize;

    int i,j, hit1, hit2, hit3, c;
    double x00,x10,x11,x01,y00,y10,y11,y01,T00,T10,T11,T01;
    double min1,max1,interp1,interp2,cstart,cend,f1=0.0,f2=0.0,f3=0.0;
    double start=0.0, cont;
    int skipx=skip,skipy=skip;
    int repeati=1, repeatj=1;
    i=0; j=0;
    rrPSBegin(win, RRPS_LINES );
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
		    if(hit1) rrPSVertex(win,x00*(1.0-f1)+x10*f1,y00*(1.0-f1)+y10*f1);
		    if(hit2) rrPSVertex(win,x10*(1.0-f2)+x01*f2,y10*(1.0-f2)+y01*f2);
		    if(hit3) rrPSVertex(win,x01*(1.0-f3)+x00*f3,y01*(1.0-f3)+y00*f3);
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
		    if(hit1) rrPSVertex(win,x10*(1.0-f1)+x11*f1,y10*(1.0-f1)+y11*f1);
		    if(hit2) rrPSVertex(win,x11*(1.0-f2)+x01*f2,y11*(1.0-f2)+y01*f2);
		    if(hit3) rrPSVertex(win,x01*(1.0-f3)+x10*f3,y01*(1.0-f3)+y10*f3);
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
    rrPSEnd(win);
}

void rrPSDrawEdge(rrPSWindow* pswin, rr_structured_grid* grid, int skip, int loopEdgeOn) {
    int nx = grid->isize;
    int ny = grid->jsize;
    rrPSBegin(pswin, RRPS_LINES );
    int skipy=skip, skipx=skip;
    int i=0, j=0;
    int repeati=1, repeatj=1;
    //for(i=0; i<nx-1; i++) {
    if(grid->type!=O_TYPE || loopEdgeOn==1) {
	while(repeati) {
	    rrPSVertex(pswin,grid->x[i], grid->y[i]);
	    rrPSVertex(pswin,grid->x[i+skipx], grid->y[i+skipx]);
	    rrPSVertex(pswin,grid->x[i+nx*(ny-1)], grid->y[i+nx*(ny-1)]);
	    rrPSVertex(pswin,grid->x[i+skipx+nx*(ny-1)], grid->y[i+skipx+nx*(ny-1)]);
	    i+=skipx;
	    if(i>nx-2) repeati=0;
	    if(i+skipx>nx-1) {
		skipx = (nx-1) - i;
	    }
	}
    }
    //for(j=0; j<ny-1; j++) {
    while(repeatj) {
	rrPSVertex(pswin,grid->x[nx*j], grid->y[nx*j]);
	rrPSVertex(pswin,grid->x[nx*(j+skipy)], grid->y[nx*(j+skipy)]);
	rrPSVertex(pswin,grid->x[nx-1+nx*j], grid->y[nx-1+nx*j]);
	rrPSVertex(pswin,grid->x[nx-1+nx*(j+skipy)], grid->y[nx-1+nx*(j+skipy)]);
	j+=skipy;
	if(j>ny-2) repeatj=0;
	if(j+skipy>ny-1) {
	    skipy = (ny-1) - j;
	}
    }
    rrPSEnd(pswin);
}

void rrPSDrawVelocity(rrPSWindow* pswin, rr_structured_grid* grid, double* u, double* v, int skip, double length) {
    double vmag,up,vp,xp,yp;
    int nx = grid->isize;
    int ny = grid->jsize;
    int i,j;
    rrPSBegin(pswin, RRPS_LINES );
    for(i=0; i<nx; i+=skip) {
	for(j=0; j<ny; j+=skip) {
	    xp=grid->x[i+nx*j];
	    yp=grid->y[i+nx*j];
	    up=u[i+nx*j];
	    vp=v[i+nx*j];
	    vmag=sqrt(up*up+vp*vp);
	    //up/=vmag;
	    //vp/=vmag;
	    rrPSVertex(pswin,xp,yp);
	    rrPSVertex(pswin,xp+up*length,yp+vp*length);
	}
    }
    rrPSEnd(pswin);
}
void rrPSDrawGrid(rrPSWindow* pswin,rr_structured_grid* grid, int skip) {
    int nx = grid->isize;
    int ny = grid->jsize;
    int i,j;
    int skipx=skip,skipy=skip;
    i=0; j=0;
    int repeati=1, repeatj=1;
    rrPSBegin(pswin, RRPS_LINES );
    while(repeati) {
	skipy=skip;
	repeatj=1;
	j=0;
	while(repeatj) {
	    rrPSVertex(pswin,grid->x[i+nx*j], grid->y[i+nx*j]);
	    rrPSVertex(pswin,grid->x[i+skipx+nx*j], grid->y[i+skipx+nx*j]);
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
    skipx=skip,skipy=skip;
    repeati=1, repeatj=1;
    i=0; j=0;
    while(repeatj) {
	i=0;
	skipx=skip;
	repeati=1;
	while(repeati) {
	    rrPSVertex(pswin,grid->x[i+nx*j], grid->y[i+nx*j]);
	    rrPSVertex(pswin,grid->x[i+nx*(j+skipy)], grid->y[i+nx*(j+skipy)]);
	    i+=skipx;
	    if(i>nx-2) repeati=0;
	    if(i+skipx>nx-1) {
		skipx = (nx-1) - i;
	    }
	}
	j+=skipy;
	if(j>ny-2) repeatj=0;
	if(j+skipy>ny-1) {
	    skipy = (ny-1) - j;
	}
    }
    rrPSEnd(pswin);
}

void rrPSDrawField(rrPSWindow* pswin, double* T, rr_structured_grid* grid, double Tmin2, double Tmax2, int skip) {
    int i,j;
    int nx=grid->isize;
    int ny=grid->jsize;
    double Tp;
    rrPSBegin(pswin, RRPS_QUADS );
    int skipx=skip,skipy=skip;
    int repeati = 1, repeatj;
    i=0;
    double r1,r2,r3,r4;
    double g1,g2,g3,g4;
    double b1,b2,b3,b4;
    int ind1,ind2,ind3,ind4;
    while(repeati) {
	skipy=skip;
	repeatj=1;
	j=0;
	while(repeatj) {
	    ind1=i+nx*j;
	    Tp=(T[ind1]-Tmin2)/(Tmax2-Tmin2);
	    hsv(240-240*Tp,1.0,1.0,&r1,&g1,&b1);

	    ind2=i+skipx+nx*j;
	    Tp=(T[ind2]-Tmin2)/(Tmax2-Tmin2);
	    hsv(240-240*Tp,1.0,1.0,&r2,&g2,&b2);

	    ind3=i+skipx+nx*(j+skipy);
	    Tp=(T[ind3]-Tmin2)/(Tmax2-Tmin2);
	    hsv(240-240*Tp,1.0,1.0,&r3,&g3,&b3);

	    ind4=i+nx*(j+skipy);
	    Tp=(T[ind4]-Tmin2)/(Tmax2-Tmin2);
	    hsv(240-240*Tp,1.0,1.0,&r4,&g4,&b4);
	    
	    rrPSSetColor(pswin,0.25*(r1+r2+r3+r4),0.25*(g1+g2+g3+g4),0.25*(b1+b2+b3+b4));
	    rrPSVertex(pswin,grid->x[ind1], grid->y[ind1]);
	    rrPSVertex(pswin,grid->x[ind2], grid->y[ind2]);
	    rrPSVertex(pswin,grid->x[ind3], grid->y[ind3]);
	    rrPSVertex(pswin,grid->x[ind4], grid->y[ind4]);
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
    rrPSEnd(pswin);
}


void rrPSDrawXY(rrPSWindow* pswin, double*x, double*y, int nx, int skip, double scalex, double scaley) {
    rrPSBegin(pswin,RRPS_LINES );
    int repeati = 1;
    int i=0;
    while(repeati) {
	rrPSVertex(pswin,x[i]*scalex,y[i]*scaley);
	rrPSVertex(pswin,x[i+skip]*scalex,y[i+skip]*scaley);
	i+=skip;
	if(i>nx-2) repeati=0;
	if(i+skip>nx-1) {
	    skip = (nx-1) - i;
	}
    }
    rrPSEnd(pswin);
}



