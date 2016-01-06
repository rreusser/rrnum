#include "drawfield.h"

void rrDrawEdge(rr_structured_grid* grid, int skip, int loopEdgeOn) {
    int nx = grid->isize;
    int ny = grid->jsize;
    glBegin( GL_LINES );
    int skipy=skip, skipx=skip;
    int i=0, j=0;
    int repeati=1, repeatj=1;
    //for(i=0; i<nx-1; i++) {
    if(grid->type!=O_TYPE || loopEdgeOn==1) {
	while(repeati) {
	    glVertex2f(grid->x[i], grid->y[i]);
	    glVertex2f(grid->x[i+skipx], grid->y[i+skipx]);
	    glVertex2f(grid->x[i+nx*(ny-1)], grid->y[i+nx*(ny-1)]);
	    glVertex2f(grid->x[i+skipx+nx*(ny-1)], grid->y[i+skipx+nx*(ny-1)]);
	    i+=skipx;
	    if(i>nx-2) repeati=0;
	    if(i+skipx>nx-1) {
		skipx = (nx-1) - i;
	    }
	}
    }
    //for(j=0; j<ny-1; j++) {
    while(repeatj) {
	glVertex2f(grid->x[nx*j], grid->y[nx*j]);
	glVertex2f(grid->x[nx*(j+skipy)], grid->y[nx*(j+skipy)]);
	glVertex2f(grid->x[nx-1+nx*j], grid->y[nx-1+nx*j]);
	glVertex2f(grid->x[nx-1+nx*(j+skipy)], grid->y[nx-1+nx*(j+skipy)]);
	j+=skipy;
	if(j>ny-2) repeatj=0;
	if(j+skipy>ny-1) {
	    skipy = (ny-1) - j;
	}
    }
    glEnd();
}

void rrDrawVelocity(rr_structured_grid* grid, double* u, double* v, int skip, double length) {
    double vmag,up,vp,xp,yp;
    int nx = grid->isize;
    int ny = grid->jsize;
    int i,j;
    glBegin( GL_LINES );
    for(i=0; i<nx; i+=skip) {
	for(j=0; j<ny; j+=skip) {
	    xp=grid->x[i+nx*j];
	    yp=grid->y[i+nx*j];
	    up=u[i+nx*j];
	    vp=v[i+nx*j];
	    vmag=sqrt(up*up+vp*vp);
	    //up/=vmag;
	    //vp/=vmag;
	    glVertex2f(xp,yp);
	    glVertex2f(xp+up*length,yp+vp*length);
	}
    }
    glEnd();
}
void rrDrawGrid(rr_structured_grid* grid, int skip) {
    int nx = grid->isize;
    int ny = grid->jsize;
    int i,j;
    int skipx=skip,skipy=skip;
    i=0; j=0;
    int repeati=1, repeatj=1;
    glBegin( GL_LINES );
    while(repeati) {
	skipy=skip;
	repeatj=1;
	j=0;
	while(repeatj) {
	    glVertex2f(grid->x[i+nx*j], grid->y[i+nx*j]);
	    glVertex2f(grid->x[i+skipx+nx*j], grid->y[i+skipx+nx*j]);
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
	    glVertex2f(grid->x[i+nx*j], grid->y[i+nx*j]);
	    glVertex2f(grid->x[i+nx*(j+skipy)], grid->y[i+nx*(j+skipy)]);
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
    glEnd();
}

void rrDrawField(double* T, rr_structured_grid* grid, double Tmin2, double Tmax2, int skip) {
    int i,j;
    int nx=grid->isize;
    int ny=grid->jsize;
    double Tp, r,g,b;
    glBegin( GL_QUADS );
    int skipx=skip,skipy=skip;
    int ind;
    int repeati = 1, repeatj;
    i=0;
    while(repeati) {
	skipy=skip;
	repeatj=1;
	j=0;
	while(repeatj) {
	    ind=i+nx*j;
	    Tp=(T[ind]-Tmin2)/(Tmax2-Tmin2);
	    hsv(240-240*Tp,1.0,1.0,&r,&g,&b);
	    glColor3f(r,g,b);
	    glVertex2f(grid->x[ind], grid->y[ind]);

	    ind=i+skipx+nx*j;
	    Tp=(T[ind]-Tmin2)/(Tmax2-Tmin2);
	    hsv(240-240*Tp,1.0,1.0,&r,&g,&b);
	    glColor3f(r,g,b);
	    glVertex2f(grid->x[ind], grid->y[ind]);

	    ind=i+skipx+nx*(j+skipy);
	    Tp=(T[ind]-Tmin2)/(Tmax2-Tmin2);
	    hsv(240-240*Tp,1.0,1.0,&r,&g,&b);
	    glColor3f(r,g,b);
	    glVertex2f(grid->x[ind], grid->y[ind]);

	    ind=i+nx*(j+skipy);
	    Tp=(T[ind]-Tmin2)/(Tmax2-Tmin2);
	    hsv(240-240*Tp,1.0,1.0,&r,&g,&b);
	    glColor3f(r,g,b);
	    glVertex2f(grid->x[ind], grid->y[ind]);
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


void rrPSDrawXY(double*x, double*y, int nx, int skip, double scalex, double scaley) {
    rrpsBegin( GL_LINES );
    int repeati = 1;
    int i=0;
    while(repeati) {
	glVertex3f(x[i]*scalex,y[i]*scaley,0);
	glVertex3f(x[i+skip]*scalex,y[i+skip]*scaley,0);
	i+=skip;
	if(i>nx-2) repeati=0;
	if(i+skip>nx-1) {
	    skip = (nx-1) - i;
	}
    }
    rrpsEnd();
}



