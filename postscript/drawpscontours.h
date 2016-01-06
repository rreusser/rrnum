#ifndef __DRAWPSCONTOURS__H__
#define __DRAWPSCONTOURS__H__

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <rrgeneral.h>
#include <rrgrid.h>
#include <stdio.h>
#include <rropengl.h>

#define RRPS_POINTS 1
#define RRPS_LINES 2
#define RRPS_TRIANGLES 3
#define RRPS_QUADS 4

typedef struct {
    rr_window2d* glwin;
    FILE* psfile;
    int state;
    int curtype;
    double* xlist;
    double* ylist;
    int voidcur;
    double r,g,b;
    double pageoffset_y, pageoffset_x;
    double rot;
} rrPSWindow;

void rrPSSetRotation(rrPSWindow* win, double rotdeg);
void rrPSContours(rrPSWindow* win, double* T, rr_structured_grid* grid, double cmin, double cmax, int num, int skip);
void rrClosePSFile(rrPSWindow* win);
void rrOpenPSFile(rrPSWindow* win, char* fname);
void rrPSBeginSection(rrPSWindow* pswin);
void rrPSEndSection(rrPSWindow* pswin);

void rrPSSetOffset(rrPSWindow* pswin, double x, double y);
void rrPSVertex(rrPSWindow* win, double x, double y);
void rrPSBegin(rrPSWindow* win, int type);
void rrPSEnd(rrPSWindow* win);
void rrPSCreateWindow(rrPSWindow *pswin, rr_window2d* glwin);
int rrIsInWindow(rr_window2d* win, double x, double y, double fudge);
int rrIsOutsideWindow(rr_window2d* win, double x, double y, double fudge);
void rrPSSetColor(rrPSWindow *pswin, double r, double g, double b);
void rrPSDrawGrid(rrPSWindow* pswin,rr_structured_grid* grid, int skip);
void rrPSSetLineWidth(rrPSWindow *pswin, double width);
void rrPSDrawEdge(rrPSWindow* pswin, rr_structured_grid* grid, int skip, int loopEdgeOn);
void rrPSDrawField(rrPSWindow* pswin, double* T, rr_structured_grid* grid, double Tmin2, double Tmax2, int skip);
void rrPSDrawXY(rrPSWindow* pswin, double*x, double*y, int nx, int skip, double scalex, double scaley);

#endif
