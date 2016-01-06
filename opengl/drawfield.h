#ifndef __DRAWFIELD_H__
#define __DRAWFIELD_H__

#include <GL/glut.h>
#include <rrgrid.h>
#include "rrgeneral.h"

void rrDrawLine(double x1, double y1, double x2, double y2);

void rrDrawEdge(rr_structured_grid* grid, int lineskip, int loopEdgeOn);

void rrDrawGrid(rr_structured_grid* grid, int lineskip);

void rrDrawField(double* T, rr_structured_grid* grid, double Tmin2, double Tmax2, int skip);

void rrDrawVelocity(rr_structured_grid* grid, double* u, double* v, int skip, double length);

void rrDrawXY(double*x, double*y, int nx, int skip, double scalex, double scaley);

#endif
