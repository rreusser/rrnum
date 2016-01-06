#include <GL/glut.h>
#include <rrgrid.h>
#include "rrgeneral.h"

void rrDrawEdge(rr_structured_grid* grid, int lineskip, int loopEdgeOn);

void rrDrawGrid(rr_structured_grid* grid, int lineskip);

void rrDrawField(double* T, rr_structured_grid* grid, double Tmin2, double Tmax2, int skip);

void rrDrawVelocity(rr_structured_grid* grid, double* u, double* v, int skip, double length);
