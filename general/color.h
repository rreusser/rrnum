#ifndef __COLOR_H__
#define __COLOR_H__

void hsv(double h, double s, double v, double* r, double* g, double* b);

typedef struct {
    double r,g,b;
} rrColor3f;

typedef struct {
    double r,g,b,a;
} rrColor4f;

#endif
