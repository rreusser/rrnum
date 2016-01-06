#include "color.h"

void hsv(double h, double s, double v, double* r, double* g, double* b) {
    int hi = ((int)(h/60.0))%6;
    double f=h/60.0-hi;
    double p=v*(1.0-s);
    double q=v*(1.0-f*s);
    double t=v*(1.0-(1.0-f)*s);
    switch(hi) {
	case 0: *r=v; *g=t; *b=p; break;
	case 1: *r=q; *g=v; *b=p; break;
	case 2: *r=p; *g=v; *b=t; break;
	case 3: *r=p; *g=q; *b=v; break;
	case 4: *r=t; *g=p; *b=v; break;
	case 5: *r=v; *g=p; *b=q; break;
    }
}
