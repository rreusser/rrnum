#ifndef __WINDOW2D_H__
#define __WINDOW2D_H__

#include <GL/glut.h>
#include <rrgeneral.h>
#include <stdio.h>

typedef struct {
    int window;
    int fps;
    rrColor4f bgColor;
    void (*redisp)(int value);
    int width, height;
    int drawOn;
    double aspect_ratio;
    double xmin, xmax;
    double ymin, ymax;
    double xcen, ycen;   //center point
    double xran2, yran2; //half of the total x-range/y-range
    double xmin0,xmax0;
    double ymin0, ymax0;
    double xcen0, ycen0;   //center point
    double xran20, yran20; //half of the total x-range/y-range

    double backplane,frontplane;  //clipping planes
} rr_window2d;

void rr_window_copy_bounds(rr_window2d* win);
void rr_window_redisp_func(rr_window2d* win, void (*funcptr)(int value));
void rr_toggle_drawOn(rr_window2d* win);
void rr_init_glut_window2d_bounds(rr_window2d* win, int* argc, char **argv,const char* name, const int width, const int height, const int x0, const int y0);
void rr_init_glut_window2d(rr_window2d* win, int* argc, char **argv);
void rr_window_bounds2d(rr_window2d* win, double xmin, double xmax, double ymin, double ymax);
void rr_set_fps2d(rr_window2d* win, int value);
void rr_inc_fps2d(rr_window2d* win);
void rr_dec_fps2d(rr_window2d* win);
void rr_reshape2d(rr_window2d* win, int w, int h);
void rr_resquare2d(rr_window2d* win);
void rr_rewindow2d(rr_window2d* win);
void rr_zoom2d(rr_window2d* win, double factor);
void rr_pan2d(rr_window2d* win, double dx, double dy);
void rr_zoomlive2d(rr_window2d* win, double factor);

double rr_xmousepoint2d(rr_window2d* win, int x);
double rr_ymousepoint2d(rr_window2d* win, int y);

#endif
