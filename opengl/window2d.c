#include "window2d.h"

void rr_init_glut_window2d(rr_window2d* win, int* argc, char **argv) {
    win->drawOn=1;
    glutInit(argc, argv);
    glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE );
    glutInitWindowSize(500,500);
    glutInitWindowPosition(100,100);
    win->window = glutCreateWindow("Window");

    win->bgColor.r = 1.0;
    win->bgColor.g = 1.0;
    win->bgColor.b = 1.0;
    win->bgColor.a = 0.0;
    glClearColor(win->bgColor.r,
		 win->bgColor.g,
		 win->bgColor.b,
		 win->bgColor.a);
    
    glClearDepth( 1.0f );
    glEnable( GL_DEPTH_TEST );
    glDepthFunc( GL_LEQUAL );
    glEnable( GL_COLOR_MATERIAL );
}

void rr_init_glut_window2d_bounds(rr_window2d* win, int* argc, char **argv,const char* name, const int width, const int height, const int x0, const int y0) {
    win->drawOn=1;
    glutInit(argc, argv);
    glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE );
    glutInitWindowSize(width,height);
    glutInitWindowPosition(x0,y0);
    win->window = glutCreateWindow(name);

    win->bgColor.r = 1.0;
    win->bgColor.g = 1.0;
    win->bgColor.b = 1.0;
    win->bgColor.a = 0.0;
    glClearColor(win->bgColor.r,
		 win->bgColor.g,
		 win->bgColor.b,
		 win->bgColor.a);
    
    glClearDepth( 1.0f );
    glEnable( GL_DEPTH_TEST );
    glDepthFunc( GL_LEQUAL );
    glEnable( GL_COLOR_MATERIAL );
}

void rr_window_copy_bounds(rr_window2d* win) {
    win->xran20= win->xran2;
    win->yran20= win->yran2;
    win->ycen0= win->ycen;
    win->xcen0= win->xcen;
    win->xmin0= win->xmin;
    win->xmax0= win->xmax;
    win->ymin0= win->ymin;
    win->ymax0= win->ymax;
}
void rr_set_fps2d(rr_window2d* win, int value) {
    win->fps = value;
    if(value < 1) {
	printf("Error: rr_set_fps2d: invalid framerate %i, setting to 1fps.\n",value);
	win->fps=1;
    }
}

void rr_dec_fps2d(rr_window2d* win) {
    win->fps--;
    if(win->fps < 1) {
	printf("Error: rr_dec_fps2d: invalid framerate %i, setting to 1fps.\n",win->fps);
	win->fps=1;
    }
}
void rr_inc_fps2d(rr_window2d* win) {
    win->fps++;
    if(win->fps > 1000) {
	printf("Error: rr_dec_fps2d: invalid framerate %i, setting to 1000fps.\n",win->fps);
	win->fps=1000;
    }
}
void rr_toggle_drawOn(rr_window2d* win) {
    toggle(&win->drawOn);
}
    
void rr_window_redisp_func(rr_window2d* win, void (*funcptr)(int value)) {
    win->redisp = funcptr;
}

void rr_window_bounds2d(rr_window2d* win, double xmin, double xmax, double ymin, double ymax) {
    win->xmin = xmin;
    win->xmax = xmax;
    win->ymin = ymin;
    win->ymax = ymax;
    win->xcen  = (xmax+xmin)/2.0;
    win->ycen  = (ymax+ymin)/2.0;
    win->xran2 = (xmax-xmin)/2.0;
    win->yran2 = (ymax-ymin)/2.0;
}
void rr_reshape2d(rr_window2d* win, int w, int h) {
    win->width = w;
    win->height = h;
    rr_resquare2d( win );
    glViewport(0, 0, w, h);
    rr_rewindow2d( win );
}

void rr_resquare2d(rr_window2d* win) {
    win->aspect_ratio = (double)win->width/(double)win->height;
    win->xran2 = win->yran2*win->aspect_ratio;
    win->xmax = win->xcen+win->xran2;
    win->xmin = win->xcen-win->xran2;
    rr_rewindow2d( win );
}


void rr_rewindow2d(rr_window2d* win) {
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(win->xmin, win->xmax, win->ymin, win->ymax);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void rr_zoomlive2d(rr_window2d* win, double factor) {
    win->xran2 = win->xran20*factor;
    win->yran2 = win->yran20*factor;
    win->xmin = win->xcen0-win->xran2;
    win->xmax = win->xcen0+win->xran2;
    win->ymin = win->ycen0-win->yran2;
    win->ymax = win->ycen0+win->yran2;
    rr_rewindow2d(win);
    win->redisp(0);
}
void rr_zoom2d(rr_window2d* win, double factor) {
    win->xran2 *= factor;
    win->yran2 *= factor;
    win->xmin = win->xcen-win->xran2;
    win->xmax = win->xcen+win->xran2;
    win->ymin = win->ycen-win->yran2;
    win->ymax = win->ycen+win->yran2;
    rr_rewindow2d(win);
    win->redisp(0);
}

void rr_pan2d(rr_window2d* win, double dx, double dy) {
    win->xcen += dx;
    win->ycen += dy;
    win->xmin = win->xcen-win->xran2;
    win->xmax = win->xcen+win->xran2;
    win->ymin = win->ycen-win->yran2;
    win->ymax = win->ycen+win->yran2;
    rr_rewindow2d( win );
    win->redisp(0);
}

double rr_xmousepoint2d(rr_window2d* win, int x) {
    return win->xmin+(win->xmax-win->xmin)*((double)x/win->width);
}
double rr_ymousepoint2d(rr_window2d* win, int y) {
    return win->ymin+(win->ymax-win->ymin)*(1.0-(double)y/win->height);
}



