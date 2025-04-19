#ifndef VISUALIZER_H
#define VISUALIZER_H

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "surface_properties.h"

extern int windowWidth;
extern int windowHeight;

void initGL();
void updateMinMaxTemps();
void renderAllFaces(double Tmin, double Tmax);
void renderColorBar(int w, int h, double Tmin, double Tmax);

struct Quaternion {
    double w, x, y, z;
};

extern bool   g_arcballDragging;
extern double g_virtualX, g_virtualY;
extern float  g_zoom;
extern bool   g_rightDown;
extern int    g_lastMouseX, g_lastMouseY;

extern double currentMinTemp;
extern double currentMaxTemp;
extern Quaternion g_sceneQuat;

void mouseMotionArcball(int x, int y,
                        bool &arcballDragging,
                        bool &rightDown,
                        int &lastMouseX,
                        int &lastMouseY,
                        double &virtualX,
                        double &virtualY,
                        float &zoom,
                        Quaternion &sceneQuat);

#endif
