module;

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

import surface_properties;
import geometry_loader;

export module visualizer;

export extern int windowWidth;
export extern int windowHeight;

export void initGL();
export void updateMinMaxTemps();
export void renderAllFaces(double Tmin, double Tmax);
export void renderColorBar(int w, int h, double Tmin, double Tmax);

export struct Color {
    float r, g, b;
    Color(float r_ = 0.0f, float g_ = 0.0f, float b_ = 0.0f) : r(r_), g(g_), b(b_) {}
};

export struct Quaternion {
    double w, x, y, z;
};

export Color temperatureToColor(double temp);

export extern bool   g_arcballDragging;
export extern double g_virtualX, g_virtualY;
export extern float  g_zoom;
export extern bool   g_rightDown;
export extern int    g_lastMouseX, g_lastMouseY;

export extern double currentMinTemp;
export extern double currentMaxTemp;
export extern Quaternion g_sceneQuat;

export void mouseMotionArcball(int x, int y,
                        bool &arcballDragging,
                        bool &rightDown,
                        int &lastMouseX,
                        int &lastMouseY,
                        double &virtualX,
                        double &virtualY,
                        float &zoom,
                        Quaternion &sceneQuat); 