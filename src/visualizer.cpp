module visualizer;

import surface_properties;
import geometry_loader;

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <cmath>
#include <algorithm>
#include <sstream>
#include <iostream>

int windowWidth  = 1280;
int windowHeight = 720;

bool   g_arcballDragging = false;
double g_virtualX=0.0, g_virtualY=0.0;
float  g_zoom=1.f;
bool   g_rightDown=false;
int    g_lastMouseX=0, g_lastMouseY=0;

double currentMinTemp = 0.0;
double currentMaxTemp = 5000.0;

Quaternion g_sceneQuat = {1.0, 0.0, 0.0, 0.0};

static void quatToMatrix(const Quaternion &q, double m[16]){
    double xx=q.x*q.x, yy=q.y*q.y, zz=q.z*q.z;
    double xy=q.x*q.y, xz=q.x*q.z, yz=q.y*q.z;
    double wx=q.w*q.x, wy=q.w*q.y, wz=q.w*q.z;
    m[0] = 1-2*(yy+zz); m[1] = 2*(xy+wz);   m[2] = 2*(xz-wy);   m[3] = 0;
    m[4] = 2*(xy-wz);   m[5] = 1-2*(xx+zz); m[6] = 2*(yz+wx);   m[7] = 0;
    m[8] = 2*(xz+wy);   m[9] = 2*(yz-wx);   m[10]= 1-2*(xx+yy); m[11]= 0;
    m[12]=0; m[13]=0; m[14]=0; m[15]=1;
}

static float clamp01(float v){
    if(v<0.f) v=0.f; if(v>1.f) v=1.f;
    return v;
}
struct ColorRGB{
    float r,g,b;
};
static ColorRGB hsvToRgb(float h, float s, float v){
    float c = v*s;
    float x = c*(1.f - std::fabs(std::fmod(h/60.f,2.f)-1.f));
    float m = v-c;
    float r_,g_,b_;
    if(h>=0 && h<60){r_=c; g_=x; b_=0;}
    else if(h>=60 && h<120){r_=x; g_=c; b_=0;}
    else if(h>=120 && h<180){r_=0; g_=c; b_=x;}
    else if(h>=180 && h<240){r_=0; g_=x; b_=c;}
    else if(h>=240 && h<300){r_=x; g_=0; b_=c;}
    else {r_=c; g_=0; b_=x;}
    return {r_+m, g_+m, b_+m};
}
static ColorRGB temperatureToColor(double T, double Tmin, double Tmax){
    if(Tmax<=Tmin) return {1,1,1};
    double alpha = (T - Tmin)/(Tmax - Tmin);
    if(alpha<0.) alpha=0.; if(alpha>1.) alpha=1.;
    float hue = 240.f*(1.f-(float)alpha);
    return hsvToRgb(hue, 1.f, 1.f);
}

Color temperatureToColor(double temp) {
    if (std::isnan(temp) || std::isinf(temp) || temp < 0) {
        std::cerr << "Warning: Invalid temperature value detected: " << temp << std::endl;
        return Color(1.0f, 0.0f, 0.0f);
    }
    
    if (temp < 0.0) temp = 0.0;
    if (temp > 5000.0) temp = 5000.0;

    double t = (temp - currentMinTemp) / (currentMaxTemp - currentMinTemp);
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;

    return Color(
        static_cast<float>(t),
        static_cast<float>(1.0 - t),
        0.0f
    );
}

void initGL(){
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.9f,0.9f,0.9f,1.f);
}

void updateMinMaxTemps(){
    if(g_allFaces.empty()) return;
    double mi=1e30, ma=-1e30;
    for(auto &gf: g_allFaces){
        double T = g_objects[gf.objectIndex].faces[gf.faceIndex].temp;
        if(T<mi) mi=T;
        if(T>ma) ma=T;
    }
    double pad=50.;
    currentMinTemp= std::max(0., mi-pad);
    currentMaxTemp= std::min(5000.,ma+pad);
}

void renderAllFaces(double Tmin, double Tmax){
    glBegin(GL_TRIANGLES);
    for(auto &gf: g_allFaces){
        int si = gf.objectIndex;
        int fi = gf.faceIndex;
        auto &sp = g_objects[si];
        auto &fc = sp.faces[fi];

        Vector3 v1= sp.vertices[fc.v1];
        Vector3 v2= sp.vertices[fc.v2];
        Vector3 v3= sp.vertices[fc.v3];

        ColorRGB col = temperatureToColor(fc.temp, Tmin, Tmax);
        glColor3f(col.r, col.g, col.b);

        glVertex3f((float)v1.x, (float)v1.y, (float)v1.z);
        glVertex3f((float)v2.x, (float)v2.y, (float)v2.z);
        glVertex3f((float)v3.x, (float)v3.y, (float)v3.z);
    }
    glEnd();
}

void renderColorBar(int w, int h, double Tmin, double Tmax){
    int barW = 50;
    int barH = h-100;
    int barX = 50;
    int barY = 50;

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0,w,0,h);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    glBegin(GL_QUADS);
    for(int i=0; i<barH; i++){
        double t = double(i)/double(barH);
        double T = Tmin + t*(Tmax - Tmin);
        ColorRGB c = temperatureToColor(T,Tmin,Tmax);
        glColor3f(c.r,c.g,c.b);
        glVertex2f(barX, barY+i);
        glVertex2f(barX+barW, barY+i);
        glVertex2f(barX+barW, barY+i+1);
        glVertex2f(barX, barY+i+1);
    }
    glEnd();

    glColor3f(0,0,0);
    glBegin(GL_LINE_LOOP);
    glVertex2f(barX,   barY);
    glVertex2f(barX+barW, barY);
    glVertex2f(barX+barW, barY+barH);
    glVertex2f(barX,   barY+barH);
    glEnd();

    void* font = GLUT_BITMAP_HELVETICA_12;
    int nLabels=5;
    for(int i=0; i<nLabels; i++){
        double f= double(i)/(nLabels-1);
        double TT = Tmin + f*(Tmax - Tmin);
        std::ostringstream oss; oss<<(int)TT<<"K";
        std::string txt = oss.str();
        int yPos = barY + (int)(f*barH);
        glRasterPos2i(barX+barW+10, yPos);
        for(char c: txt){
            glutBitmapCharacter(font, c);
        }
    }

    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
}

static Quaternion quatMultiply(const Quaternion &a, const Quaternion &b){
    return {
        a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z,
        a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y,
        a.w*b.y - a.x*b.z + a.y*b.w + a.z*b.x,
        a.w*b.z + a.x*b.y - a.y*b.x + a.z*b.w
    };
}
static Quaternion quatNormalize(const Quaternion &q){
    double len= std::sqrt(q.w*q.w + q.x*q.x + q.y*q.y + q.z*q.z);
    return {q.w/len, q.x/len, q.y/len, q.z/len};
}
static Quaternion quatFromAxisAngle(const Vector3 &axis, double angle){
    Vector3 na = normalize(axis);
    double s = std::sin(angle/2.0);
    return {std::cos(angle/2.0), na.x*s, na.y*s, na.z*s};
}
static Vector3 getArcballVectorFromVirtual(double vx, double vy){
    double radius  = std::min(windowWidth, windowHeight)/2.0;
    double centerX = windowWidth/2.0;
    double centerY = windowHeight/2.0;
    double dx = (vx - centerX)/radius;
    double dy = (centerY - vy)/radius;
    double r  = std::sqrt(dx*dx + dy*dy);
    double z  = (r<=1.0)? std::sqrt(1.0-r*r) : 0.0;
    Vector3 v={dx,dy,z};
    return normalize(v);
}

void mouseMotionArcball(int x, int y,
                        bool &arcballDragging,
                        bool &rightDown,
                        int &lastMouseX,
                        int &lastMouseY,
                        double &virtualX,
                        double &virtualY,
                        float &zoom,
                        Quaternion &sceneQuat)
{
    int simY = windowHeight - y;
    static int prevSimY = windowHeight - lastMouseY;
    int dy = simY - prevSimY;
    prevSimY = simY;

    if(rightDown){
        zoom*= (1.f + dy*0.01f);
        if(zoom<0.1f) zoom=0.1f;
    }
    if(arcballDragging){
        int deltaX = x - lastMouseX;
        int deltaY = y - lastMouseY;
        virtualX+= deltaX;
        virtualY+= deltaY;

        Vector3 lastVec = getArcballVectorFromVirtual(virtualX - deltaX, virtualY - deltaY);
        Vector3 curVec  = getArcballVectorFromVirtual(virtualX, virtualY);
        Vector3 axis= cross(lastVec, curVec);
        double dotVal= dot(lastVec, curVec);
        if(dotVal<-1.) dotVal=-1.; if(dotVal>1.) dotVal=1.;
        double angle = std::acos(dotVal);
        double lAxis = length(axis);
        if(lAxis>1e-6 && angle!=0.0){
            Quaternion dq= quatFromAxisAngle(axis,angle);
            sceneQuat= quatMultiply(dq, sceneQuat);
            sceneQuat= quatNormalize(sceneQuat);
        }
    }
    lastMouseX= x;
    lastMouseY= y;
}
