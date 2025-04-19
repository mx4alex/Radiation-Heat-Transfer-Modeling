#ifndef SURFACE_PROPERTIES_H
#define SURFACE_PROPERTIES_H

#include <vector>
#include <string>
#include <mutex>
#include <cmath>

struct Vector3 {
    double x, y, z;
};

inline Vector3 operator-(const Vector3& a, const Vector3& b){
    return { a.x - b.x, a.y - b.y, a.z - b.z };
}
inline Vector3 operator+(const Vector3& a, const Vector3& b){
    return { a.x + b.x, a.y + b.y, a.z + b.z };
}
inline Vector3 operator*(const Vector3& a, double s){
    return { a.x * s, a.y * s, a.z * s };
}
inline double dot(const Vector3& a, const Vector3& b){
    return a.x*b.x + a.y*b.y + a.z*b.z;
}
inline Vector3 cross(const Vector3& a, const Vector3& b){
    return { a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x };
}
inline double length(const Vector3& a){
    return std::sqrt(dot(a,a));
}
inline Vector3 normalize(const Vector3& a){
    double l = length(a);
    if(l < 1e-14) return {0,0,0};
    return { a.x/l, a.y/l, a.z/l };
}

struct Face {
    int v1, v2, v3;
    Vector3 normal;
    double area;
    double temp;
    std::vector<int> neighbors;
};

enum class BoundaryConditionType {
    First,
    Second
};

struct Sphere {
    std::vector<Vector3> vertices;
    std::vector<Face>    faces;
    double emissivity[3];
    double reflectivity[3];
    BoundaryConditionType boundaryCondition;
    double boundaryValue;
    std::string stlFileName;
    Vector3 position;
    double initialTemperature;
    double mass;
    double cp;
    double totalArea;
};

struct GlobalFace {
    int sphereIndex;
    int faceIndex;
};

extern std::vector<Sphere>              g_spheres;
extern std::vector<GlobalFace>          g_allFaces;
extern std::vector<std::vector<double>> g_viewFactors;
extern std::mutex                       simMutex;
extern bool                             simulationRunning;

#endif 
