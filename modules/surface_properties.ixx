module;

#include <vector>
#include <string>
#include <mutex>
#include <cmath>

export module surface_properties;

export struct Vector3 {
    double x, y, z;
    Vector3(double x_ = 0.0, double y_ = 0.0, double z_ = 0.0) : x(x_), y(y_), z(z_) {}
};

export inline Vector3 operator-(const Vector3& a, const Vector3& b){
    return { a.x - b.x, a.y - b.y, a.z - b.z };
}
export inline Vector3 operator+(const Vector3& a, const Vector3& b){
    return { a.x + b.x, a.y + b.y, a.z + b.z };
}
export inline Vector3 operator*(const Vector3& a, double s){
    return { a.x * s, a.y * s, a.z * s };
}
export inline double dot(const Vector3& a, const Vector3& b){
    return a.x*b.x + a.y*b.y + a.z*b.z;
}
export inline Vector3 cross(const Vector3& a, const Vector3& b){
    return { a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x };
}
export inline double length(const Vector3& a){
    return std::sqrt(dot(a,a));
}
export inline Vector3 normalize(const Vector3& a){
    double l = length(a);
    if(l < 1e-14) return {0,0,0};
    return { a.x/l, a.y/l, a.z/l };
}

export struct Face {
    int v1, v2, v3;
    Vector3 normal;
    double area;
    double temp;
    std::vector<int> neighbors;
};

export enum class BoundaryConditionType {
    First,
    Second
};

export struct Object {
    std::vector<Vector3> vertices;
    std::vector<Face>    faces;
    std::vector<double>  emissivity;
    std::vector<double>  reflectivity;
    BoundaryConditionType boundaryCondition;
    double boundaryValue;
    std::string stlFileName;
    Vector3 position;
    double initialTemperature;
    double mass;
    double cp;
    double totalArea;
};

export struct GlobalFace {
    int objectIndex;
    int faceIndex;
};

export extern std::vector<Object>              g_objects;
export extern std::vector<GlobalFace>          g_allFaces;
export extern std::vector<std::vector<double>> g_viewFactors;
export extern std::mutex                       simMutex;
export extern bool                             simulationRunning;
export extern int                              g_numWavelengths;

export void updateGlobalData();
export void setObjectProperties(int objectIndex, double temperature, double emissivity);