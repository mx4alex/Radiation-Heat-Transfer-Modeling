module surface_properties;

import spectral_constants;

#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>
#include <mutex>

std::vector<Object> g_objects;
std::vector<GlobalFace> g_allFaces;
std::vector<std::vector<double>> g_viewFactors;
std::mutex simMutex;
bool simulationRunning = false;

static double dot(const Vector3& a, const Vector3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static Vector3 cross(const Vector3& a, const Vector3& b) {
    return Vector3(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}

static double length(const Vector3& v) {
    return std::sqrt(dot(v, v));
}

static Vector3 normalize(const Vector3& v) {
    double len = length(v);
    return Vector3(v.x / len, v.y / len, v.z / len);
}

void updateGlobalData() {
    g_allFaces.clear();
    for (size_t i = 0; i < g_objects.size(); i++) {
        for (size_t j = 0; j < g_objects[i].faces.size(); j++) {
            g_allFaces.push_back({static_cast<int>(i), static_cast<int>(j)});
        }
    }
}

void setObjectProperties(int objectIndex, double temperature, double emissivity) {
    if (objectIndex < 0 || objectIndex >= static_cast<int>(g_objects.size())) {
        throw std::out_of_range("Invalid object index");
    }

    Object& obj = g_objects[objectIndex];
    
    for (auto& face : obj.faces) {
        face.temp = temperature;
    }
    
    for (int i = 0; i < NUM_BANDS; i++) {
        obj.emissivity[i] = emissivity;
        obj.reflectivity[i] = 1.0 - emissivity;
    }
    
    obj.initialTemperature = temperature;
}

