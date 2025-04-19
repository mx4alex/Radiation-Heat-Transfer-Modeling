#include "surface_properties.h"

std::vector<Sphere>              g_spheres;
std::vector<GlobalFace>          g_allFaces;
std::vector<std::vector<double>> g_viewFactors;
std::mutex                       simMutex;
bool                             simulationRunning = false;
