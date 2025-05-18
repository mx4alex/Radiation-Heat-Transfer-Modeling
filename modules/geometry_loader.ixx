module;

#include <string>
#include <vector>
#include "surface_properties.h"

export module geometry_loader;

export {
    bool loadSTL(const std::string &filename,
                 std::vector<Vector3> &vertices,
                 std::vector<Face>    &faces);
        
    void buildGlobalFaceList();
    void buildFaceAdjacency();
    void updateGlobalData();
} 