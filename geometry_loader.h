#ifndef GEOMETRY_LOADER_H
#define GEOMETRY_LOADER_H

#include <string>
#include <vector>

#include "surface_properties.h"

bool loadASCIISTL(const std::string &filename,
                  std::vector<Vector3> &vertices,
                  std::vector<Face> &faces);

void buildGlobalFaceList();
void buildFaceAdjacency();
void updateGlobalData();

#endif 
