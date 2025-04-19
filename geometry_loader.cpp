#include "geometry_loader.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <cctype>

bool loadASCIISTL(const std::string &filename,
                  std::vector<Vector3> &vertices,
                  std::vector<Face> &faces)
{
    std::ifstream ifs(filename);
    if(!ifs){
        std::cerr << "ERROR: cannot open file: " << filename << "\n";
        return false;
    }
    vertices.clear();
    faces.clear();

    std::string line;
    bool inFacet = false;
    std::vector<Vector3> facetVerts;

    while(std::getline(ifs, line)){
        std::string lowerLine = line;
        std::transform(lowerLine.begin(), lowerLine.end(), lowerLine.begin(),
                       [](unsigned char c){ return std::tolower(c); });

        auto trim = [](std::string &s){
            while(!s.empty() && std::isspace((unsigned char)s.front())) s.erase(s.begin());
            while(!s.empty() && std::isspace((unsigned char)s.back()))  s.pop_back();
        };
        trim(lowerLine);

        if(lowerLine.rfind("facet normal", 0) == 0){
            inFacet = true;
            facetVerts.clear();
        }
        else if(lowerLine.rfind("vertex", 0) == 0 && inFacet){
            Vector3 v;
            std::istringstream iss(line);
            std::string vertexWord;
            iss >> vertexWord >> v.x >> v.y >> v.z;
            facetVerts.push_back(v);
        }
        else if(lowerLine.rfind("endfacet", 0) == 0 && inFacet){
            inFacet = false;
            if(facetVerts.size() == 3){
                int baseIndex = (int)vertices.size();
                vertices.push_back(facetVerts[0]);
                vertices.push_back(facetVerts[1]);
                vertices.push_back(facetVerts[2]);

                Face fc;
                fc.v1 = baseIndex;
                fc.v2 = baseIndex+1;
                fc.v3 = baseIndex+2;

                Vector3 e1 = facetVerts[1] - facetVerts[0];
                Vector3 e2 = facetVerts[2] - facetVerts[0];
                Vector3 n = cross(e1, e2);
                double ln = length(n);

                if(ln > 1e-14){
                    Vector3 centroid = (facetVerts[0] + facetVerts[1] + facetVerts[2])*(1./3.);
                    double dotVal = dot(n, centroid);
                    if(dotVal < 0.0){
                        n = { -n.x, -n.y, -n.z };
                    }
                    fc.normal = normalize(n);
                    fc.area   = 0.5 * ln;
                } else {
                    fc.normal = {0,0,1};
                    fc.area   = 0.0;
                }
                fc.temp = 300.0;

                faces.push_back(fc);
            }
        }
    }

    if(faces.empty()){
        std::cerr << "No faces found in ASCII STL: " << filename << "\n";
        return false;
    }
    return true;
}

void buildGlobalFaceList(){
    g_allFaces.clear();
    for(int si = 0; si < (int)g_spheres.size(); si++){
        for(int fi = 0; fi < (int)g_spheres[si].faces.size(); fi++){
            GlobalFace gf;
            gf.sphereIndex = si;
            gf.faceIndex   = fi;
            g_allFaces.push_back(gf);
        }
    }
    std::cout << "Total faces (all spheres) = " << g_allFaces.size() << "\n";
}

void buildFaceAdjacency(){
    typedef std::pair<int,int> Edge;
    auto edgeHash = [&](const Edge &a){
        int mn = (a.first < a.second) ? a.first : a.second;
        int mx = (a.first < a.second) ? a.second : a.first;
        return std::to_string(mn) + "_" + std::to_string(mx);
    };

    std::unordered_map<std::string, std::vector<int>> edgeToFaces;

    int globalFaceIndex = 0;
    for(int si = 0; si < (int)g_spheres.size(); si++){
        for(int fi = 0; fi < (int)g_spheres[si].faces.size(); fi++){
            auto &fc = g_spheres[si].faces[fi];
            Edge e1 = { fc.v1, fc.v2 };
            Edge e2 = { fc.v2, fc.v3 };
            Edge e3 = { fc.v3, fc.v1 };
            edgeToFaces[edgeHash(e1)].push_back(globalFaceIndex);
            edgeToFaces[edgeHash(e2)].push_back(globalFaceIndex);
            edgeToFaces[edgeHash(e3)].push_back(globalFaceIndex);
            globalFaceIndex++;
        }
    }

    globalFaceIndex = 0;
    for(int si = 0; si < (int)g_spheres.size(); si++){
        for(int fi = 0; fi < (int)g_spheres[si].faces.size(); fi++){
            auto &fc = g_spheres[si].faces[fi];

            auto addNeighbors = [&](const Edge &E){
                auto &vec = edgeToFaces[edgeHash(E)];
                for(int fidx : vec){
                    if(fidx != globalFaceIndex)
                        fc.neighbors.push_back(fidx);
                }
            };

            Edge e1 = { fc.v1, fc.v2 };
            Edge e2 = { fc.v2, fc.v3 };
            Edge e3 = { fc.v3, fc.v1 };

            addNeighbors(e1);
            addNeighbors(e2);
            addNeighbors(e3);

            globalFaceIndex++;
        }
    }
    std::cout << "Face adjacency built.\n";
}

void updateGlobalData(){
    for(auto &s : g_spheres){
        for(auto &f : s.faces){
            f.neighbors.clear();
        }
    }
    buildGlobalFaceList();
    buildFaceAdjacency();
}
