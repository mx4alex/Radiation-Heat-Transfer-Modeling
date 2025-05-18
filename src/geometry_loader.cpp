module geometry_loader;

import surface_properties;
import radiation_solver;
import spectral_constants;

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <cctype>

std::unordered_map<uint64_t, int> faceLookup;

bool loadASCIISTL(std::istream &ifs,
                  std::vector<Vector3> &vertices,
                  std::vector<Face>    &faces)
{
    vertices.clear(); faces.clear();
    std::string line; bool inFacet = false;
    std::vector<Vector3> facetVerts;

    auto trim = [](std::string &s){
        while(!s.empty() && std::isspace(static_cast<unsigned char>(s.front()))) s.erase(s.begin());
        while(!s.empty() && std::isspace(static_cast<unsigned char>(s.back())))  s.pop_back();
    };

    while (std::getline(ifs, line)) {
        std::string lower = line;
        std::transform(lower.begin(), lower.end(), lower.begin(),
                       [](unsigned char c){ return std::tolower(c); });
        trim(lower);

        if (lower.rfind("facet normal", 0) == 0) {
            inFacet = true; facetVerts.clear();
        }
        else if (lower.rfind("vertex", 0) == 0 && inFacet) {
            Vector3 v; std::istringstream iss(line); std::string tmp;
            iss >> tmp >> v.x >> v.y >> v.z;
            facetVerts.push_back(v);
        }
        else if (lower.rfind("endfacet", 0) == 0 && inFacet) {
            inFacet = false;
            if (facetVerts.size() == 3) {
                int base = vertices.size();
                vertices.push_back(facetVerts[0]);
                vertices.push_back(facetVerts[1]);
                vertices.push_back(facetVerts[2]);

                Face f; f.v1 = base; f.v2 = base+1; f.v3 = base+2;

                Vector3 n = cross(facetVerts[1]-facetVerts[0],
                                  facetVerts[2]-facetVerts[0]);
                double ln = length(n);
                if (ln > 1e-14) {
                    Vector3 centroid = (facetVerts[0]+facetVerts[1]+facetVerts[2]) * (1.0/3.0);
                    if (dot(n, centroid) < 0.0) n = {-n.x,-n.y,-n.z};
                    f.normal = normalize(n);
                    f.area   = 0.5 * ln;
                } else {
                    f.normal = {0,0,1}; f.area = 0.0;
                }
                f.temp = 300.0;
                faces.push_back(f);
            }
        }
    }
    return !faces.empty();
}

static constexpr std::uint32_t MAX_TRIANGLES = 12'000'000;

static bool isBinarySTL(std::ifstream &ifs, std::uint64_t fileSize,
                        std::uint32_t &triCountOut)
{
    if (fileSize < 84) return false;

    ifs.seekg(80, std::ios::beg);
    std::uint32_t nHeader = 0;
    ifs.read(reinterpret_cast<char*>(&nHeader), 4);
    if (!ifs) return false;

    const std::uint64_t bytesAfter = fileSize - 84;
    const bool ok       = (bytesAfter % 50ull) == 0ull &&
                          (bytesAfter / 50ull) == nHeader;

    if (ok && nHeader != 0 && nHeader <= MAX_TRIANGLES) {
        triCountOut = nHeader;
        return true;
    }
    return false;
}

static bool loadBinarySTL(std::istream &ifs, std::uint32_t triCount,
                          std::vector<Vector3> &vertices,
                          std::vector<Face>    &faces)
{
    try {
        vertices.clear();  faces.clear();
        vertices.reserve(static_cast<std::size_t>(triCount) * 3u);
        faces.reserve   (static_cast<std::size_t>(triCount));
    } catch (const std::bad_alloc&) {
        std::cerr << "File contains " << triCount
                  << " triangles — insufficient memory.\n";
        return false;
    }

    for (std::uint32_t i = 0; i < triCount; ++i)
    {
        float n[3], v[9];  std::uint16_t attr;
        ifs.read(reinterpret_cast<char*>(n), 12);
        ifs.read(reinterpret_cast<char*>(v), 36);
        ifs.read(reinterpret_cast<char*>(&attr), 2);
        if (!ifs) {
            std::cerr << "Unexpected end of file at triangle "
                      << i << ".\n";
            return false;
        }

        const int base = static_cast<int>(vertices.size());
        vertices.emplace_back(Vector3{v[0], v[1], v[2]});
        vertices.emplace_back(Vector3{v[3], v[4], v[5]});
        vertices.emplace_back(Vector3{v[6], v[7], v[8]});

        Face f{ base, base + 1, base + 2 };
        Vector3 nvec{ n[0], n[1], n[2] };

        double ln = length(nvec);
        if (ln < 1e-14) {
            nvec = cross(vertices[base + 1] - vertices[base],
                         vertices[base + 2] - vertices[base]);
            ln = length(nvec);
        }
        f.normal = (ln > 1e-14) ? normalize(nvec) : Vector3{0,0,1};
        f.area   = (ln > 1e-14) ? 0.5 * ln : 0.0;
        f.temp   = 300.0;

        faces.push_back(f);
    }
    return true;
}

bool loadSTL(const std::string &filename,
             std::vector<Vector3> &vertices,
             std::vector<Face>    &faces)
{
    std::ifstream ifs(filename, std::ios::binary);
    if (!ifs) {
        std::cerr << "Cannot open file " << filename << '\n';
        return false;
    }

    ifs.seekg(0, std::ios::end);
    const std::uint64_t fileSize = static_cast<std::uint64_t>(ifs.tellg());
    ifs.seekg(0, std::ios::beg);

    std::uint32_t triCount = 0;
    const bool binary = isBinarySTL(ifs, fileSize, triCount);

    if (binary) {
        ifs.seekg(84, std::ios::beg);
        return loadBinarySTL(ifs, triCount, vertices, faces);
    }
    return loadASCIISTL(ifs, vertices, faces);
}

void buildGlobalFaceList()
{
    g_allFaces.clear();
    faceLookup.clear();
    int globalFaceIndex = 0;
    for (int si = 0; si < static_cast<int>(g_objects.size()); si++) {
        for (int fi = 0; fi < static_cast<int>(g_objects[si].faces.size()); fi++) {
            auto &fc = g_objects[si].faces[fi];
            GlobalFace gf{si, fi};
            g_allFaces.push_back(gf);
            faceLookup[makeKey(si, fi)] = globalFaceIndex++;
        }
    }
    std::cout << "Total faces (all objects) = " << g_allFaces.size() << "\n";
}

void buildFaceAdjacency()
{
    typedef std::pair<int,int> Edge;
    auto edgeHash = [&](const Edge &a) {
        int mn = std::min(a.first, a.second);
        int mx = std::max(a.first, a.second);
        return std::to_string(mn) + "_" + std::to_string(mx);
    };
    std::unordered_map<std::string, std::vector<int>> edgeToFaces;
    int globalFaceIndex = 0;
    for (int si = 0; si < static_cast<int>(g_objects.size()); si++) {
        for (int fi = 0; fi < static_cast<int>(g_objects[si].faces.size()); fi++) {
            auto &fc = g_objects[si].faces[fi];
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
    for (int si = 0; si < static_cast<int>(g_objects.size()); si++) {
        for (int fi = 0; fi < static_cast<int>(g_objects[si].faces.size()); fi++) {
            auto &fc = g_objects[si].faces[fi];
            auto addNeighbors = [&](const Edge &E) {
                for (int fidx : edgeToFaces[edgeHash(E)]) {
                    if (fidx != globalFaceIndex)
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
    for(auto &s : g_objects)
        for(auto &f : s.faces)
            f.neighbors.clear();
    buildGlobalFaceList();
    buildFaceAdjacency();
    std::vector<int> faceIndices(g_allFaces.size());
    for (int i = 0; i < static_cast<int>(g_allFaces.size()); i++){
        faceIndices[i] = i;
    }
    if(g_bvhRoot) {
        deleteBVH(g_bvhRoot);
        g_bvhRoot = nullptr;
    }
    g_bvhRoot = buildBVH(faceIndices);
}
