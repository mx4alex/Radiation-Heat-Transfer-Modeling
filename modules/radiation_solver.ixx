module;

import surface_properties;
import spectral_constants;

#include <vector>
#include <unordered_map>
#include <utility>
#include <cstdint>
#include <array>

export module radiation_solver;

export enum class RadiationModel {
    MONTE_CARLO,
    DISCRETE_ORDINATES,
    S2S
};

export struct Vector3;
export struct Object;
export struct Face;
export struct GlobalFace;
export struct BVHNode;
export struct TriInfo;

export extern std::vector<Object> g_objects;
export extern std::vector<GlobalFace> g_allFaces;
export extern std::unordered_map<uint64_t, int> faceLookup;
export extern std::vector<std::vector<double>> g_viewFactors;
export extern std::vector<TriInfo> g_triInfos;
export extern RadiationModel g_model;
export extern BVHNode* g_bvhRoot;

export extern int monteCarloSamples;
export extern int discreteOrdinates_nTheta;
export extern int discreteOrdinates_nPhi;
export extern int s2sSubDiv;

export inline uint64_t makeKey(int objectIndex, int faceIndex) {
    return (static_cast<uint64_t>(objectIndex) << 32) | (static_cast<uint32_t>(faceIndex));
}

export void buildBlackbodyFractionTable();
export double getBlackbodyFraction(int band, double temperature);
export void computeDONodes(int nTheta, int nPhi);
export void computeTriInfos(int subDiv);
export BVHNode* buildBVH(std::vector<int>& indices);
export void deleteBVH(BVHNode* node);

export void computeViewFactors_S2S(int subDiv);
export void radiationFull(double alpha);
export void buildFaceGeomCache(); 