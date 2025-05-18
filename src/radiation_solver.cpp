#define _USE_MATH_DEFINES

module radiation_solver;

import surface_properties;
import spectral_constants;

#include <random>
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <thread>
#include <atomic>
#include <mutex>
#include <algorithm>
#include <fstream>
#include <iomanip>

int monteCarloSamples       = 1000;
int discreteOrdinates_nTheta= 8;
int discreteOrdinates_nPhi  = 16;
int s2sSubDiv               = 4;

RadiationModel g_model = RadiationModel::S2S;

struct Ray {
    Vector3 origin;
    Vector3 dir;
};

std::pair<int,int> traverseBVH(BVHNode* node, const Ray &ray, double &closest_t);

static std::pair<int,int> traceAllObjects(const Ray &ray) {
    double closest_t = 1e30;
    return traverseBVH(g_bvhRoot, ray, closest_t);
}

struct TriInfo {
    std::vector<Vector3> samples;
    Vector3 normal;
    double area;
};
std::vector<TriInfo> g_triInfos;

void computeTriInfos(int subDiv) {
    int NF = static_cast<int>(g_allFaces.size());
    g_triInfos.resize(NF);
    for (int i = 0; i < NF; i++){
        auto &gf = g_allFaces[i];
        auto &sp = g_objects[gf.objectIndex];
        auto &fc = sp.faces[gf.faceIndex];
        TriInfo info;
        info.normal = fc.normal;
        info.area = fc.area;
        std::vector<Vector3> pts;
        pts.reserve(subDiv * subDiv);
        Vector3 A = sp.vertices[fc.v1];
        Vector3 B = sp.vertices[fc.v2];
        Vector3 C = sp.vertices[fc.v3];
        for (int iu = 0; iu < subDiv; iu++){
            for (int iv = 0; iv < subDiv; iv++){
                double u = (iu + 0.5) / subDiv;
                double v = (iv + 0.5) / subDiv;
                if(u + v > 1.0) continue;
                double w = 1.0 - u - v;
                Vector3 p = { A.x * w + B.x * u + C.x * v,
                              A.y * w + B.y * u + C.y * v,
                              A.z * w + B.z * u + C.z * v };
                pts.push_back(p);
            }
        }
        info.samples = std::move(pts);
        g_triInfos[i] = std::move(info);
    }
}

struct DONode {
    double ldx, ldy, ldz;
    double weight;
};
std::vector<DONode> g_DO_nodes;
void computeDONodes(int nTheta, int nPhi) {
    g_DO_nodes.clear();
    double dTheta = (M_PI / 2.0) / nTheta;
    double dPhi = (2.0 * M_PI) / nPhi;
    for (int it = 0; it < nTheta; it++){
        double theta = (it + 0.5) * dTheta;
        double sinTheta = std::sin(theta);
        double cosTheta = std::cos(theta);
        for (int ip = 0; ip < nPhi; ip++){
            double phi = (ip + 0.5) * dPhi;
            double ldx = sinTheta * std::cos(phi);
            double ldy = sinTheta * std::sin(phi);
            double ldz = cosTheta;
            double weight = sinTheta * dTheta * dPhi;
            g_DO_nodes.push_back({ ldx, ldy, ldz, weight });
        }
    }
}

struct AABB {
    Vector3 min;
    Vector3 max;
    AABB() {
        min = { 1e30, 1e30, 1e30 };
        max = { -1e30, -1e30, -1e30 };
    }
    void expand(const Vector3 &p) {
        min.x = std::min(min.x, p.x);
        min.y = std::min(min.y, p.y);
        min.z = std::min(min.z, p.z);
        max.x = std::max(max.x, p.x);
        max.y = std::max(max.y, p.y);
        max.z = std::max(max.z, p.z);
    }
    bool intersect(const Ray &ray, double t_min, double t_max) const {
        for (int i = 0; i < 3; i++) {
            double origin, dir, minVal, maxVal;
            if(i == 0){ origin = ray.origin.x; dir = ray.dir.x; minVal = min.x; maxVal = max.x; }
            else if(i == 1){ origin = ray.origin.y; dir = ray.dir.y; minVal = min.y; maxVal = max.y; }
            else { origin = ray.origin.z; dir = ray.dir.z; minVal = min.z; maxVal = max.z; }
            double invD = 1.0 / dir;
            double t0 = (minVal - origin) * invD;
            double t1 = (maxVal - origin) * invD;
            if(invD < 0) std::swap(t0, t1);
            t_min = std::max(t0, t_min);
            t_max = std::min(t1, t_max);
            if(t_max <= t_min)
                return false;
        }
        return true;
    }
};

struct BVHNode {
    AABB box;
    BVHNode *left = nullptr;
    BVHNode *right = nullptr;
    std::vector<int> faceIndices;
};

AABB computeAABBForFace(int faceIndex) {
    AABB box;
    auto &gf = g_allFaces[faceIndex];
    auto &sp = g_objects[gf.objectIndex];
    auto &fc = sp.faces[gf.faceIndex];
    Vector3 v1 = sp.vertices[fc.v1];
    Vector3 v2 = sp.vertices[fc.v2];
    Vector3 v3 = sp.vertices[fc.v3];
    box.expand(v1);
    box.expand(v2);
    box.expand(v3);
    return box;
}

BVHNode* buildBVH(std::vector<int>& indices, int start, int end) {
    BVHNode* node = new BVHNode();
    int count = end - start;
    AABB box;
    for (int i = start; i < end; i++){
        AABB b = computeAABBForFace(indices[i]);
        box.expand(b.min);
        box.expand(b.max);
    }
    node->box = box;
    const int threshold = 4;
    if(count <= threshold) {
        for (int i = start; i < end; i++){
            node->faceIndices.push_back(indices[i]);
        }
    } else {
        AABB centroidBox;
        for (int i = start; i < end; i++){
            int idx = indices[i];
            auto &gf = g_allFaces[idx];
            auto &sp = g_objects[gf.objectIndex];
            auto &fc = sp.faces[gf.faceIndex];
            Vector3 v1 = sp.vertices[fc.v1];
            Vector3 v2 = sp.vertices[fc.v2];
            Vector3 v3 = sp.vertices[fc.v3];
            Vector3 centroid = (v1 + v2 + v3) * (1.0 / 3.0);
            centroidBox.expand(centroid);
        }
        Vector3 extent = { centroidBox.max.x - centroidBox.min.x,
                           centroidBox.max.y - centroidBox.min.y,
                           centroidBox.max.z - centroidBox.min.z };
        int axis = (extent.y > extent.x && extent.y > extent.z) ? 1 : (extent.z > extent.x ? 2 : 0);
        auto comparator = [axis](int a, int b) {
            auto &gfA = g_allFaces[a];
            auto &spA = g_objects[gfA.objectIndex];
            auto &fcA = spA.faces[gfA.faceIndex];
            Vector3 v1A = spA.vertices[fcA.v1];
            Vector3 v2A = spA.vertices[fcA.v2];
            Vector3 v3A = spA.vertices[fcA.v3];
            Vector3 centroidA = (v1A + v2A + v3A) * (1.0 / 3.0);
            auto &gfB = g_allFaces[b];
            auto &spB = g_objects[gfB.objectIndex];
            auto &fcB = spB.faces[gfB.faceIndex];
            Vector3 v1B = spB.vertices[fcB.v1];
            Vector3 v2B = spB.vertices[fcB.v2];
            Vector3 v3B = spB.vertices[fcB.v3];
            Vector3 centroidB = (v1B + v2B + v3B) * (1.0 / 3.0);
            if(axis == 0) return centroidA.x < centroidB.x;
            if(axis == 1) return centroidA.y < centroidB.y;
            return centroidA.z < centroidB.z;
        };
        std::sort(indices.begin() + start, indices.begin() + end, comparator);
        int mid = start + count / 2;
        node->left = buildBVH(indices, start, mid);
        node->right = buildBVH(indices, mid, end);
    }
    return node;
}

BVHNode* buildBVH(std::vector<int>& indices) {
    return buildBVH(indices, 0, indices.size());
}

void deleteBVH(BVHNode* node) {
    if (!node) return;
    deleteBVH(node->left);
    deleteBVH(node->right);
    delete node;
}

BVHNode* g_bvhRoot = nullptr;

std::pair<int,int> traverseBVH(BVHNode* node, const Ray &ray, double &closest_t) {
    std::pair<int,int> best(-1,-1);
    if (!node) return best;
    if (!node->box.intersect(ray, 0.0, closest_t))
        return best;
    if (node->left == nullptr && node->right == nullptr) {
        for (int idx : node->faceIndices) {
            auto &gf = g_allFaces[idx];
            auto &sp = g_objects[gf.objectIndex];
            auto &fc = sp.faces[gf.faceIndex];
            Vector3 v1 = sp.vertices[fc.v1];
            Vector3 v2 = sp.vertices[fc.v2];
            Vector3 v3 = sp.vertices[fc.v3];
            Vector3 e1 = v2 - v1;
            Vector3 e2 = v3 - v1;
            Vector3 pvec = cross(ray.dir, e2);
            double det = dot(e1, pvec);
            if (std::fabs(det) < 1e-14) continue;
            double invDet = 1.0 / det;
            Vector3 tvec = ray.origin - v1;
            double u = dot(tvec, pvec) * invDet;
            if (u < 0.0 || u > 1.0) continue;
            Vector3 qvec = cross(tvec, e1);
            double v = dot(ray.dir, qvec) * invDet;
            if (v < 0.0 || (u+v) > 1.0) continue;
            double t = dot(e2, qvec) * invDet;
            if (t > 1e-14 && t < closest_t) {
                closest_t = t;
                best = { gf.objectIndex, gf.faceIndex };
            }
        }
        return best;
    } else {
        std::pair<int,int> best_left = traverseBVH(node->left, ray, closest_t);
        std::pair<int,int> best_right = traverseBVH(node->right, ray, closest_t);
        return (best_left.first != -1) ? best_left : best_right;
    }
}

bool traverseBVH_occlusion(BVHNode* node, const Ray &ray, double max_t, int iFace, int jFace) {
    if (!node) return false;
    if (!node->box.intersect(ray, 0.0, max_t))
        return false;
    if (node->left == nullptr && node->right == nullptr) {
        for (int idx : node->faceIndices) {
            if (idx == iFace || idx == jFace) continue;
            auto &gf = g_allFaces[idx];
            auto &sp = g_objects[gf.objectIndex];
            auto &fc = sp.faces[gf.faceIndex];
            Vector3 v1 = sp.vertices[fc.v1];
            Vector3 v2 = sp.vertices[fc.v2];
            Vector3 v3 = sp.vertices[fc.v3];
            Vector3 e1 = v2 - v1;
            Vector3 e2 = v3 - v1;
            Vector3 pvec = cross(ray.dir, e2);
            double det = dot(e1, pvec);
            if (std::fabs(det) < 1e-14) continue;
            double invDet = 1.0 / det;
            Vector3 tvec = ray.origin - v1;
            double u = dot(tvec, pvec) * invDet;
            if (u < 0.0 || u > 1.0) continue;
            Vector3 qvec = cross(tvec, e1);
            double v = dot(ray.dir, qvec) * invDet;
            if (v < 0.0 || (u+v) > 1.0) continue;
            double t = dot(e2, qvec) * invDet;
            if (t > 1e-10 && t < max_t - 1e-10)
                return true;
        }
        return false;
    } else {
        if (traverseBVH_occlusion(node->left, ray, max_t, iFace, jFace))
            return true;
        if (traverseBVH_occlusion(node->right, ray, max_t, iFace, jFace))
            return true;
        return false;
    }
}

bool isOccluded(const Vector3 &start, const Vector3 &end, int iFace, int jFace) {
    Vector3 dir = end - start;
    double dist = length(dir);
    if(dist < 1e-14)
         return false;
    dir = { dir.x / dist, dir.y / dist, dir.z / dist };
    Ray ray { start, dir };
    return traverseBVH_occlusion(g_bvhRoot, ray, dist, iFace, jFace);
}

struct FaceGeom {
    int objectIndex, faceIndex;
    Vector3 centroid, N, t, b;
    std::array<double,NUM_BANDS> emissivity;
};
std::vector<FaceGeom> g_faceGeom;

void buildFaceGeomCache() {
    g_faceGeom.clear();
    g_faceGeom.reserve(g_allFaces.size());
    for (auto &gf : g_allFaces) {
        const auto &sp = g_objects[gf.objectIndex];
        const auto &fc = sp.faces[gf.faceIndex];
        FaceGeom c;
        c.objectIndex = gf.objectIndex;
        c.faceIndex   = gf.faceIndex;

        auto &v1 = sp.vertices[fc.v1];
        auto &v2 = sp.vertices[fc.v2];
        auto &v3 = sp.vertices[fc.v3];
        c.centroid = {(v1.x+v2.x+v3.x)/3.0,
                      (v1.y+v2.y+v3.y)/3.0,
                      (v1.z+v2.z+v3.z)/3.0};

        c.N = fc.normal;
        Vector3 up = (std::fabs(c.N.x)<0.9)? Vector3{1,0,0} : Vector3{0,1,0};
        c.t = normalize(cross(c.N, up));
        c.b = cross(c.N, c.t);

        for (int band=0; band<NUM_BANDS; ++band)
            c.emissivity[band] = sp.emissivity[band];

        g_faceGeom.emplace_back(std::move(c));
    }
}

void logRadiationData(const std::string& method, int faceIndex, double oldTemp, double flux, double newTemp) {
    static std::mutex logMutex;
    static std::ofstream logFile("radiation_log.csv", std::ios::app);
    
    std::lock_guard<std::mutex> lock(logMutex);
    
    if (logFile.tellp() == 0) {
        logFile << "Method,FaceIndex,OldTemperature,RadiationFlux,NewTemperature\n";
    }
    
    logFile << method << ","
            << faceIndex << ","
            << std::fixed << std::setprecision(2) << oldTemp << ","
            << std::scientific << std::setprecision(6) << flux << ","
            << std::fixed << std::setprecision(2) << newTemp << "\n";
    logFile.flush();
}

void radiationMonteCarlo(double alpha, int samplesPerFace)
{
    int NF = static_cast<int>(g_allFaces.size());
    std::vector<double> newT(NF);

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < NF; i++) {
        std::mt19937 gen(12345u + static_cast<unsigned>(i));
        std::uniform_real_distribution<double> dist01(0.0, 1.0);

        auto &gf_i = g_allFaces[i];
        auto &sp_i = g_objects[gf_i.objectIndex];
        auto &fc_i = sp_i.faces[gf_i.faceIndex];

        double Ti = fc_i.temp;
        double E_i[NUM_BANDS];
        double sigmaT4_i = STEFAN_BOLTZMANN * std::pow(Ti, 4.0);
        for (int band = 0; band < NUM_BANDS; band++) {
            E_i[band] = sp_i.emissivity[band] * sigmaT4_i * getBlackbodyFraction(band, Ti);
        }

        Vector3 v1 = sp_i.vertices[fc_i.v1];
        Vector3 v2 = sp_i.vertices[fc_i.v2];
        Vector3 v3 = sp_i.vertices[fc_i.v3];
        Vector3 centroid = {(v1.x+v2.x+v3.x)/3.0,
                            (v1.y+v2.y+v3.y)/3.0,
                            (v1.z+v2.z+v3.z)/3.0};
        Vector3 Nf = fc_i.normal;
        Vector3 up = (std::fabs(Nf.x) < 0.9) ? Vector3{1,0,0} : Vector3{0,1,0};
        Vector3 t = normalize(cross(Nf, up));
        Vector3 b = cross(Nf, t);

        double flux_sum = 0.0;
        for (int s = 0; s < samplesPerFace; s++) {
            double r1 = std::sqrt(dist01(gen));
            double r2 = dist01(gen);
            double u = 1 - r1, v = r1 * (1 - r2), w = r1 * r2;
            Vector3 samplePoint = {
                v1.x*u + v2.x*v + v3.x*w,
                v1.y*u + v2.y*v + v3.y*w,
                v1.z*u + v2.z*v + v3.z*w
            };

            double phi      = 2.0 * M_PI * dist01(gen);
            double cosTheta = std::sqrt(dist01(gen));
            double sinTheta = std::sqrt(1.0 - cosTheta*cosTheta);
            Vector3 localDir = { sinTheta*std::cos(phi),
                                 sinTheta*std::sin(phi),
                                 cosTheta };
            Vector3 dir = {
                localDir.x*t.x + localDir.y*b.x + localDir.z*Nf.x,
                localDir.x*t.y + localDir.y*b.y + localDir.z*Nf.y,
                localDir.x*t.z + localDir.y*b.z + localDir.z*Nf.z
            };

            Ray ray{samplePoint, dir};
            auto hit = traceAllObjects(ray);
            if (hit.first < 0) continue;

            auto it = faceLookup.find(makeKey(hit.first, hit.second));
            if (it == faceLookup.end()) continue;
            int jGlobal = it->second;

            auto &gf_j = g_allFaces[jGlobal];
            auto &sp_j = g_objects[gf_j.objectIndex];
            auto &fc_j = sp_j.faces[gf_j.faceIndex];
            double Tj = fc_j.temp;
            double sigmaT4_j = STEFAN_BOLTZMANN * std::pow(Tj,4.0);
            double E_j[NUM_BANDS];
            for (int band = 0; band < NUM_BANDS; band++) {
                E_j[band] = sp_j.emissivity[band] * sigmaT4_j * getBlackbodyFraction(band, Tj);
            }

            double diffE = 0.0;
            for (int band = 0; band < NUM_BANDS; band++)
                diffE += (E_j[band] - E_i[band]);

            flux_sum += diffE;
        }

        double Q_i = M_PI * flux_sum / double(samplesPerFace);
        newT[i]   = Ti + alpha * Q_i;
        logRadiationData("MonteCarlo", i, Ti, Q_i, newT[i]);
    }

    #pragma omp parallel for
    for (int i = 0; i < NF; i++) {
        auto &gf = g_allFaces[i];
        g_objects[gf.objectIndex].faces[gf.faceIndex].temp = newT[i];
    }

    std::cout << "[Monte Carlo] Radiation step completed.\n";
}

void radiationDO(double alpha, int nTheta, int nPhi)
{
    static int prev_nTheta = 0, prev_nPhi = 0;
    if(nTheta != prev_nTheta || nPhi != prev_nPhi) {
        computeDONodes(nTheta, nPhi);
        prev_nTheta = nTheta;
        prev_nPhi = nPhi;
    }
    int NF = static_cast<int>(g_allFaces.size());
    std::vector<double> newT(NF, 300.0);
    
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < NF; i++){
        auto &gf_i = g_allFaces[i];
        auto &sp_i = g_objects[gf_i.objectIndex];
        auto &fc_i = sp_i.faces[gf_i.faceIndex];
        double Ti = fc_i.temp;
        double E_i[NUM_BANDS];
        for (int band = 0; band < NUM_BANDS; band++){
             double sigmaT4 = STEFAN_BOLTZMANN * std::pow(Ti, 4.0);
             E_i[band] = (1.0 - sp_i.reflectivity[band]) * sigmaT4 * getBlackbodyFraction(band, Ti) * sp_i.emissivity[band];
        }
        Vector3 v1 = sp_i.vertices[fc_i.v1];
        Vector3 v2 = sp_i.vertices[fc_i.v2];
        Vector3 v3 = sp_i.vertices[fc_i.v3];
        Vector3 centroid = { (v1.x+v2.x+v3.x)/3.0, (v1.y+v2.y+v3.y)/3.0, (v1.z+v2.z+v3.z)/3.0 };
        Vector3 Nf = fc_i.normal;
        Vector3 up = (std::fabs(Nf.x) < 0.9) ? Vector3{1,0,0} : Vector3{0,1,0};
        Vector3 t = normalize(cross(Nf, up));
        Vector3 b = cross(Nf, t);
        double Qi = 0.0;
        for (auto &node : g_DO_nodes) {
            Vector3 dir = { t.x * node.ldx + b.x * node.ldy + Nf.x * node.ldz,
                            t.y * node.ldx + b.y * node.ldy + Nf.y * node.ldz,
                            t.z * node.ldx + b.z * node.ldy + Nf.z * node.ldz };
            Ray ray { centroid, dir };
            auto hit = traceAllObjects(ray);
            if(hit.first < 0)
                continue;
            int jGlobal = -1;
            auto it = faceLookup.find(makeKey(hit.first, hit.second));
            if(it != faceLookup.end())
                jGlobal = it->second;
            else
                continue;
            auto &gf_j = g_allFaces[jGlobal];
            auto &sp_j = g_objects[gf_j.objectIndex];
            auto &fc_j = sp_j.faces[gf_j.faceIndex];
            double Tj = fc_j.temp;
            double E_j[NUM_BANDS];
            for (int band = 0; band < NUM_BANDS; band++){
                 double sigmaT4 = STEFAN_BOLTZMANN * std::pow(Tj, 4.0);
                 E_j[band] = (1.0 - sp_j.reflectivity[band]) * sigmaT4 * getBlackbodyFraction(band, Tj) * sp_j.emissivity[band];
            }
            double flux_dir = 0.0;
            for (int band = 0; band < NUM_BANDS; band++){
                 flux_dir += (E_j[band] - E_i[band]);
            }
            Qi += node.weight * flux_dir;
        }
        newT[i] = Ti + alpha * Qi;
        
        logRadiationData("DiscreteOrdinates", i, Ti, Qi, newT[i]);
    }
    for (int i = 0; i < NF; i++){
         g_objects[g_allFaces[i].objectIndex].faces[g_allFaces[i].faceIndex].temp = newT[i];
    }
    std::cout << "[DO] Radiation step completed.\n";
}

void radiationS2S(double alpha, int subDiv) {
    int NF = static_cast<int>(g_allFaces.size());
    std::vector<double> newT(NF, 300.0);

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < NF; i++) {
        auto &gf_i = g_allFaces[i];
        Object &sp_i = g_objects[gf_i.objectIndex];
        Face &fc_i = sp_i.faces[gf_i.faceIndex];
        double Ti = fc_i.temp;

        double E_i[NUM_BANDS] = {0.0};
        for (int band = 0; band < NUM_BANDS; band++) {
            double sigmaT4 = STEFAN_BOLTZMANN * std::pow(Ti, 4.0);
            E_i[band] = (1.0 - sp_i.reflectivity[band]) *
                        sp_i.emissivity[band] *
                        sigmaT4 *
                        getBlackbodyFraction(band, Ti);
        }

        double netFlux = 0.0;
        for (int j = 0; j < NF; j++) {
            if (i == j) continue;
            
            auto &gf_j = g_allFaces[j];
            Object &sp_j = g_objects[gf_j.objectIndex];
            Face &fc_j = sp_j.faces[gf_j.faceIndex];
            double Tj = fc_j.temp;

            double E_j[NUM_BANDS] = {0.0};
            for (int band = 0; band < NUM_BANDS; band++) {
                double sigmaT4 = STEFAN_BOLTZMANN * std::pow(Tj, 4.0);
                E_j[band] = (1.0 - sp_j.reflectivity[band]) *
                            sp_j.emissivity[band] *
                            sigmaT4 *
                            getBlackbodyFraction(band, Tj);
            }

            double diffE = 0.0;
            for (int band = 0; band < NUM_BANDS; band++) {
                diffE += (E_j[band] - E_i[band]);
            }

            netFlux += g_viewFactors[i][j] * diffE;
        }

        double Q_i = netFlux / fc_i.area;
        double damping = 1.0 / (1.0 + std::abs(Q_i) * 0.1);
        newT[i] = Ti + alpha * Q_i * damping;

        if (newT[i] < 0.0) newT[i] = 0.0;
        if (newT[i] > 5000.0) newT[i] = 5000.0;

        logRadiationData("S2S", i, Ti, Q_i, newT[i]);
    }

    #pragma omp parallel for
    for (int i = 0; i < NF; i++) {
        g_objects[g_allFaces[i].objectIndex].faces[g_allFaces[i].faceIndex].temp = newT[i];
    }

    std::cout << "[S2S] Radiation step completed.\n";
}

void radiationFull(double alpha){
    buildBlackbodyFractionTable();
    switch(g_model) {
       case RadiationModel::MONTE_CARLO:
           radiationMonteCarlo(alpha, monteCarloSamples);
           break;
       case RadiationModel::DISCRETE_ORDINATES:
           radiationDO(alpha, discreteOrdinates_nTheta, discreteOrdinates_nPhi);
           break;
       case RadiationModel::S2S:
           radiationS2S(alpha, s2sSubDiv);
           break;
    }
}

void computeViewFactors_S2S(int subDiv) {
    int NF = static_cast<int>(g_allFaces.size());
    g_viewFactors.assign(NF, std::vector<double>(NF, 0.0));

    struct TriInfo {
        std::vector<Vector3> samples;
        Vector3 normal;
        double area;
    };
    std::vector<TriInfo> triInfos(NF);

    auto sampleTriangle = [&](const Vector3 &A, const Vector3 &B, const Vector3 &C, int subdiv) -> std::vector<Vector3> {
        std::vector<Vector3> pts;
        pts.reserve((subdiv * (subdiv + 1)) / 2);
        for (int i = 0; i < subdiv; i++) {
            for (int j = 0; j < subdiv - i; j++) {
                double u = static_cast<double>(i) / subdiv;
                double v = static_cast<double>(j) / subdiv;
                double w = 1.0 - u - v;
                if (w < 0.0) continue;
                pts.push_back({ A.x * u + B.x * v + C.x * w,
                                A.y * u + B.y * v + C.y * w,
                                A.z * u + B.z * v + C.z * w });
            }
        }
        return pts;
    };

    for (int i = 0; i < NF; i++) {
        auto &gf_i = g_allFaces[i];
        auto &sp_i = g_objects[gf_i.objectIndex];
        auto &fc_i = sp_i.faces[gf_i.faceIndex];
        
        Vector3 v1 = sp_i.vertices[fc_i.v1];
        Vector3 v2 = sp_i.vertices[fc_i.v2];
        Vector3 v3 = sp_i.vertices[fc_i.v3];
        
        triInfos[i].samples = sampleTriangle(v1, v2, v3, subDiv);
        triInfos[i].normal = fc_i.normal;
        triInfos[i].area = fc_i.area;
    }

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < NF; i++) {
        auto &tri_i = triInfos[i];
        for (int j = 0; j < NF; j++) {
            if (i == j) continue;
            
            auto &tri_j = triInfos[j];
            double Fij = 0.0;
            int valid_samples = 0;
            
            for (const auto &pi : tri_i.samples) {
                for (const auto &pj : tri_j.samples) {
                    Vector3 r = {pj.x - pi.x, pj.y - pi.y, pj.z - pi.z};
                    double r2 = r.x*r.x + r.y*r.y + r.z*r.z;
                    if (r2 < 1e-10) continue;
                    
                    double r_len = std::sqrt(r2);
                    Vector3 r_norm = {r.x/r_len, r.y/r_len, r.z/r_len};
                    
                    double cos_theta_i = dot(tri_i.normal, r_norm);
                    double cos_theta_j = dot(tri_j.normal, {-r_norm.x, -r_norm.y, -r_norm.z});
                    
                    if (cos_theta_i <= 0.0 || cos_theta_j <= 0.0) continue;
                    
                    Ray ray = {pi, r_norm};
                    auto hit = traceAllObjects(ray);
                    if (hit.first != g_allFaces[j].objectIndex || hit.second != g_allFaces[j].faceIndex) continue;
                    
                    Fij += cos_theta_i * cos_theta_j / (M_PI * r2);
                    valid_samples++;
                }
            }
            
            if (valid_samples > 0) {
                g_viewFactors[i][j] = Fij / valid_samples;
            }
        }
    }

    for (int i = 0; i < NF; i++) {
        double sum = 0.0;
        for (int j = 0; j < NF; j++) {
            sum += g_viewFactors[i][j];
        }
        if (sum > 0.0) {
            for (int j = 0; j < NF; j++) {
                g_viewFactors[i][j] /= sum;
            }
        }
    }

    std::cout << "[S2S] View factors computed and normalized.\n";
}