module conduction_solver;

import surface_properties;

#include <vector>
#include <cmath>
#include <algorithm>

static const double TMIN_LIMIT=0.0;
static const double TMAX_LIMIT=5000.0;

void conductionOneStepAll(double kappa, double dt) {
    int NF = static_cast<int>(g_allFaces.size());
    std::vector<double> deltaEnergy(NF, 0.0);

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < NF; i++) {
        auto &gf_i = g_allFaces[i];
        Object &sp_i = g_objects[gf_i.objectIndex];
        Face &fc_i = sp_i.faces[gf_i.faceIndex];
        double Ti = fc_i.temp;

        int neighborCount = static_cast<int>(fc_i.neighbors.size());
        if (neighborCount == 0)
            continue;
        double A_interface = fc_i.area / neighborCount;

        for (int n : fc_i.neighbors) {
            auto &gf_j = g_allFaces[n];
            Object &sp_j = g_objects[gf_j.objectIndex];
            Face &fc_j = sp_j.faces[gf_j.faceIndex];
            double Tj = fc_j.temp;
            double dQ = kappa * (Tj - Ti) * A_interface;
            deltaEnergy[i] += dQ;
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < NF; i++) {
        auto &gf_i = g_allFaces[i];
        Object &sp_i = g_objects[gf_i.objectIndex];
        Face &fc_i = sp_i.faces[gf_i.faceIndex];

        double m_face = (sp_i.totalArea > 0.0) ? (fc_i.area / sp_i.totalArea) * sp_i.mass : sp_i.mass;
        double dT = deltaEnergy[i] * dt / (m_face * sp_i.cp);
        double T_new = fc_i.temp + dT;
        fc_i.temp = std::max(0.0, std::min(5000.0, T_new));
    }
}
void applyBoundaryConditions(double dtFactor){
    for (int si = 0; si < static_cast<int>(g_objects.size()); si++){
        auto &object = g_objects[si];
        if (object.boundaryCondition == BoundaryConditionType::Second){
            double flux = object.boundaryValue;
            double dT = (flux * dtFactor * object.totalArea) / (object.mass * object.cp);
            for (auto &fc : object.faces){
                fc.temp += dT;
                fc.temp = std::max(0.0, std::min(5000.0, fc.temp));
            }
        }
    }
}
