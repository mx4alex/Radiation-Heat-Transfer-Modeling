#include "conduction_solver.h"
#include <algorithm>
#include "surface_properties.h" 

static const double TMIN_LIMIT=200.0;
static const double TMAX_LIMIT=2000.0;

void conductionOneStepAll(double kappa){
    int NF = (int)g_allFaces.size();
    std::vector<double> deltaT(NF,0.0);

    #pragma omp parallel for
    for(int i=0; i<NF; i++){
        auto &gf_i = g_allFaces[i];
        auto &sp_i = g_spheres[gf_i.sphereIndex];
        auto &fc_i = sp_i.faces[gf_i.faceIndex];

        double Ti = fc_i.temp;
        for(int n : fc_i.neighbors){
            auto &gf_j = g_allFaces[n];
            auto &sp_j = g_spheres[gf_j.sphereIndex];
            auto &fc_j = sp_j.faces[gf_j.faceIndex];
            double Tj = fc_j.temp;

            double dQ = kappa*(Tj - Ti)*fc_i.area;
            #pragma omp atomic
            deltaT[i]+= dQ;
        }
    }

    for(int i=0; i<NF; i++){
        auto &gf_i = g_allFaces[i];
        auto &sp_i = g_spheres[gf_i.sphereIndex];
        auto &fc_i = sp_i.faces[gf_i.faceIndex];

        double T_old = fc_i.temp;
        double T_new = T_old + deltaT[i];
        if(T_new<TMIN_LIMIT) T_new=TMIN_LIMIT;
        if(T_new>TMAX_LIMIT) T_new=TMAX_LIMIT;
        fc_i.temp = T_new;
    }
}

void applyBoundaryConditions(double dtFactor){
    for(int si=0; si<(int)g_spheres.size(); si++){
        auto &sphere = g_spheres[si];
        if(sphere.boundaryCondition==BoundaryConditionType::Second){
            double flux = sphere.boundaryValue;
            double dT = (flux* dtFactor * sphere.totalArea)/(sphere.mass*sphere.cp);
            for(auto &fc : sphere.faces){
                double Tnew = fc.temp + dT;
                if(Tnew<TMIN_LIMIT) Tnew=TMIN_LIMIT;
                if(Tnew>TMAX_LIMIT) Tnew=TMAX_LIMIT;
                fc.temp = Tnew;
            }
        }
    }
}
