#define _USE_MATH_DEFINES
#include "radiation_solver.h"
#include "spectral_solver.h" 
#include <random>
#include <cmath>
#include <iostream>
#include <thread>
#include <atomic>
#include <mutex>

#include "surface_properties.h"

int monteCarloSamples       = 1000;
int discreteOrdinates_nTheta= 8;
int discreteOrdinates_nPhi  = 16;
int s2sSubDiv               = 4;

static const double STEFAN_BOLTZMANN = 5.670374419e-8;

struct Ray {
    Vector3 origin;
    Vector3 dir;
};

static std::pair<int,int> traceAllSpheres(const Ray &ray){
    double tMin = 1e30;
    std::pair<int,int> best(-1, -1);

    for(int si=0; si<(int)g_spheres.size(); si++){
        auto &sp = g_spheres[si];
        for(int fi=0; fi<(int)sp.faces.size(); fi++){
            auto &fc = sp.faces[fi];
            Vector3 v1 = sp.vertices[fc.v1];
            Vector3 v2 = sp.vertices[fc.v2];
            Vector3 v3 = sp.vertices[fc.v3];

            Vector3 e1 = v2 - v1;
            Vector3 e2 = v3 - v1;
            Vector3 pvec = cross(ray.dir,e2);
            double det = dot(e1,pvec);
            if(std::fabs(det)<1e-14) continue;
            double invDet = 1.0/det;
            Vector3 tvec = ray.origin - v1;
            double u = dot(tvec,pvec)*invDet;
            if(u<0.0 || u>1.0) continue;
            Vector3 qvec = cross(tvec,e1);
            double v = dot(ray.dir,qvec)*invDet;
            if(v<0.0 || u+v>1.0) continue;
            double t = dot(e2,qvec)*invDet;
            if(t>1e-14 && t<tMin){
                tMin = t;
                best = {si, fi};
            }
        }
    }
    return best;
}

static void computeViewFactors_MonteCarlo(int samplesPerFace){
    int NF = (int)g_allFaces.size();
    g_viewFactors.assign(NF, std::vector<double>(NF,0.0));

    std::mt19937 gen(12345);
    std::uniform_real_distribution<double> dist01(0.0,1.0);

    #pragma omp parallel
    {
        std::vector<std::vector<double>> localFF(NF, std::vector<double>(NF,0.0));
        #pragma omp for
        for(int fi=0; fi<NF; fi++){
            auto &gf_i = g_allFaces[fi];
            auto &sp_i = g_spheres[gf_i.sphereIndex];
            auto &fc_i = sp_i.faces[gf_i.faceIndex];

            Vector3 A = sp_i.vertices[fc_i.v1];
            Vector3 B = sp_i.vertices[fc_i.v2];
            Vector3 C = sp_i.vertices[fc_i.v3];
            Vector3 Nf= fc_i.normal;

            for(int s=0; s<samplesPerFace; s++){
                double r1 = std::sqrt(dist01(gen));
                double r2 = dist01(gen);
                double uu = 1.0 - r1;
                double vv = r1*(1.0 - r2);
                double ww = r1*r2;
                Vector3 point = {
                    A.x*uu + B.x*vv + C.x*ww,
                    A.y*uu + B.y*vv + C.y*ww,
                    A.z*uu + B.z*vv + C.z*ww
                };

                double phi  = 2.0*M_PI*dist01(gen);
                double cosT = dist01(gen);
                double sinT = std::sqrt(1.0 - cosT*cosT);

                Vector3 up = (std::fabs(Nf.x)<0.9)? Vector3{1,0,0}: Vector3{0,1,0};
                Vector3 t  = normalize(cross(Nf, up));
                Vector3 b  = cross(Nf, t);

                Vector3 dir = {
                    t.x*(std::cos(phi)*sinT) + b.x*(std::sin(phi)*sinT) + Nf.x*cosT,
                    t.y*(std::cos(phi)*sinT) + b.y*(std::sin(phi)*sinT) + Nf.y*cosT,
                    t.z*(std::cos(phi)*sinT) + b.z*(std::sin(phi)*sinT) + Nf.z*cosT
                };

                Ray ray{ point, dir };
                auto hit = traceAllSpheres(ray);
                if(hit.first >= 0){
                    int jGlobal = -1;
                    for(int jj=0; jj<NF; jj++){
                        if(g_allFaces[jj].sphereIndex==hit.first &&
                           g_allFaces[jj].faceIndex  ==hit.second){
                            jGlobal = jj;
                            break;
                        }
                    }
                    if(jGlobal>=0 && jGlobal!=fi){
                        localFF[fi][jGlobal]+=1.0;
                    }
                }
            }
        }
        #pragma omp critical
        {
            for(int i=0; i<NF; i++){
                for(int j=0; j<NF; j++){
                    g_viewFactors[i][j]+= localFF[i][j];
                }
            }
        }
    }
    for(int i=0; i<NF; i++){
        double sum=0.0;
        for(int j=0; j<NF; j++){
            sum+= g_viewFactors[i][j];
        }
        if(sum>1e-14){
            for(int j=0; j<NF; j++){
                g_viewFactors[i][j]/= sum;
            }
        }
    }
    std::cout<<"[MonteCarlo] View factors computed.\n";
}

static void computeViewFactors_DiscreteOrdinates(int nTheta, int nPhi){
    int NF = (int)g_allFaces.size();
    g_viewFactors.assign(NF, std::vector<double>(NF,0.0));

    std::vector<Vector3> directions;
    std::vector<double>  weights;

    for(int it=0; it<nTheta; it++){
        double alphaT = (it+0.5)/double(nTheta);
        double theta = alphaT*(M_PI/2.0);
        for(int ip=0; ip<nPhi; ip++){
            double alphaP = (ip+0.5)/double(nPhi);
            double phi = alphaP*(2.0*M_PI);
            double sinT = std::sin(theta);
            double cosT = std::cos(theta);
            Vector3 dir = {sinT*std::cos(phi), sinT*std::sin(phi), cosT};
            double dOmega= 2.0*M_PI/(nTheta*nPhi);
            directions.push_back(dir);
            weights.push_back(dOmega);
        }
    }

    #pragma omp parallel
    {
        std::vector<std::vector<double>> localFF(NF, std::vector<double>(NF,0.0));
        #pragma omp for
        for(int fi=0; fi<NF; fi++){
            auto &gf_i = g_allFaces[fi];
            auto &sp_i = g_spheres[gf_i.sphereIndex];
            auto &fc_i = sp_i.faces[gf_i.faceIndex];

            Vector3 v1 = sp_i.vertices[fc_i.v1];
            Vector3 v2 = sp_i.vertices[fc_i.v2];
            Vector3 v3 = sp_i.vertices[fc_i.v3];
            Vector3 centroid = {
                (v1.x+v2.x+v3.x)/3.0,
                (v1.y+v2.y+v3.y)/3.0,
                (v1.z+v2.z+v3.z)/3.0
            };
            Vector3 Nf = fc_i.normal;

            Vector3 up = (std::fabs(Nf.x)<0.9)?Vector3{1,0,0}:Vector3{0,1,0};
            Vector3 t  = normalize(cross(Nf, up));
            Vector3 b  = cross(Nf,t);

            for(size_t id=0; id<directions.size(); id++){
                double w = weights[id];
                Vector3 dLoc = directions[id];
                Vector3 dir_global = {
                    t.x*dLoc.x + b.x*dLoc.y + Nf.x*dLoc.z,
                    t.y*dLoc.x + b.y*dLoc.y + Nf.y*dLoc.z,
                    t.z*dLoc.x + b.z*dLoc.y + Nf.z*dLoc.z
                };
                Ray ray{centroid, dir_global};
                auto hit = traceAllSpheres(ray);
                if(hit.first>=0){
                    int jGlobal=-1;
                    for(int jj=0; jj<NF; jj++){
                        if(g_allFaces[jj].sphereIndex==hit.first &&
                           g_allFaces[jj].faceIndex  ==hit.second){
                            jGlobal=jj;break;
                        }
                    }
                    if(jGlobal>=0 && jGlobal!=fi){
                        localFF[fi][jGlobal]+= w;
                    }
                }
            }
        }
        #pragma omp critical
        {
            for(int i=0; i<NF; i++){
                for(int j=0; j<NF; j++){
                    g_viewFactors[i][j]+= localFF[i][j];
                }
            }
        }
    }

    for(int i=0; i<NF; i++){
        double sum=0.0;
        for(int j=0; j<NF; j++){
            sum+= g_viewFactors[i][j];
        }
        if(sum>1e-14){
            for(int j=0; j<NF; j++){
                g_viewFactors[i][j]/= sum;
            }
        }
    }
    std::cout<<"[DO] View factors computed.\n";
}

static void computeViewFactors_S2S(int subDiv){
    int NF = (int)g_allFaces.size();
    g_viewFactors.assign(NF, std::vector<double>(NF,0.0));

    auto sampleTriangle = [&](const Vector3&A, const Vector3&B, const Vector3&C, int subdiv){
        std::vector<Vector3> pts;
        pts.reserve(subdiv*subdiv);
        for(int i=0; i<subdiv; i++){
            for(int j=0; j<subdiv-i; j++){
                double u = double(i)/subdiv;
                double v = double(j)/subdiv;
                double w = 1.0-u-v;
                if(w<0.0) continue;
                Vector3 p={
                    A.x*u + B.x*v + C.x*w,
                    A.y*u + B.y*v + C.y*w,
                    A.z*u + B.z*v + C.z*w
                };
                pts.push_back(p);
            }
        }
        return pts;
    };

    struct TriInfo {
        std::vector<Vector3> samples;
        Vector3 normal;
        double area;
    };
    std::vector<TriInfo> triInfos(NF);
    for(int i=0; i<NF; i++){
        auto &gf = g_allFaces[i];
        auto &sp = g_spheres[gf.sphereIndex];
        auto &fc = sp.faces[gf.faceIndex];
        Vector3 A = sp.vertices[fc.v1];
        Vector3 B = sp.vertices[fc.v2];
        Vector3 C = sp.vertices[fc.v3];

        triInfos[i].samples = sampleTriangle(A,B,C,subDiv);
        triInfos[i].normal  = fc.normal;
        triInfos[i].area    = fc.area;
    }

    #pragma omp parallel for schedule(dynamic)
    for(int iFace=0; iFace<NF; iFace++){
        auto &ti_i = triInfos[iFace];
        double Ai  = ti_i.area;
        Vector3 Ni = ti_i.normal;
        for(int jFace=0; jFace<NF; jFace++){
            if(jFace==iFace) continue;
            auto &ti_j = triInfos[jFace];
            double Aj  = ti_j.area;
            Vector3 Nj = ti_j.normal;
            double sum_ij=0.0;
            for(auto &pi: ti_i.samples){
                for(auto &pj: ti_j.samples){
                    Vector3 rij = pj - pi;
                    double r2 = dot(rij,rij);
                    double r  = std::sqrt(r2);
                    double cos_i = dot(Ni,rij)/r;
                    double cos_j = -dot(Nj,rij)/r;
                    if(cos_i>0.0 && cos_j>0.0 && r>1e-14){
                        sum_ij+=(cos_i*cos_j)/(M_PI*r2);
                    }
                }
            }
            double subCount = (subDiv*(subDiv+1))/2.0;
            double subArea_i= Ai/subCount;
            double subArea_j= Aj/subCount;
            g_viewFactors[iFace][jFace]= sum_ij* subArea_i*subArea_j / Ai;
        }
    }

    for(int i=0; i<NF; i++){
        double sum=0.0;
        for(int j=0; j<NF; j++){
            sum+= g_viewFactors[i][j];
        }
        if(sum>1e-14){
            for(int j=0; j<NF; j++){
                g_viewFactors[i][j]/= sum;
            }
        }
    }
    std::cout<<"[S2S] View factors computed.\n";
}

void computeAllViewFactors(RadiationModel model){
    switch(model){
    case RadiationModel::MONTE_CARLO:
        computeViewFactors_MonteCarlo(monteCarloSamples);
        break;
    case RadiationModel::DISCRETE_ORDINATES:
        computeViewFactors_DiscreteOrdinates(discreteOrdinates_nTheta, discreteOrdinates_nPhi);
        break;
    case RadiationModel::S2S:
        computeViewFactors_S2S(s2sSubDiv);
        break;
    }
}

void radiationOneStep(double alpha){
    int NF = (int)g_allFaces.size();
    std::vector<double> newT(NF, 300.0);

    buildBlackbodyFractionTable();

    #pragma omp parallel for
    for(int i=0; i<NF; i++){
        auto &gf_i = g_allFaces[i];
        auto &sp_i = g_spheres[gf_i.sphereIndex];
        auto &fc_i = sp_i.faces[gf_i.faceIndex];

        double Ti = fc_i.temp;
        double Ti4= Ti*Ti*Ti*Ti;
        double net=0.0;

        for(int j=0; j<NF; j++){
            if(j==i) continue;
            auto &gf_j = g_allFaces[j];
            auto &sp_j = g_spheres[gf_j.sphereIndex];
            auto &fc_j = sp_j.faces[gf_j.faceIndex];

            double Tj = fc_j.temp;
            double Tj4= Tj*Tj*Tj*Tj;

            for(int band=0; band<3; band++){
                double eps_i = sp_i.emissivity[band];
                double rho_i = sp_i.reflectivity[band];
                double Ai = 1.0 - rho_i;

                double eps_j = sp_j.emissivity[band];
                double rho_j = sp_j.reflectivity[band];
                double Aj = 1.0 - rho_j;

                double Bi = STEFAN_BOLTZMANN*Ti4* getBlackbodyFraction(band, Ti) * eps_i;
                double Bj = STEFAN_BOLTZMANN*Tj4* getBlackbodyFraction(band, Tj) * eps_j;

                net += g_viewFactors[i][j]*(Aj*Bj - Ai*Bi);
            }
        }
        newT[i] = Ti + alpha*net;
    }

    #pragma omp parallel for
    for(int i=0; i<NF; i++){
        g_spheres[g_allFaces[i].sphereIndex].faces[g_allFaces[i].faceIndex].temp = newT[i];
    }
}
