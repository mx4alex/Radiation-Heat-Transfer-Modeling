#include "spectral_solver.h"

#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include "surface_properties.h" 

static const double STEFAN_BOLTZMANN = 5.670374419e-8;
static const double h  = 6.62607015e-34;
static const double c_ = 2.99792458e8;
static const double kB = 1.380649e-23;

static const int NUM_BANDS = 3;
static const double bandEdges[NUM_BANDS + 1] = {
    0.1e-6,
    2.0e-6,
    5.0e-6,
    20.0e-6
};

static const double T_MIN = 200.0;
static const double T_MAX = 2000.0;
static const int    N_T_POINTS = 100;

static std::vector<std::vector<double>> fractionTable;
static std::vector<double> Tgrid;

static double planckEmission_lambda(double lambda, double T){
    if(lambda <= 0.0 || T <= 0.0) return 0.0;
    double exponent = (h*c_)/(lambda*kB*T);
    if(exponent > 700.) return 0.0;
    double denom = std::exp(exponent) - 1.0;
    if(denom < 1e-30) return 0.0;
    double a = 2.*h*c_*c_;
    return a/(std::pow(lambda,5)*denom);
}

static double integratePlanck(double l1, double l2, double T, int N=10){
    if(l2 <= l1 || T <= 0.) return 0.0;
    double step = (l2 - l1)/double(N);
    double sum = 0.0;
    for(int i = 0; i <= N; i++){
        double x = l1 + i*step;
        double f = planckEmission_lambda(x, T);
        double w = 1.0;
        if(i>0 && i<N) w = (i%2==0)?2.0:4.0;
        sum += w*f;
    }
    sum *= (step/3.0);
    return sum;
}

static double integratePlanckTotal(double T){
    return integratePlanck(0.1e-6, 100e-6, T, 20);
}

static double computeFractionBand(int bandIndex, double T){
    if(bandIndex<0 || bandIndex>=NUM_BANDS) return 0.0;
    double l1 = bandEdges[bandIndex];
    double l2 = bandEdges[bandIndex+1];
    double valBand  = integratePlanck(l1, l2, T, 10);
    double valTotal = integratePlanckTotal(T);
    if(valTotal < 1e-30) return 0.0;
    return valBand / valTotal;
}

void buildBlackbodyFractionTable(){
    fractionTable.resize(N_T_POINTS);
    Tgrid.resize(N_T_POINTS);
    for(int i = 0; i < N_T_POINTS; i++){
        double alpha = double(i)/(N_T_POINTS-1);
        double Tval = T_MIN + alpha*(T_MAX - T_MIN);
        Tgrid[i] = Tval;
    }

    for(int i = 0; i < N_T_POINTS; i++){
        fractionTable[i].resize(NUM_BANDS);
        for(int band = 0; band < NUM_BANDS; band++){
            fractionTable[i][band] = computeFractionBand(band, Tgrid[i]);
        }
    }
}

static double interpolateFraction(int bandIndex, double T){
    if(T<=T_MIN) return fractionTable[0][bandIndex];
    if(T>=T_MAX) return fractionTable[N_T_POINTS-1][bandIndex];

    double scale = (T - T_MIN)/(T_MAX - T_MIN)*(N_T_POINTS-1);
    int iLow = (int)std::floor(scale);
    int iHigh= iLow+1;
    if(iHigh>=N_T_POINTS) iHigh = N_T_POINTS-1;
    double frac = scale - iLow;
    double vL = fractionTable[iLow][bandIndex];
    double vH = fractionTable[iHigh][bandIndex];
    return vL + (vH - vL)*frac;
}

double getBlackbodyFraction(int bandIndex, double T){
    return interpolateFraction(bandIndex,T);
}
