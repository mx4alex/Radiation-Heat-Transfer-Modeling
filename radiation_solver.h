#ifndef RADIATION_SOLVER_H
#define RADIATION_SOLVER_H

#include "surface_properties.h"

enum class RadiationModel {
    MONTE_CARLO,
    DISCRETE_ORDINATES,
    S2S
};

extern int monteCarloSamples;
extern int discreteOrdinates_nTheta;
extern int discreteOrdinates_nPhi;
extern int s2sSubDiv;

void computeAllViewFactors(RadiationModel model);
void radiationOneStep(double alpha);

#endif 
