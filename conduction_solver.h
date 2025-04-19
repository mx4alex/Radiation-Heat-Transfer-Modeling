#ifndef CONDUCTION_SOLVER_H
#define CONDUCTION_SOLVER_H

#include "surface_properties.h"

void conductionOneStepAll(double kappa);
void applyBoundaryConditions(double dtFactor);

#endif 
