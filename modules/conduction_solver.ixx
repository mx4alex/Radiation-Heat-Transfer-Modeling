module;

import surface_properties;

export module conduction_solver;

export void conductionOneStepAll(double kappa, double dt);
export void applyBoundaryConditions(double dtFactor); 