module;

export module spectral_constants;

export static const double h  = 6.62607015e-34;
export static const double c_ = 2.99792458e8; 
export static const double kB = 1.380649e-23;
export static const double STEFAN_BOLTZMANN = 5.670374419e-8;

export static const double T_MIN = 0.0;
export static const double T_MAX = 5000.0;
export static const int N_T_POINTS = 100;

export static const int NUM_BANDS = 3;
export static const double bandEdges[NUM_BANDS + 1] = { 0.1e-6, 2.0e-6, 5.0e-6, 20.0e-6 }; 