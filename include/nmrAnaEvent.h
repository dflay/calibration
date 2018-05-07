#ifndef NMR_ANA_EVENT_H
#define NMR_ANA_EVENT_H

// a data structure for NMR-ANA results
#define NUM_ANA_METHODS 6  

enum nmrAnaType{
   kMidpoint          = 0,
   kLinearInterp      = 1,
   kLeastSquares      = 2,
   kMidpointPhase     = 3,
   kLinearInterpPhase = 4,
   kLeastSquaresPhase = 5
};

typedef struct nmrAnaEvent{
   double ampl;                  // Maximum amplitude (V) 
   double noise;                 // RMS noise (V) 
   double freq[NUM_ANA_METHODS]; // extracted zero-crossing frequency (midpoint, linear interp, least-squares, then phase versions) 
   double nc;                    // number of cycles 
   int run;                      // NMR-DAQ run 
   int pulse;                    // pulse number
   int zc;                       // number of zero crossings
} nmrAnaEvent_t; 

#endif
