#ifndef NMR_ANA_EVENT_H
#define NMR_ANA_EVENT_H

// a data structure for NMR-ANA results
#define NUM_ANA_METHODS 6  

namespace plungingProbeAnalysis { 
   enum nmrAnaType{
      kMidpoint          = 0,
      kLinearInterp      = 1,
      kLeastSquares      = 2,
      kMidpointPhase     = 3,
      kLinearInterpPhase = 4,
      kLeastSquaresPhase = 5
   };
}

typedef struct nmrAnaEvent{
   unsigned long long time;      // UTC time (ns)  
   double temp;                  // temperature (deg C) 
   double freq_LO;               // LO freq (Hz) 
   double freq_pi2;              // pi/2 freq (Hz) 
   double ampl;                  // Maximum amplitude (V) 
   double noise;                 // RMS noise (V) 
   double t2;                    // T2 time of signal (sec) 
   double freq[NUM_ANA_METHODS]; // extracted zero-crossing frequency (midpoint, linear interp, least-squares, then phase versions) 
   double nc;                    // number of cycles 
   int run;                      // NMR-DAQ run 
   int pulse;                    // pulse number
   int ch;                       // channel number 
   int zc;                       // number of zero crossings
} nmrAnaEvent_t; 

#endif
