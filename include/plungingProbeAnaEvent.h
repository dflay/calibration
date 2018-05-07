#ifndef PLUNGING_PROBE_ANA_EVENT_H 
#define PLUNGING_PROBE_ANA_EVENT_H

// a plunging probe event where we group the N NMR-DAQ events 
// into a single "analysis" event. We'll never have more than 25 shots 
// per NMR-DAQ run, so we'll cap the arrays at 25 

typedef struct plungingProbeAnaEvent{
   unsigned long long time[25];     // average time of the event  
   double r[25];                    // radial coordinate
   double y[25];                    // vertical coordinate 
   double phi[25];                  // azimuthal coordinate 
   double temp[25];                 // average temperature (deg C) 
   double temp_err[25];             // temperature uncertainty (deg C) 
   double freq[25];                 // average frequency (Hz) 
   double freq_err[25];             // frequency uncertainty (Hz)
   double freq_LO[25];              // local oscillator (Hz) 
   double freq_RF[25];              // pi/2 frequency (Hz)
   int run;                         // NMR-DAQ run   
   int numTraces;                   // number of traces for the run 
} plungingProbeAnaEvent_t;  

#endif  
