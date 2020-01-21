#ifndef PLUNGING_PROBE_ANA_EVENT_H 
#define PLUNGING_PROBE_ANA_EVENT_H

// a plunging probe event where we group the N NMR-DAQ events 
// into a single "analysis" event.  

#define PP_MAX_EVENTS 5000

typedef struct plungingProbeAnaEvent{
   unsigned long long time[PP_MAX_EVENTS];     // average time of the event  
   double r[PP_MAX_EVENTS];                    // radial coordinate
   double y[PP_MAX_EVENTS];                    // vertical coordinate 
   double phi[PP_MAX_EVENTS];                  // azimuthal coordinate 
   double temp[PP_MAX_EVENTS];                 // average temperature (deg C) 
   double temp_err[PP_MAX_EVENTS];             // temperature uncertainty (deg C) 
   double ampl[PP_MAX_EVENTS];                 // signal amplitude 
   double freq[PP_MAX_EVENTS];                 // average frequency (Hz) 
   double freq_err[PP_MAX_EVENTS];             // frequency uncertainty (Hz)
   double freq_LO[PP_MAX_EVENTS];              // local oscillator (Hz) 
   double freq_RF[PP_MAX_EVENTS];              // pi/2 frequency (Hz)
   double t2Time[PP_MAX_EVENTS];               // T2 time (sec) 
   int run;                                    // NMR-DAQ run   
   int midasRun;                               // MIDAS run    
   int numTraces;                              // number of traces for the (NMR-DAQ) run 
   int nzc[PP_MAX_EVENTS];                     // number of zero crossings  
   int traceNumber[PP_MAX_EVENTS];             // signal trace number   
   int channelNumber[PP_MAX_EVENTS];           // signal channel number   
} plungingProbeAnaEvent_t;  

#endif  
