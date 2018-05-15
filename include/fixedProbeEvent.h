#ifndef FIXED_PROBE_EVENT_H
#define FIXED_PROBE_EVENT_H 

// a fixed probe data point that corresponds to 
// the average field at the given time 

typedef struct fixedProbeEvent{
   double time;
   double freq;
   double freqErr; 
} fixedProbeEvent_t; 

// a data struct that contains the probe's physical location

typedef struct fxpr{
   char yokeID;  // A, B, C,...
   char layerID; // top (T), bottom (B)
   char radID;   // inner (I), middle (M), outer (O) 
   int aziID;    // 1, 2, 3,...
   int probeID;  // 0--377
   double freq; 
} fxpr_t; 

#endif 
