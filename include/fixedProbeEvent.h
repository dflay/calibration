#ifndef FIXED_PROBE_EVENT_H
#define FIXED_PROBE_EVENT_H 

// a fixed probe data point that corresponds to 
// the average field at the given time 
typedef struct averageFixedProbeEvent {
   unsigned long long time;
   double freq;
   double freqErr; 
} averageFixedProbeEvent_t; 

// a data struct that contains the probe's physical location

typedef struct fixedProbeEvent {
   unsigned long long time;  // UTC timestamp 
   double freq;              // frequency (Hz)  
   int aziID;                // 1, 2, 3,...
   int probeID;              // 0--377
   int stationID;            // 0--71 
   char yokeID;              // A, B, C,...
   char layerID;             // top (T), bottom (B)
   char radID;               // inner (I), middle (M), outer (O) 
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
