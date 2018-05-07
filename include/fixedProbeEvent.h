#ifndef FIXED_PROBE_EVENT_H
#define FIXED_PROBE_EVENT_H 

// a fixed probe data point that corresponds to 
// the average field at the given time 

typedef struct fixedProbeEvent{
   double time;
   double freq;
   double freqErr; 
} fixedProbeEvent_t; 

#endif 
