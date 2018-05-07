#ifndef TROLLEY_ANA_EVENT_H
#define TROLLEY_ANA_EVENT_H

// a trolley analysis event

#define NUM_TRLY 17   

typedef struct trolleyAnaEvent{
   unsigned long long time[NUM_TRLY];  // time of measurement 
   Double_t freq[NUM_TRLY];            // frequency for each trolley probe (Hz)
   Double_t r[NUM_TRLY];               // radial position of probe
   Double_t y[NUM_TRLY];               // vertical position of probe
   Double_t phi[NUM_TRLY];             // azimuthal position of probe
   Double_t temp[NUM_TRLY];            // temperature (deg C)
} trolleyAnaEvent_t;

#endif 
