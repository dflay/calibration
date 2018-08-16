#ifndef TROLLEY_ANA_EVENT_H
#define TROLLEY_ANA_EVENT_H

// a trolley analysis event

#define NUM_TRLY 17   

typedef struct trolleyAnaEvent{
   unsigned long long time[NUM_TRLY];  // time of measurement (ns since epoch)  
   Double_t freq[NUM_TRLY];            // frequency for each trolley probe (Hz)
   Double_t r[NUM_TRLY];               // radial position of probe
   Double_t y[NUM_TRLY];               // vertical position of probe
   Double_t phi[NUM_TRLY];             // azimuthal position of probe
   Double_t temp[NUM_TRLY];            // temperature (deg C)
} trolleyAnaEvent_t;

// trolley swap event 
typedef struct trolleySwapEvent { 
   Double_t time;  // sec since epoch  
   Double_t freq;
   Double_t freq_err;
   Double_t temp;
   Double_t temp_err;
   Double_t r;
   Double_t r_err;
   Double_t y;
   Double_t y_err;
   Double_t phi;
   Double_t phi_err;
} trolleySwapEvent_t; 

// trolley probe positions 

typedef struct trolleyProbePosition{   
   Double_t r[NUM_TRLY]; 
   Double_t dr[NUM_TRLY]; 
   Double_t y[NUM_TRLY]; 
   Double_t dy[NUM_TRLY]; 
   Double_t phi[NUM_TRLY]; 
   Double_t dphi[NUM_TRLY]; 
} trolleyProbePosition_t; 

// for Delta-B measurements 
typedef struct trolleyDeltaB{
   Double_t r[NUM_TRLY];
   Double_t y[NUM_TRLY];
   Double_t phi[NUM_TRLY];
   Double_t normQuad[NUM_TRLY]; 
   Double_t normQuad_err[NUM_TRLY];
   Double_t skewQuad[NUM_TRLY]; 
   Double_t skewQuad_err[NUM_TRLY];
   Double_t azi[NUM_TRLY]; 
   Double_t azi_err[NUM_TRLY];
   Int_t probeID[NUM_TRLY]; 
} trolleyDeltaB_t; 

#endif 
