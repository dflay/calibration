#ifndef NMR_MEAS_H
#define NMR_MEAS_H

// a struct for nmr measurements 

typedef struct nmr_meas{
   std::string name;
   double freq;                // extracted frequency (no drift corr) 
   double freq_err;            // shot error  
   double freq_fxpr;           // extracted frequency (fxpr drift corr) 
   double freq_fxpr_err;       // shot error 
   double freq_trly;           // extracted frequency (trly drift corr) 
   double freq_trly_err;       // shot error 
   double p2p_err;             // point-to-point drift correction error 
   double r2r_err;             // run-to-run drift correction error  
   double T;                   // temperature  
} nmr_meas_t;

// a calibration multiswap event 

typedef struct calibSwap{
   double time;                // UTC time stamp 
   double freq;                // frequency (Hz) 
   double freqErr;             // frequency uncertainty (Hz) 
   double temp;                // temperature (deg C) 
   double tempErr;             // temperature uncertainty (deg C)
   double r;
   double rErr; 
   double y;
   double yErr; 
   double phi;
   double phiErr; 
   int type;                   // 0 = PP, 1 = TRLY 
} calibSwap_t; 

#endif 
