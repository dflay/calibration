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

#endif 
