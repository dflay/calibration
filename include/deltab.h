#ifndef DELTAB_H
#define DELTAB_H

// a struct to store deltaB results 

typedef struct deltab { 
   std::string name;
   double dB;                // no drift correction
   double dB_err;            // shot error 
   double dB_fxpr;           // drift corrected using FXPR  
   double dB_fxpr_err;       // shot error 
   double dB_trly;           // drift corrected using TRLY 
   double dB_trly_err;       // shot error 
   double drift_fxpr;        // drift correction using FXPR (mostly just R2R) 
   double drift_fxpr_err;    // drift error using fxpr (P2P + R2R) 
   double drift_trly;        // drift correction using TRLY (mostly just R2R) 
   double drift_trly_err;    // drift error using trly (P2P + R2R)
   double p2p_err;           // point-to-point (P2P) drift error  
   double r2r_err;           // run-to-run (R2R) drift error 
} deltab_t; 

// another struct
// index 0 = x, 1 = y, z = 2 
typedef struct deltab_prod {
   double dB[3];
   double dB_err[3];
   double dB_aba[3];
   double dB_aba_err[3];
   double dB_opt[3];
   double dB_opt_err[3];
} deltab_prod_t; 

#endif 
