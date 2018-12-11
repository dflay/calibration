#ifndef RESULT_H
#define RESULT_H

// a data struct to store the results 

typedef struct result{
   // PP (raw) 
   double ppRaw;
   double ppRaw_err;
   double ppRaw_fxpr;
   double ppRaw_fxpr_err;
   double ppRaw_trly;
   double ppRaw_trly_err;
   // PP (free) 
   double ppFree;
   double ppFree_err;
   double ppFree_fxpr;
   double ppFree_fxpr_err;
   double ppFree_trly;
   double ppFree_trly_err;
   // TRLY  
   double trly;
   double trly_err;
   double trly_fxpr;
   double trly_fxpr_err;
   double trly_trly;
   double trly_trly_err;
   // R2R drift correction   
   double driftShim;
   double driftShim_err;
   double driftShim_fxpr;
   double driftShim_fxpr_err;
   double driftShim_trly;
   double driftShim_trly_err;
   // PP-TRLY+driftCor   
   double diff;
   double diff_err;
   double diff_fxpr;
   double diff_fxpr_err;
   double diff_trly;
   double diff_trly_err;
} result_t; 

typedef struct result_prod{
   double diff;
   double diffFree;
   double diffErr;         // shot error 
   double mErr;            // misalignment error 
   double pErr;            // free-proton error  
   double diff_aba;
   double diffFree_aba;
   double diffErr_aba;
   double mErr_aba;
   double pErr_aba;       
} result_prod_t; 

// fill this struct and write to a ROOT file; 
// this is for a single trolley probe. 
// after reading in from the ROOT file,
// we should have a vector of size 17  
typedef struct calib_result { 
   double calibCoeff;             // [raw] calibration coefficient (PP-TRLY) 
   double calibCoeffErr;          // [raw] shot uncertainty for calib coeff 
   double calibCoeff_aba;         // [ABA] calib coeff
   double calibCoeffErr_aba;      // [ABA] shot uncertainty  
   double calibCoeffFree;         // [raw] free-proton calib coeff  
   double calibCoeffFreeErr;      // [raw] shot uncertainty  
   double calibCoeffFree_aba;     // [ABA] free-proton calib coeff   
   double calibCoeffFreeErr_aba;  // [ABA] shot uncertainty  
   double freeErr;                // error from proton corrections (added in quadrature)  
   double dr;                     // misalignment error: combines (x,y,z) errors added in quadrature 
   double deltaB_x;               // deltaB(PP-TRLY) for x gradient 
   double deltaB_y;               // deltaB(PP-TRLY) for y gradient
   double deltaB_z;               // deltaB(PP-TRLY) for z gradient 
   double dBdx_imp;               // imposed gradient along x (radial)    
   double dBdy_imp;               // imposed gradient along y (vertical)  
   double dBdz_imp;               // imposed gradient along z (azimuthal) 
   double dBdx_shim;              // shimmed gradient along x (radial)    
   double dBdy_shim;              // shimmed gradient along y (vertical)  
   double dBdz_shim;              // shimmed gradient along z (azimuthal) 
} calib_result_t; 

#endif 
