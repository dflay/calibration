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
   // raw 
   double diff;
   double diffFree;
   double diffErr;         // shot error 
   double mErr;            // misalignment error 
   double pErr;            // free-proton error 
   // ABA 
   double diff_aba;
   double diffFree_aba;
   double diffErr_aba;
   double mErr_aba;
   double pErr_aba;      
   // optimized (mix of raw and ABA for Delta-B values; affects misalignment error)  
   double diff_opt;
   double diffFree_opt;
   double diffErr_opt;
   double mErr_opt;
   double pErr_opt;      
   double systErr;  
} result_prod_t; 

// fill this struct and write to a ROOT file; 
// this is for a single trolley probe. 
// after reading in from the ROOT file,
// we should have a vector of size 17 
// all values are necessarily stored in Hz
// types of results 
// raw: straight difference of values for swaps or Delta-B 
// ABA: use ABA algorithm to remove drift 
// opt: combination of raw and ABA for if/when ABA is missing 
//      - Always use opt for Delta-B, misalignment results! (opt = ABA if ABA is available) 
typedef struct calib_result { 
   double calibCoeff;             // [raw] calibration coefficient (PP-TRLY) 
   double calibCoeffErr;          // [raw] shot uncertainty for calib coeff 
   double calibCoeff_aba;         // [ABA] calib coeff
   double calibCoeffErr_aba;      // [ABA] shot uncertainty  
   double calibCoeff_opt;         // [opt] calib coeff
   double calibCoeffErr_opt;      // [opt] shot uncertainty  

   double calibCoeffFree;         // [raw] free-proton calib coeff  
   double calibCoeffFreeErr;      // [raw] shot uncertainty  
   double calibCoeffFree_aba;     // [ABA] free-proton calib coeff   
   double calibCoeffFreeErr_aba;  // [ABA] shot uncertainty  
   double calibCoeffFree_opt;     // [opt] free-proton calib coeff   
   double calibCoeffFreeErr_opt;  // [opt] shot uncertainty  

   double freeErr;                // error from proton corrections (added in quadrature) 

   double systErr;                // the total systematic uncertainty from other sources (TRLY footprint, cuts, etc)  

   // misCor will always be opt results 
   double misCor;                 // misalignment correction along x axis in Hz        
   double misCor_err;             // misalignment correction uncertainty along x axis in Hz 

   double dx;                     // misalignment error along x axis in Hz 
   double dy;                     // misalignment error along y axis in Hz 
   double dz;                     // misalignment error along z axis in Hz 
   double dr;                     // misalignment error: combines (x,y,z) errors added in quadrature 

   double deltaB_pp_x;            // deltaB(PP) for x gradient 
   double deltaB_pp_y;            // deltaB(PP) for y gradient
   double deltaB_pp_z;            // deltaB(PP) for z gradient 

   double deltaB_pp_xErr;         // deltaB(PP) for x gradient [uncertainty]  
   double deltaB_pp_yErr;         // deltaB(PP) for y gradient [uncertainty] 
   double deltaB_pp_zErr;         // deltaB(PP) for z gradient [uncertainty]  

   double deltaB_tr_x;            // deltaB(TRLY) for x gradient 
   double deltaB_tr_y;            // deltaB(TRLY) for y gradient
   double deltaB_tr_z;            // deltaB(TRLY) for z gradient
 
   double deltaB_tr_xErr;         // deltaB(TRLY) for x gradient [uncertainty]  
   double deltaB_tr_yErr;         // deltaB(TRLY) for y gradient [uncertainty] 
   double deltaB_tr_zErr;         // deltaB(TRLY) for z gradient [uncertainty]  

   double dBdx_imp;               // imposed gradient along x (radial)    
   double dBdy_imp;               // imposed gradient along y (vertical)  
   double dBdz_imp;               // imposed gradient along z (azimuthal) 

   double dBdx_impErr;            // imposed gradient along x (radial)    [uncertainty]  
   double dBdy_impErr;            // imposed gradient along y (vertical)  [uncertainty]  
   double dBdz_impErr;            // imposed gradient along z (azimuthal) [uncertainty]  

   double dBdx_shim;              // shimmed gradient along x (radial)    
   double dBdy_shim;              // shimmed gradient along y (vertical)  
   double dBdz_shim;              // shimmed gradient along z (azimuthal) 

   double dBdx_shimErr;           // shimmed gradient along x (radial)    [uncertainty]  
   double dBdy_shimErr;           // shimmed gradient along y (vertical)  [uncertainty]  
   double dBdz_shimErr;           // shimmed gradient along z (azimuthal) [uncertainty]  

} calib_result_t;

const char * const calib_result_str = "calibCoeff/D:calibCoeffErr/D:calibCoeff_aba/D:calibCoeffErr_aba/D:calibCoeff_opt/D:calibCoeffErr_opt/D:calibCoeffFree/D:calibCoeffFreeErr/D:calibCoeffFree_aba/D:calibCoeffFreeErr_aba/D:calibCoeffFree_opt/D:calibCoeffFreeErr_opt/D:freeErr/D:misCor/D:misCor_err/D:systErr/D:dx/D:dy/D:dz/D:dr/D:deltaB_pp_x/D:deltaB_pp_y/D:deltaB_pp_z/D:deltaB_pp_xErr/D:deltaB_pp_yErr/D:deltaB_pp_zErr/D:deltaB_tr_x/D:deltaB_tr_y/D:deltaB_tr_z/D:deltaB_tr_xErr/D:deltaB_tr_yErr/D:deltaB_tr_zErr/D:dBdx_imp/D:dBdy_imp/D:dBdz_imp/D:dBdx_impErr/D:dBdy_impErr/D:dBdz_impErr/D:dBdx_shim/D:dBdy_shim/D:dBdz_shim/D:dBdx_shimErr/D:dBdy_shimErr/D:dBdz_shimErr/D"; 

#endif 
