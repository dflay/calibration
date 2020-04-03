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
   double freeProtCor;    // free proton correction size (ppb) 
   double freeProtCorErr; // free proton correction size (ppb) uncertainty  
   double ppTemp;           // temperature in deg C 
   double ppTempErr;        // temperature uncertainty in deg C   
   double trTemp;           // temperature in deg C 
   double trTempErr;        // temperature uncertainty in deg C  
   // cor = misalignment corrected 
   // raw 
   double diff;
   double diffCor;   
   double diffErr;         // shot error 
   double diffCorErr;      // shot error 
   double mErr;            // misalignment error 
   double pErr;            // free-proton error 
   // ABA 
   double diff_aba;
   double diffCor_aba;
   double diffErr_aba;
   double diffCorErr_aba;
   double mErr_aba;
   double pErr_aba;      
   // optimized (mix of raw and ABA for Delta-B values; affects misalignment error) 
   double diff_opt;
   double diffCor_opt;
   double diffErr_opt;
   double diffCorErr_opt;
   double mErr_opt;
   double pErr_opt;      
   double systErr;  

   // old stuff
   // double diffFree;
   // double diffCorFree;
   // double diffFree_aba;
   // double diffCorFree_aba;
   // double diffFree_opt;
   // double diffCorFree_opt;

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
   // without free proton, without misalignment
   double calibCoeff;             // [raw] calibration coefficient (PP-TRLY) 
   double calibCoeffErr;          // [raw] shot uncertainty for calib coeff 
   double calibCoeff_aba;         // [ABA] calib coeff
   double calibCoeffErr_aba;      // [ABA] shot uncertainty  
   double calibCoeff_opt;         // [opt] calib coeff
   double calibCoeffErr_opt;      // [opt] shot uncertainty  

   // without free proton, with misalignment
   double calibCoeffCor;             // [raw] calibration coefficient (PP-TRLY) 
   double calibCoeffCorErr;          // [raw] shot uncertainty for calib coeff 
   double calibCoeffCor_aba;         // [ABA] calib coeff
   double calibCoeffCorErr_aba;      // [ABA] shot uncertainty  
   double calibCoeffCor_opt;         // [opt] calib coeff
   double calibCoeffCorErr_opt;      // [opt] shot uncertainty  

   // with free proton, without misalignment
   double calibCoeffFree;         // [raw] free-proton calib coeff  
   double calibCoeffFreeErr;      // [raw] shot uncertainty  
   double calibCoeffFree_aba;     // [ABA] free-proton calib coeff   
   double calibCoeffFreeErr_aba;  // [ABA] shot uncertainty  
   double calibCoeffFree_opt;     // [opt] free-proton calib coeff   
   double calibCoeffFreeErr_opt;  // [opt] shot uncertainty  

   // with free proton, with misalignment
   double calibCoeffCorFree;         // [raw] free-proton calib coeff  
   double calibCoeffCorFreeErr;      // [raw] shot uncertainty  
   double calibCoeffCorFree_aba;     // [ABA] free-proton calib coeff   
   double calibCoeffCorFreeErr_aba;  // [ABA] shot uncertainty  
   double calibCoeffCorFree_opt;     // [opt] free-proton calib coeff   
   double calibCoeffCorFreeErr_opt;  // [opt] shot uncertainty  

   double freeErr;                // error from proton corrections (added in quadrature) 

   double systErr;                // the total systematic uncertainty from other sources (TRLY footprint, cuts, etc)  

   // misCor will always be opt results 
   double misCor;                 // misalignment correction along x axis in Hz        
   double misCor_err;             // misalignment correction uncertainty along x axis in Hz 

   double dB_x;                   // misalignment error along x axis in Hz 
   double dB_y;                   // misalignment error along y axis in Hz 
   double dB_z;                   // misalignment error along z axis in Hz 
   double dB_r_tot;               // misalignment error: combines (x,y,z) errors added in quadrature 

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

   double shim_x_a;               // shimmed field x a parameter (x2)  
   double shim_x_aErr;            // shimmed field x a parameter (x2) [uncertainty]    
   double shim_x_b;               // shimmed field x b parameter (x1)   
   double shim_x_bErr;            // shimmed field x b parameter (x1) [uncertainty]     
   double shim_x_c;               // shimmed field x c parameter (x0)  
   double shim_x_cErr;            // shimmed field x c parameter (x0) [uncertainty]     
   
   double shim_y_a;               // shimmed field y a parameter (y2)   
   double shim_y_aErr;            // shimmed field y a parameter (y2)  [uncertainty]  
   double shim_y_b;               // shimmed field y b parameter (y1)    
   double shim_y_bErr;            // shimmed field y b parameter (y1)  [uncertainty]     
   double shim_y_c;               // shimmed field y c parameter (y0)                  
   double shim_y_cErr;            // shimmed field y c parameter (y0)  [uncertainty]      
 
   double shim_z_a;               // shimmed field z a parameter (z2)     
   double shim_z_aErr;            // shimmed field z a parameter (z2)  [uncertainty]      
   double shim_z_b;               // shimmed field z b parameter (z1)                 
   double shim_z_bErr;            // shimmed field z b parameter (z1)  [uncertainty]   
   double shim_z_c;               // shimmed field z c parameter (z0)                
   double shim_z_cErr;            // shimmed field z c parameter (z0)  [uncertainty]  

   double misCor_x;               // misalignment correction along x  
   double misCor_xErr;            // misalignment correction along x [uncertainty]   
   double misCor_y;               // misalignment correction along y  
   double misCor_yErr;            // misalignment correction along y [uncertainty]   
   double misCor_z;               // misalignment correction along z  
   double misCor_zErr;            // misalignment correction along z [uncertainty]   
   double misCor_z_bar;           // misalignment correction along z (barcode corrected) 
   double misCor_zErr_bar;        // misalignment correction along z (barcode corrected) [uncertainty]  

   double mis_x;                  // misalignment along x in mm  
   double mis_xErr;               // misalignment along x in mm [uncertainty]   
   double mis_y;                  // misalignment along y in mm  
   double mis_yErr;               // misalignment along y in mm [uncertainty]   
   double mis_z;                  // misalignment along z in mm  
   double mis_zErr;               // misalignment along z in mm [uncertainty]  
   double mis_z_bar;              // misalignment along z in mm (barcode corrected) 
   double mis_zErr_bar;           // misalignment along z in mm (barcode corrected) [uncertainty]  

   double fpCor;                  // free-proton correction in ppb 
   double fpCorErr;               // free-proton correction in ppb [uncertainty]  

   double ppTemp;                 // PP temperature (avg) in deg C  
   double ppTempErr;              // PP temperature (avg) in deg C [uncertainty]  

   double trTemp;                 // TRLY temperature (avg) in deg C  
   double trTempErr;              // TRLY temperature (avg) in deg C [uncertainty] 

} calib_result_t;

const char * const calib_result_str = "calibCoeff/D:calibCoeffErr/D:calibCoeff_aba/D:calibCoeffErr_aba/D:calibCoeff_opt/D:calibCoeffErr_opt/D:calibCoeffFree/D:calibCoeffFreeErr/D:calibCoeffFree_aba/D:calibCoeffFreeErr_aba/D:calibCoeffFree_opt/D:calibCoeffFreeErr_opt/D:calibCoeffCor/D:calibCoeffCorErr/D:calibCoeffCor_aba/D:calibCoeffCorErr_aba/D:calibCoeffCor_opt/D:calibCoeffCorErr_opt/D:calibCoeffCorFree/D:calibCoeffCorFreeErr/D:calibCoeffCorFree_aba/D:calibCoeffCorFreeErr_aba/D:calibCoeffCorFree_opt/D:calibCoeffCorFreeErr_opt/D:freeErr/D:misCor/D:misCor_err/D:systErr/D:dx/D:dy/D:dz/D:dr/D:deltaB_pp_x/D:deltaB_pp_y/D:deltaB_pp_z/D:deltaB_pp_xErr/D:deltaB_pp_yErr/D:deltaB_pp_zErr/D:deltaB_tr_x/D:deltaB_tr_y/D:deltaB_tr_z/D:deltaB_tr_xErr/D:deltaB_tr_yErr/D:deltaB_tr_zErr/D:dBdx_imp/D:dBdy_imp/D:dBdz_imp/D:dBdx_impErr/D:dBdy_impErr/D:dBdz_impErr/D:dBdx_shim/D:dBdy_shim/D:dBdz_shim/D:dBdx_shimErr/D:dBdy_shimErr/D:dBdz_shimErr/D:shim_x_a/D:shim_x_aErr/D:shim_x_b/D:shim_x_bErr/D:shim_x_c/D:shim_x_cErr/D:shim_y_a/D:shim_y_aErr/D:shim_y_b/D:shim_y_bErr/D:shim_y_c/D:shim_y_cErr/D:shim_z_a/D:shim_z_aErr/D:shim_z_b/D:shim_z_bErr/D:shim_z_c/D:shim_z_cErr/D:misCor_x/D:misCor_xErr/D:misCor_y/D:misCor_yErr/D:misCor_z/D:misCor_zErr/D:mis_x/D:mis_xErr/D:mis_y/D:mis_yErr/D:mis_z/D:mis_zErr/D:fpCor/D:fpCorErr/D:ppTemp/D:ppTempErr/D:trTemp/D:trTempErr/D"; 

#endif 
