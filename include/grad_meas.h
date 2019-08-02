#ifndef GRAD_MEAS_H
#define GRAD_MEAS_H 

// a gradient measurement

typedef struct grad_meas {
   std::string name;         // rad, vert, or azi 
   double grad;              // measured gradient 
   double grad_err;
   double grad_fxpr;         // drift corrected (fxpr) 
   double grad_fxpr_err;    
   double grad_trly;         // drift corrected (trly) 
   double grad_trly_err;
   double drift_fxpr;        // drift corrected (fxpr) 
   double drift_fxpr_err;    
   double drift_trly;        // drift corrected (trly) 
   double drift_trly_err;
   double p2p_err;           // point-to-point drift correction error 
   double r2r_err;           // run-to-run drift correction error 
} grad_meas_t;

// imposed gradient 

typedef struct imposed_gradient{
   double pos;
   double grad;
   double grad_err;
} imposed_gradient_t; 

// generic gradient 

typedef struct gradMeas {
   std::string name;
   double val;
   double val_err;
   double val_drift_cor; 
   double val_drift_cor_err; 
} gradMeas_t; 

#endif  
