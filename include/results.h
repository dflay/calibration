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

#endif 
