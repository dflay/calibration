#ifndef PERTURBATION_H
#define PERTURBATION_H

// a struct to keep track of calibration probe perturbation and correction terms  

typedef struct perturbaton{
   double sigma;                // diamagnetic shielding 
   double sigma_err;           
   double chi;                  // sample magnetic susceptibility 
   double chi_err;               
   double eps;                  // shape factor 
   double eps_err;         
   double delta_m;              // material effect
   double delta_m_err; 
   double delta_eps;            // shape asymmetry 
   double delta_eps_err;
   double delta_mag;            // magnetic image effect  
   double delta_mag_err;        
} perturbation_t; 

#endif 
