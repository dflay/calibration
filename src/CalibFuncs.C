#include "../include/CalibFuncs.h"
//______________________________________________________________________________
int GetOmegaP_free(perturbation_t pert,double freq,double freqErr,double temp,double tempErr,
                   double &freqFree,double &freqFreeErr){
   // calculate delta_t 
   double eps      = pert.eps; 
   // diamagnetic shielding
   double sigma=0,sigma_err=0;
   GetDiamagneticShielding(pert.sigma,pert.sigma_err,temp,sigma,sigma_err); 
   // magnetic sucseptibility  
   double chi=0;
   GetMagneticSusceptibility(pert.chi,temp,chi);
   // put it all together -- NOTE MINUS SIGNS! 
   double delta_t  = sigma*1E-9 - pert.delta_m*1E-9 
                   + (eps-4.*TMath::Pi()/3.)*chi*1E-9 - pert.delta_eps*1E-9 - pert.delta_mag*1E-9; 
   // error on delta_t  
   double err_sq   = pert.sigma_err*pert.sigma_err + pert.delta_m_err*pert.delta_m_err 
                   + pert.delta_eps_err*pert.delta_eps_err + pert.delta_mag_err*pert.delta_mag_err;
   double err      = TMath::Sqrt(err_sq)*0.06179;  // convert to Hz 
   // calculate omega_p_free  
   freqFree    = freq/(1.-delta_t); 
   freqFreeErr = TMath::Sqrt(err*err + freqErr*freqErr);  
   return 0;
}
//______________________________________________________________________________
int GetOmegaP_free(nmr_meas_t pp,perturbation_t pert,double *freq_free,double *freq_free_err){
   // calculate delta_t 
   double eps      = pert.eps; 
   // diamagnetic shielding
   double sigma=0,sigma_err=0;
   GetDiamagneticShielding(pert.sigma,pert.sigma_err,pp.T,sigma,sigma_err); 
   // magnetic sucseptibility  
   double chi=0;
   GetMagneticSusceptibility(pert.chi,pp.T,chi);
   // put it all together -- NOTE MINUS SIGNS! 
   double delta_t  = sigma*1E-9 - pert.delta_m*1E-9 
                   + (eps-4.*TMath::Pi()/3.)*chi*1E-9 - pert.delta_eps*1E-9 - pert.delta_mag*1E-9; 
   // error on delta_t  
   double err_sq   = pert.sigma_err*pert.sigma_err + pert.delta_m_err*pert.delta_m_err 
                   + pert.delta_eps_err*pert.delta_eps_err + pert.delta_mag_err*pert.delta_mag_err;
   double err      = TMath::Sqrt(err_sq)*0.06179;  // convert to Hz 
   // get the measured frequency -- USING ABSOLUTE FREQUENCY  
   double freq[3]     = {pp.freq    ,pp.freq_fxpr    ,pp.freq_trly};
   double freq_err[3] = {pp.freq_err,pp.freq_fxpr_err,pp.freq_trly_err};
   // calculate omega_p_free  
   for(int i=0;i<3;i++){
      freq_free[i]     = freq[i]/(1.-delta_t);  
      freq_free_err[i] = TMath::Sqrt(err*err + freq_err[i]*freq_err[i]);  
   }
   return 0;
}
//______________________________________________________________________________
int GetDiamagneticShielding(double sigma,double dsigma,double T,double &SIG,double &ERR){
   // compute diamagnetic shielding with temperature dependence
   // input units: ppb; output units: ppb  
   SIG = sigma + 10.36*(T-25);
   ERR = TMath::Sqrt(dsigma*dsigma + 0.30*0.30);   
   return 0; 
}
//______________________________________________________________________________
int GetMagneticSusceptibility(double chi,double T,double &CHI){
   // compute magnetic susceptibility when accounting for temperature dependence 
   double a[3] = {1.38810E-4,-1.2685E-7,8.09E-10};
   CHI = chi*( 1. + a[0]*(T-20.) + a[1]*TMath::Power(T-20.,2.) + a[2]*TMath::Power(T-20.,3.) );
   return 0;
}
