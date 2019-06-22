#include "../include/CalibFuncs.h"
//______________________________________________________________________________
int GetOmegaP_err(perturbation_t pert,double &err){
   // error on delta_t  
   // double err_sq = pert.sigma_err*pert.sigma_err     + pert.chi_err*pert.chi_err 
   //               + pert.delta_m_err*pert.delta_m_err + pert.delta_eps_err*pert.delta_eps_err 
   //               + pert.delta_mag_err*pert.delta_mag_err;
   // UPDATE: accounts for rad damping, rotational effects, sag, oxygen, etc 
   double err_sq = pert.sigma_err*pert.sigma_err + pert.chi_err*pert.chi_err 
                 + pert.delta_t_err*pert.delta_t_err;
   err = TMath::Sqrt(err_sq)*0.06179;  // convert to Hz
   return 0;
}
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
   double delta_tot  = GetDeltaTerm(sigma,pert.delta_m,chi,eps,pert.delta_eps,pert.delta_mag); 
   // error on delta_t 
   double err=0;
   int rc = GetOmegaP_err(pert,err);  
   // calculate omega_p_free  
   freqFree    = freq/(1.-delta_tot);
   freqFreeErr = 0; // don't add in systematic uncertainties yet! TMath::Sqrt(err*err + freqErr*freqErr); 
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
   double delta_t  = GetDeltaTerm(sigma,pert.delta_m,chi,eps,pert.delta_eps,pert.delta_mag);
   // error on delta_t  
   double err=0;
   int rc = GetOmegaP_err(pert,err);  
   // get the measured frequency -- USING ABSOLUTE FREQUENCY  
   double freq[3]     = {pp.freq    ,pp.freq_fxpr    ,pp.freq_trly};
   double freq_err[3] = {pp.freq_err,pp.freq_fxpr_err,pp.freq_trly_err};
   // calculate omega_p_free  
   for(int i=0;i<3;i++){
      freq_free[i]     = freq[i]/(1.-delta_t);  
      freq_free_err[i] = 0; // don't add in systematic uncertainties yet! TMath::Sqrt(err*err + freqErr*freqErr); 
   }
   return 0;
}
//______________________________________________________________________________
double GetDeltaTerm(double sigma,double delta_m,double chi,double eps,double delta_eps,double delta_mag){
   // compute the delta_t term.   
   // in the published formula, these signs are all positive 
   // input file carries the sign of the perturbation as measured (delta_m, delta_eps, delta_mag).  
   // if these delta terms are negative, minus sign in GetOmegaP_free corrects that. 
   // the factor of 1E-9 converts to 'absolute' scale, since the input is in ppb.
   double delta_t = 1E-9*(sigma + delta_m + (eps-4.*TMath::Pi()/3.)*chi + delta_eps + delta_mag);
   return delta_t;
}
//______________________________________________________________________________
int GetDiamagneticShielding(double sigma,double dsigma,double T,double &SIG,double &ERR){
   // compute diamagnetic shielding with temperature dependence
   // input units: ppb; output units: ppb  
   SIG = sigma - 10.36*(T-25);
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
