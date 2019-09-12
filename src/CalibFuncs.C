#include "../include/CalibFuncs.h"
//______________________________________________________________________________
int GetOmegaP_err(perturbation_t pert,double T,double &err){
   // error on delta_t 
   int rc=0; 
   double SIG=0,DSIG=0;
   double delta_b=0,delta_b_err=0;
   rc = GetDiamagneticShielding(pert.sigma,pert.sigma_err,T,SIG,DSIG);
   rc = GetBulkMagneticSusceptibility(pert.eps,pert.eps_err,T,delta_b,delta_b_err);
   double err_sq   = DSIG*DSIG + delta_b_err*delta_b_err  
                   + pert.delta_s_err*pert.delta_s_err   + pert.delta_p_err*pert.delta_p_err
                   + pert.delta_rd_err*pert.delta_rd_err + pert.delta_d_err*pert.delta_d_err
                   + pert.delta_v_err*pert.delta_v_err;
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
   double delta_b=0,delta_b_err=0;
   rc = GetBulkMagneticSusceptibility(pert.eps,pert.eps_err,temp,delta_b,delta_b_err);
   double delta_tot  = GetDeltaTerm(sigma,delta_b,pert.delta_s,pert.delta_p,pert.delta_rd,pert.delta_d); 
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
double GetDeltaTerm(double sigma,double delta_b,double delta_s,double delta_p,
                    double delta_rd,double delta_d){
  // Compute the delta_tot term.   
  // - Input file carries the sign of the perturbation as measured (delta_s, delta_p, etc).  
  // - The factor of 1E-9 converts to 'absolute' scale, since the input is in ppb
  // - Must flip the sign on delta terms.  if these are negative, we need to INCREASE 
  //   the field.  Since we divide by (1-delta_tot), we need to ensure the field goes UP 
  //   because of this negative perturbation.  The situation is reversed for positive perturbations.
  double DELTA_S  = (-1.)*delta_s;   // material
  double DELTA_P  = (-1.)*delta_p;   // sample paramagnetic impurities 
  double DELTA_B  =  (1.)*delta_b;   // bulk magnetic susceptibility NOTE NO NEGATIVE SIGN HERE  
  double DELTA_RD = (-1.)*delta_rd;  // radiation damping  
  double DELTA_D  = (-1.)*delta_d;   // proton dipolar field   
  double delta_t  = 1E-9*(sigma + DELTA_B + DELTA_S + DELTA_P + DELTA_RD + DELTA_D);
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
int GetBulkMagneticSusceptibility(double eps,double deps,double chi,double dchi,double T,
                                  double &DELTA_B,double &DELTA_B_ERR){
   // central value; this is in SI formalism 
   // units are ppb 
   double CHI=0; 
   int rc      = GetMagneticSusceptibility(chi,T,CHI);  
   DELTA_B     = ( eps/(4.*TMath::Pi()) - 1./3. )*CHI;
   // calculate the uncertainty 
   double DEPS = ( CHI/(4.*TMath::Pi()) )*deps;
   double DCHI = ( eps/(4.*TMath::Pi()) - 1./3. )*dchi;
   DELTA_B_ERR = TMath::Sqrt( DEPS*DEPS + DCHI*DCHI);
   return 0; 
}
//______________________________________________________________________________
int GetMagneticSusceptibility(double chi,double dchi,double T,double &CHI,double &CHI_ERR){
   // compute magnetic susceptibility when accounting for temperature dependence 
   double a[3] = {1.38810E-4,-1.2685E-7,8.09E-10};
   double ARG  = 1. + a[0]*(T-20.) + a[1]*TMath::Power(T-20.,2.) + a[2]*TMath::Power(T-20.,3.); 
   CHI     = chi*ARG;
   CHI_ERR = dchi*ARG;
   return 0;
}
