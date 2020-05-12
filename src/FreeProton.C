#include "../include/FreeProton.h"
//______________________________________________________________________________
FreeProton::FreeProton(const char *inpath){
   Clear();
   int rc=0;
   std::string thePath = inpath;
   if(thePath.compare("NONE")!=0){
      rc = LoadData(inpath);
   }
   fT0_sigma = 25.0;
   fT0_chi   = 20.0;
}
//______________________________________________________________________________
FreeProton::~FreeProton(){

}
//______________________________________________________________________________
void FreeProton::Clear(){
   fProbeID      = "NONE"; 
   fsigma        = 0;
   feps          = 0;
   fchi          = 0;
   fdelta_s      = 0;
   fdelta_p      = 0;
   fdelta_b      = 0;
   fdelta_rd     = 0;
   fdelta_d      = 0;
   fdelta_v      = 0;
   fsigma_err    = 0;
   feps_err      = 0;
   fchi_err      = 0;
   fdelta_s_err  = 0;
   fdelta_p_err  = 0;
   fdelta_b_err  = 0;
   fdelta_rd_err = 0;
   fdelta_d_err  = 0;
   fdelta_v_err  = 0;
}
//______________________________________________________________________________
int FreeProton::LoadData(const char *inpath){
   json obj;
   std::string inpath_str = inpath;
   int rc = gm2fieldUtil::Import::ImportJSON(inpath_str,obj);
   if(rc!=0){
      std::cout << "[FreeProton::LoadData]: Cannot open the file: " << inpath_str << std::endl;
      return 1;
   }else{
      fProbeID       = obj["probe-id"]; 
      // convert all to absolute values! 
      fsigma         = obj["sigma"];
      fsigma_err     = obj["sigma-err"];
      fchi           = obj["chi"];
      fchi_err       = obj["chi-err"];
      feps           = obj["eps"];  
      feps_err       = obj["eps-err"];
      fdelta_s       = obj["delta-s"];
      fdelta_s_err   = obj["delta-s-err"];
      fdelta_p       = obj["delta-p"];
      fdelta_p_err   = obj["delta-p-err"];
      fdelta_rd      = obj["delta-rd"];
      fdelta_rd_err  = obj["delta-rd-err"];
      fdelta_d       = obj["delta-d"];
      fdelta_d_err   = obj["delta-d-err"];
      fdelta_v       = obj["delta-v"];
      fdelta_v_err   = obj["delta-v-err"];
      // convert to absolute numbers 
      fsigma         *= 1E-9; 
      fsigma_err     *= 1E-9; 
      fchi           *= 1E-9; 
      fchi_err       *= 1E-9; 
      feps           *= 1.; // this is already absolute 
      feps_err       *= 1.; // this is already absolute 
      fdelta_s       *= 1E-9; 
      fdelta_s_err   *= 1E-9; 
      fdelta_p       *= 1E-9; 
      fdelta_p_err   *= 1E-9; 
      fdelta_rd      *= 1E-9; 
      fdelta_rd_err  *= 1E-9; 
      fdelta_d       *= 1E-9; 
      fdelta_d_err   *= 1E-9; 
      fdelta_v       *= 1E-9; 
      fdelta_v_err   *= 1E-9; 
   }

   return 0;
}
//______________________________________________________________________________
double FreeProton::GetOmegaP_free(double freq,double T){
   double delta_t   = GetDelta_t(T); 
   double freq_free = freq/(1.-delta_t); 
   return freq_free;
}
//______________________________________________________________________________
double FreeProton::GetDelta_t_err(double T){
   // error on delta_t 
   int rc=0; 
   double SIG=0,DSIG=0;
   double delta_b=0,delta_b_err=0;
   // temperature-dependent terms 
   CalculateDiamagneticShielding(T,SIG,DSIG);
   CalculateBulkMagneticSusceptibility(T,delta_b,delta_b_err);
   // full calculation 
   double err_sq = DSIG*DSIG + delta_b_err*delta_b_err   
                 + fdelta_s_err*fdelta_s_err   + fdelta_p_err*fdelta_p_err
                 + fdelta_rd_err*fdelta_rd_err + fdelta_d_err*fdelta_d_err
                 + fdelta_v_err*fdelta_v_err;
   double err = TMath::Sqrt(err_sq);  
   return err;
}
//______________________________________________________________________________
double FreeProton::GetDelta_t(double T){
  // Compute the delta_tot term.   
  // - Input file carries the sign of the perturbation as measured (delta_s, delta_p, etc).  
  // - Must flip the sign on delta terms.  if these are negative, we need to INCREASE 
  //   the field.  Since we divide by (1-delta_tot), we need to ensure the field goes UP 
  //   because of this negative perturbation.  The situation is reversed for positive perturbations.

  // first evaluate temperature-dependent terms 
  double sigma=0,sigma_err=0; 
  CalculateDiamagneticShielding(T,sigma,sigma_err);

  double delta_b=0,delta_b_err=0; 
  CalculateBulkMagneticSusceptibility(T,delta_b,delta_b_err); 

  // now evaluate 
  double DELTA_S  = (-1.)*fdelta_s;   // material
  double DELTA_P  = (-1.)*fdelta_p;   // sample paramagnetic impurities 
  double DELTA_B  =  (1.)*delta_b;    // bulk magnetic susceptibility NOTE NO NEGATIVE SIGN HERE  
  double DELTA_RD = (-1.)*fdelta_rd;  // radiation damping  
  double DELTA_D  = (-1.)*fdelta_d;   // proton dipolar field   
  double delta_t  = sigma + DELTA_B + DELTA_S + DELTA_P + DELTA_RD + DELTA_D;
  return delta_t;
}
//______________________________________________________________________________
double FreeProton::GetChi(double T){
   double CHI=0,CHI_ERR=0;
   CalculateMagneticSusceptibility(T,CHI,CHI_ERR);
   return CHI;
}
//______________________________________________________________________________
double FreeProton::GetChi_err(double T){
   double CHI=0,CHI_ERR=0;
   CalculateMagneticSusceptibility(T,CHI,CHI_ERR);
   return CHI_ERR;
}
//______________________________________________________________________________
double FreeProton::GetSigma(double T){
  double sigma=0,sigma_err=0; 
  CalculateDiamagneticShielding(T,sigma,sigma_err);
  return sigma;
}
//______________________________________________________________________________
double FreeProton::GetSigma_err(double T){
  double sigma=0,sigma_err=0; 
  CalculateDiamagneticShielding(T,sigma,sigma_err);
  return sigma_err;
}
//______________________________________________________________________________
double FreeProton::GetDelta_b(double T){
  double delta_b=0,delta_b_err=0; 
  CalculateBulkMagneticSusceptibility(T,delta_b,delta_b_err); 
  return delta_b;
}
//______________________________________________________________________________
double FreeProton::GetDelta_b_err(double T){
  double delta_b=0,delta_b_err=0; 
  CalculateBulkMagneticSusceptibility(T,delta_b,delta_b_err); 
  return delta_b_err;
}
//______________________________________________________________________________
void FreeProton::CalculateBulkMagneticSusceptibility(double T,double &delta_b,double &delta_b_err){
   // central value; this is in SI formalism 
   // units are ppb 
   double CHI=0,CHI_ERR=0;
   CalculateMagneticSusceptibility(T,CHI,CHI_ERR);
   delta_b     = (feps - 1./3.)*CHI;
   // calculate the uncertainty 
   double DEPS = feps_err/(feps - 1./3.);
   double DCHI = CHI_ERR/CHI;
   delta_b_err = TMath::Abs(delta_b)*TMath::Sqrt(DEPS*DEPS + DCHI*DCHI);
}
//______________________________________________________________________________
void FreeProton::CalculateMagneticSusceptibility(double T,double &CHI,double &CHI_ERR){
   // compute magnetic susceptibility when accounting for temperature dependence
   double DT   = T - fT0_chi;  // NOTE: SIGN IS DIFFERENT FROM DIAMAGNETIC SHIELDING!
   // where does this come from?
   double a[3] = {1.38810E-4,-1.2685E-7,8.09E-10};
   // from J Chem Phys 72, 4434 
   // double a[3]  = {0.997218,1.0213,-14.12E-5};
   // double da[3] = {0.000005,0.0007,0.03E-5};
   double ARG  = 1. + a[0]*DT + a[1]*TMath::Power(DT,2.) + a[2]*TMath::Power(DT,3.);
   CHI     = fchi*ARG;
   CHI_ERR = fchi_err*ARG;
}
//______________________________________________________________________________
void FreeProton::CalculateDiamagneticShielding(double T,double &SIG,double &ERR){
   // compute diamagnetic shielding with temperature dependence
   // sigma_T0 = sigma @ T = 25 deg C
   double DT         = fT0_sigma - T;  
   double dsigdT     = -10.36E-9; 
   double dsigdT_err =   0.30E-9;
   SIG = fsigma + dsigdT*DT;
   ERR = TMath::Sqrt(fsigma_err*fsigma_err + DT*DT*dsigdT_err*dsigdT_err);  
}
//______________________________________________________________________________
void FreeProton::Print(std::string units){
 
   double sf=1;
   std::string unitStr = ""; 
   if( units.compare("ppb")==0 ){
      sf      = 1E-9;
      unitStr = "ppb"; 
   }

   std::cout << "-------------------------------------------------" << std::endl; 
   std::cout << "Probe: " << fProbeID  << std::endl;
   std::cout << Form("delta_s             = %.1lf ± %.1lf %s",fdelta_s /sf,fdelta_s_err /sf,unitStr.c_str() ) << std::endl;
   std::cout << Form("delta_p             = %.1lf ± %.1lf %s",fdelta_p /sf,fdelta_p_err /sf,unitStr.c_str() ) << std::endl;
   std::cout << Form("delta_rd            = %.1lf ± %.1lf %s",fdelta_rd/sf,fdelta_rd_err/sf,unitStr.c_str() ) << std::endl;
   std::cout << Form("delta_d             = %.1lf ± %.1lf %s",fdelta_d /sf,fdelta_d_err /sf,unitStr.c_str() ) << std::endl;
   std::cout << Form("delta_v             = %.1lf ± %.1lf %s",fdelta_v /sf,fdelta_v_err /sf,unitStr.c_str() ) << std::endl;
   std::cout << Form("eps                 = %.8lf ± %.8lf"   ,feps,feps_err)                                  << std::endl;  
   std::cout << Form("sigma(T = 25 deg C) = %.1lf ± %.1lf %s",fsigma/sf,fsigma_err/sf,unitStr.c_str() )       << std::endl;  
   std::cout << Form("chi(T = 20 deg C)   = %.1lf ± %.1lf %s",fchi/sf,fchi_err/sf,unitStr.c_str() )           << std::endl;  
   std::cout << "-------------------------------------------------" << std::endl;
}
