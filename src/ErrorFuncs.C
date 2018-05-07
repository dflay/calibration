#include "../include/ErrorFuncs.h"
//______________________________________________________________________________
int GetDriftError(std::vector<deltab_t> pp,std::vector<deltab_t> trly,
                  std::vector<grad_meas_t> gradImposed,std::vector<grad_meas_t> gradShim,
                  std::vector<double> &errFXPR,std::vector<double> &errTRLY){

   // compute the error due to FIELD DRIFT 
   // index on all vectors is x,y,z
   double dq=0,Delta=0;
   double delta_dq_fxpr=0,delta_dq_trly=0;
   double arg1_fxpr=0,arg2_fxpr=0,argFXPR=0;
   double arg1_trly=0,arg2_trly=0,argTRLY=0;
   const int N = pp.size();
   for(int i=0;i<N;i++){
      // some generic terms we need 
      Delta         = pp[i].dB - trly[i].dB;                   // Delta(pp) - Delta(trly) 
      dq            = Delta/gradImposed[i].grad;               // error in alignment (mm) 
      delta_dq_fxpr = GetDeltaDQ("fxpr",pp[i],trly[i],gradImposed[i]);   
      delta_dq_trly = GetDeltaDQ("trly",pp[i],trly[i],gradImposed[i]);   
      // now compute terms 
      // fxpr  
      arg1_fxpr = gradShim[i].drift_fxpr_err*dq;           // error in shimmed-field gradient  
      arg2_fxpr = gradShim[i].grad*delta_dq_fxpr;
      argFXPR   = TMath::Sqrt(arg1_fxpr*arg1_fxpr + arg2_fxpr*arg2_fxpr); 
      // trly  
      arg1_trly = gradShim[i].drift_trly_err*dq;           // error in shimmed-field gradient  
      arg2_trly = gradShim[i].grad*delta_dq_trly;
      argTRLY   = TMath::Sqrt(arg1_trly*arg1_trly + arg2_trly*arg2_trly);
      // fill vectors 
      errFXPR.push_back(argFXPR); 
      errTRLY.push_back(argTRLY); 
   }
   return 0;
}
//______________________________________________________________________________
double GetDeltaDQ(std::string type,deltab_t pp,deltab_t trly,
                  grad_meas_t gradImposed){
   // get error in coordinate 
   double T1_num=0;
   double T2_num=0;
   if( type.compare("fxpr")==0 ){
      // use FXPR drift error 
      T1_num = TMath::Sqrt( TMath::Power(pp.drift_fxpr,2.) + TMath::Power(trly.drift_fxpr,2.) ); 
      T2_num = (-1.)*(pp.dB-trly.dB)*(gradImposed.drift_fxpr); 
   }else{
      // use TRLY drift err 
      T1_num = TMath::Sqrt( TMath::Power(pp.drift_trly,2.) + TMath::Power(trly.drift_trly,2.) ); 
      T2_num = (-1.)*(pp.dB-trly.dB)*(gradImposed.drift_fxpr); 
   }
   double T1_den = gradImposed.grad; 
   double T2_den = TMath::Power(gradImposed.grad,2.);  
   double T1     = T1_num/T1_den;
   double T2     = T2_num/T2_den;
   double res    = TMath::Sqrt(T1*T1 + T2*T2);
   return res;
}
