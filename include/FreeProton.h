#ifndef FREE_PROTON_H
#define FREE_PROTON_H

// a class for calculating free proton corrections 

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <string> 

#include "TMath.h"

#include "gm2fieldImport.h"

class FreeProton { 

   private:

      std::string fProbeID; 
 
      double fsigma,feps,fchi,fdelta_s,fdelta_p,fdelta_b,fdelta_rd,fdelta_d,fdelta_v;
      double fsigma_err,feps_err,fchi_err,fdelta_s_err,fdelta_p_err,fdelta_b_err,fdelta_rd_err,fdelta_d_err,fdelta_v_err;
      double fsigma_T,fchi_T; 
      double fsigma_T_err,fchi_T_err; 

      void CalculateDiamagneticShielding(double T,double &SIG,double &ERR);
      void CalculateBulkMagneticSusceptibility(double T,double &delta_b,double &delta_b_err); 
      void CalculateMagneticSusceptibility(double T,double &CHI,double &CHI_ERR); 

   public:

      FreeProton(std::string probeName="NONE",int runPeriod=1); 
      ~FreeProton(); 

      void Clear(); 

      int LoadData(std::string probeName,int runPeriod); 

      double GetOmegaP_free(double freq,double T); 

      double GetDelta_t(double T);
      double GetDelta_t_err(double T);

      // These don't change as a function of T
      double GetDelta_s()          const { return fdelta_s;  } 
      double GetDelta_p()          const { return fdelta_p;  } 
      double GetDelta_rd()         const { return fdelta_rd; } 
      double GetDelta_d()          const { return fdelta_d;  } 
      double GetDelta_v()          const { return fdelta_v;  }

      double GetDelta_s_err()      const { return fdelta_s_err;  } 
      double GetDelta_p_err()      const { return fdelta_p_err;  } 
      double GetDelta_rd_err()     const { return fdelta_rd_err; } 
      double GetDelta_d_err()      const { return fdelta_d_err;  } 
      double GetDelta_v_err()      const { return fdelta_v_err;  }

      // functions of T
      double GetDelta_b(double T); 
      double GetSigma(double T);
      double GetChi(double T);  

      double GetDelta_b_err(double T); 
      double GetSigma_err(double T);
      double GetChi_err(double T);  

      std::string GetProbeID()     const { return fProbeID; }   

}; 

#endif 
