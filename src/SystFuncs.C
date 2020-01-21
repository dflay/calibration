#include "../include/SystFuncs.h"
#include "myRand.C"
namespace systFunc { 
   //______________________________________________________________________________
   int RandomizeTimeValues(double delta,std::vector<double> &x){
      // randomize the value in x[i] on the range (x[i]-delta,x[i]) 
      // specific to randomized times for event selection. 
      // since we count events going *backwards* in time, we randomize 
      // on the range of values leading up to x[i].  
      double rval=0;
      const int N = x.size();
      for(int i=0;i<N;i++){
	 rval = random_df::getRandomNumber(x[i]-delta,x[i]);
	 x[i] = rval; 
      }
      return 0;
   }
   //______________________________________________________________________________
   int RandomizeFitValues(std::vector<double> &x,std::vector<double> dx){
      // randomize the value in x[i] on the range Gaus(x[i],dx[i]) 
      // specific to varying fit parameters within their 1-sigma uncertainties.  
      double rval=0;
      const int N = x.size();
      for(int i=0;i<N;i++){
	 rval = random_df::getRandomNumber_gaus(x[i],dx[i]);
	 x[i] = rval; 
      }
      return 0;
   }
   //______________________________________________________________________________
   int RandomizeFitValues(int NP,double *x,double *dx){
      // randomize the value in x[i] on the range Gaus(x[i],dx[i]) 
      // specific to varying fit parameters within their 1-sigma uncertainties.  
      double rval=0;
      for(int i=0;i<NP;i++){
	 rval = random_df::getRandomNumber_gaus(x[i],dx[i]);
	 x[i] = rval; 
      }
      return 0;
   }
}
