#include "../include/myRand.h"
//______________________________________________________________________________
namespace random_df {
   //______________________________________________________________________________
   double getRandomNumber(double range){
      TRandom3 *R = new TRandom3(0);
      double r = R->Rndm();
      double val = (-1.)*range + 2.*r*range;
      delete R;
      return val;
   }
   //______________________________________________________________________________
   double getRandomNumber(double min,double max){
      TRandom3 *R = new TRandom3(0);
      double r = R->Rndm();
      double val = min + r*(max-min);
      delete R;
      return val;
   }
   //______________________________________________________________________________
   double getRandomNumber_gaus(double mean,double sigma){
      TRandom3 *R = new TRandom3(0);
      double r = R->Gaus(mean,sigma);
      delete R;
      return r;
   }
}
