#ifndef MYRANDOM_H
#define MYRANDOM_H 

// a namespace for generating random numbers

#include "TRandom3.h"

namespace random_df {
   double getRandomNumber(double range); 
   double getRandomNumber(double min,double max); 
   double getRandomNumber_gaus(double mean,double sigma); 
} 

#endif 
