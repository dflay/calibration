#ifndef SYSTEMATIC_FUNCTIONS_H 
#define SYSTEMATIC_FUNCTIONS_H

#include <cstdlib> 
#include <iostream>
#include <vector>

// a namespace for systematic study functions
namespace systFunc {  
   int RandomizeTimeValues(double delta,std::vector<double> &x); 
   int RandomizeFitValue(double &x,double dx); 
   int RandomizeFitValues(std::vector<double> &x,std::vector<double> dx);
   int RandomizeFitValues(int NP,double *x,double *dx); 
} // ::systFunc 

#endif  
