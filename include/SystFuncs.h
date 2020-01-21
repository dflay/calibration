#ifndef SYSTEMATIC_FUNCTIONS_H 
#define SYSTEMATIC_FUNCTIONS_H

#include <cstdlib> 
#include <iostream>
#include <vector>

// a namespace for systematic study functions
namespace systFunc {  
   int RandomizeTimeValues(double delta,std::vector<double> &x); 
   int RandomizeFitValues(std::vector<double> &x,std::vector<double> dx); 
} // ::systFunc 

#endif  
