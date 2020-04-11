// test passing a function as an argument  

#include <cstdlib> 
#include <iostream>
#include <vector>
#include <string> 

#include "./src/myRand.C"

double testFunc( double (*f)(double x),double v); 
double myFunc(double x);

int Test(){

   double v = testFunc(&myFunc,5); 
   std::cout << v << std::endl; 

   double mean = 61.79E+6; 
   double sig  = 1E+3; 

   double freq = random_df::getRandomNumber_gaus(mean,sig);

   std::cout << Form("%.3lf",freq) << std::endl;

   return 0;
}
//______________________________________________________________________________
double testFunc( double (*f)(double x),double v){
   double val = (*f)(v);
   return val; 
} 
//______________________________________________________________________________
double myFunc(double x){
   return x;
}
