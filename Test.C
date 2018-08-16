// test passing a function as an argument  

#include <cstdlib> 
#include <iostream>
#include <vector>
#include <string> 

double testFunc( double (*f)(double x),double v); 
double myFunc(double x);

int Test(){

   double v = testFunc(&myFunc,5); 
   std::cout << v << std::endl; 

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
