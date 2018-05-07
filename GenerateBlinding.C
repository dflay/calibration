// generate a random number for blinding 

#include <cstdlib> 
#include <iostream>
#include <fstream>

#include "TRandom3.h"

#include "./src/CustomExport.C"

int GenerateBlinding(){

   double range    = 1.*61.79;     // 1 ppm
   double rangeErr = 20.*0.06179;  // 20 ppb 

   TRandom3 *myRand = new TRandom3(0); 

   double r = myRand->Rndm();
   double val = (-1.)*range + 2.*range*r;  
   r = myRand->Rndm();
   double val2 = (-1.)*range + 2.*range*r;  
   
   r = myRand->Rndm();
   double x_err = rangeErr*r;  
   r = myRand->Rndm();
   double y_err = rangeErr*r;  
   r = myRand->Rndm();
   double z_err = rangeErr*r;  
   r = myRand->Rndm();
   double f1_err = rangeErr*r;  
   r = myRand->Rndm();
   double f2_err = rangeErr*r;  

   double V[7] = {val,val2,x_err,y_err,z_err,f1_err,f2_err}; 

   char outpath[200];
   sprintf(outpath,"./misc/blind/values.csv"); 
   PrintToFile(outpath,V); 

   return 0;
}
