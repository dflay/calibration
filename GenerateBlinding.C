// generate random numbers for blinding 

#include <cstdlib> 
#include <iostream>
#include <fstream>

#include "Blinder.h"

int GenerateBlinding(){

   int units    = gm2fieldUtil::Constants::ppb; 
   double range = 100;

   gm2fieldUtil::Blinder *myBlind = new gm2fieldUtil::Blinder("flay",range,units); 
   myBlind->UpdateBlinding();  

   delete myBlind; 

   return 0;
}
