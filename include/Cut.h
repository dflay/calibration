#ifndef CUT_H
#define CUT_H

// a class to keep track of cuts in the data 

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <string> 

#include "gm2fieldImport.h"

class Cut { 

   private:
      int fVerbosity;  
      json fData;

   public: 
      Cut(std::string filepath="UNKNOWN.json",int verbosity=0);
      ~Cut();
     
      void SetVerbosity(int v) { fVerbosity = v; } 
 
      int Print();
      int LoadData(std::string filepath);

      bool CheckEvent_trly(); 
      bool CheckEvent_nmrAna(int run,int trace,int nzc,double ampl,double freq);  

}; 

#endif 
