#ifndef CUT_H
#define CUT_H

// a class to keep track of cuts in the data 

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <string> 

#include "gm2fieldImport.h"
#include "gm2fieldFunc.h"

#include "plungingProbeAnaEvent.h"
#include "fixedProbeEvent.h"

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

      // standard cut functions 
      bool CheckEvent_trly(); 
      bool CheckEvent_nmrAna(int run,int trace,int nzc,double ampl,double freq);  

      // specialized function 
      int FilterPPData(int runPeriod,int probe,std::string type,std::string axis,
                       std::vector<plungingProbeAnaEvent_t> in,std::vector<plungingProbeAnaEvent_t> &out,
                       std::string inpath);

      int FilterFXPRData(int runPeriod,int probe,
	                 std::vector<averageFixedProbeEvent_t> in,std::vector<averageFixedProbeEvent_t> &out,
	                 std::string inpath); 

}; 

#endif 
