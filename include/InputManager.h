#ifndef INPUT_MANAGER_H
#define INPUT_MANAGER_H

// a manager to keep track of all our input parameters 

#include <cstdlib> 
#include <vector> 
#include <iostream> 
#include <string>

#include "gm2fieldImport.h"

class InputManager{

   private:
      bool fIsFullAnalysis,fIsBlind,fUseP2PFit,fIsFinalLocation;
      int fTrolleyProbe,fAxis;  
      std::string fAnaDate,fFitFunc;
      std::vector<int> fRunList; 
      std::vector<std::string> fRunLabel; 

   public: 
      InputManager();
      ~InputManager();

      bool IsFullAnalysis()  const { return fIsFullAnalysis;  } 
      bool IsBlind()         const { return fIsBlind;         } 
      bool IsFinalLocation() const { return fIsFinalLocation; } 
      bool UseP2PFit()       const { return fUseP2PFit;       } 
     
      int GetTrolleyProbe()  const { return fTrolleyProbe;    } 
      int GetAxis()          const { return fAxis;            }
 
      int Init(); 
      int ClearVectors(); 
      int Print();  
      int Load(std::string inpath);
      int GetRunList(std::vector<int> &v);        
      int GetRunLabels(std::vector<std::string> &v);       

      std::string GetAnalysisDate() const { return fAnaDate; } 
      std::string GetFitFunction()  const { return fFitFunc; } 

}; 

#endif 
