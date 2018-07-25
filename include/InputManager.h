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
      bool fIsSimple,fIsFullAnalysis,fIsBlind,fUseP2PFit,fIsFinalLocation,fUseAxis;
      int fTrolleyProbe,fAxis,fFXPRListTag;  
      std::string fType,fDevice,fAnaDate,fFitFunc;
      std::vector<int> fRunList; 
      std::vector<std::string> fRunLabel; 

   public: 
      InputManager();
      ~InputManager();

      void UseAxis(bool q=true) { fUseAxis = q; } 

      int Init(); 
      int ClearVectors(); 
      int Print();  
      int Load(std::string inpath);

      int GetRunList(std::vector<int> &v);        
      int GetRunLabels(std::vector<std::string> &v);       

      bool IsFullAnalysis()         const { return fIsFullAnalysis;  } 
      bool IsBlind()                const { return fIsBlind;         } 
      bool IsFinalLocation()        const { return fIsFinalLocation; } 
      bool UseP2PFit()              const { return fUseP2PFit;       } 
     
      int GetTrolleyProbe()         const { return fTrolleyProbe;    } 
      int GetAxis()                 const { return fAxis;            }
      int GetFixedProbeListTag()    const { return fFXPRListTag;     } 

      std::string GetType()         const { return fType;            } 
      std::string GetDevice()       const { return fDevice;          } 
      std::string GetAnalysisDate() const { return fAnaDate;         } 
      std::string GetFitFunction()  const { return fFitFunc;         } 

}; 

#endif 
