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
      json fParams; 
      bool fIsSimple,fIsFullAnalysis,fIsBlind,fUseP2PFit,fIsFinalLocation;
      bool fUseAxis,fIsFreeProton,fLoadSwapTime,fLoadSCCTime,fUseTimeWeight,fUseTempCor;
      bool fUseOscCor; 
      int fTrolleyProbe,fAxis,fFXPRListTag,fBlindUnits,fRunPeriod,fNumEventsToAvg,fNumEventsTimeWindow; 
      double fBlindScale,fTrolleyAngle,fDBZCurrent; 
      std::string fType,fDevice,fRunDate,fFitFunc,fBlindLabel,fProdTag,fNMRANATag,fCutFile;
      std::vector<int> fRunList,fFXPRList; 
      std::vector<std::string> fRunLabel; 

      int Parse();
      int LoadFXPRList(); 

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
      int GetFXPRList(std::vector<int> &v);       

      bool DoesKeyExist(std::string keyName); 

      bool IsFullAnalysis()          const { return fIsFullAnalysis;  } 
      bool IsBlind()                 const { return fIsBlind;         } 
      bool IsFinalLocation()         const { return fIsFinalLocation; } 
      bool UseP2PFit()               const { return fUseP2PFit;       } 
      bool GetFreeProtonStatus()     const { return fIsFreeProton;    }
      bool GetSwapTimeStatus()       const { return fLoadSwapTime;    }  
      bool GetSCCTimeStatus()        const { return fLoadSCCTime;     }  
      bool GetTimeWeightStatus()     const { return fUseTimeWeight;   } 
      bool GetTempCorStatus()        const { return fUseTempCor;      } 
      bool GetOscCorStatus()         const { return fUseOscCor;       }  
     
      int GetTrolleyProbe()          const { return fTrolleyProbe;    } 
      int GetAxis()                  const { return fAxis;            }
      int GetFixedProbeListTag()     const { return fFXPRListTag;     }
      int GetBlindUnits()            const { return fBlindUnits;      }
      int GetRunPeriod()             const { return fRunPeriod;       } 
      int GetNumEventsToAvg()        const { return fNumEventsToAvg;  }  
      int GetNumEventsTimeWindow()   const { return fNumEventsTimeWindow; }  

      double GetBlindScale()         const { return fBlindScale;      }  
      double GetTrolleyAngle()       const { return fTrolleyAngle;    }  
      double GetDBZCurrent()         const { return fDBZCurrent;      }  

      std::string GetType()          const { return fType;            } 
      std::string GetDevice()        const { return fDevice;          } 
      std::string GetAnalysisDate()  const { return fRunDate;         } 
      std::string GetFitFunction()   const { return fFitFunc;         } 
      std::string GetBlindLabel()    const { return fBlindLabel;      } 
      std::string GetProductionTag() const { return fProdTag;         }
      std::string GetNMRANATag()     const { return fNMRANATag;       }
      std::string GetCutFile()       const { return fCutFile;         } 

      std::string GetValue(std::string key)                       const { return fParams[key]; }  
      std::string GetValue(std::string key,std::string subKey)    const { return fParams[key][subKey]; } 

      // identical to above, but more descriptive function name and templated 
      template <typename T> T GetValueFromKey(std::string key)                       const { return fParams[key]; }  
      template <typename T> T GetValueFromSubKey(std::string key,std::string subKey) const { return fParams[key][subKey]; } 

}; 

#endif 
