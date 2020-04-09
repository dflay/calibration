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
      bool fUseOscCor,fRemoveFXPRDrift,fUseMisalignCor,fUseTempCor_pp;
      bool fSyst,fVaryDBTime_tr,fVarySwapTime_tr,fVaryShimFit,fVaryImpGradFit,fShimGradAltEnable;
      bool fTRLYFootprintStatus; 
      int fSystDirNum; 
      int fTrolleyProbe,fAxis,fFXPRListTag,fBlindUnits,fRunPeriod,fNumEventsToAvg,fNumEventsTimeWindow,fImpGradFitDim,fImpGradFitOrder;
      double fTempCor_tr,fTempCor_pp; 
      double fBlindScale,fTrolleyAngle,fDBZCurrent,fDBDeltaTime_tr,fSwapDeltaTime_tr,fTRLYFootprint,fTRLYFootprintErr; 
      std::string fType,fDevice,fRunDate,fFitFunc,fBlindLabel,fProdTag,fNMRANATag,fCutFile,fPPID,fOscCorType,fAnaDate;
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
      bool GetOscCorStatus()         const { return fUseOscCor;       } 
      bool GetFXPRDriftStatus()      const { return fRemoveFXPRDrift; }  
      bool GetMisalignCorStatus()    const { return fUseMisalignCor;  }  
      bool GetTempCorStatus()        const { return fUseTempCor;      } 
      bool GetTempCorStatus_pp()     const { return fUseTempCor_pp;   }
      bool GetSystStatus()           const { return fSyst;            }
      bool GetShimGradAltStatus()    const { return fShimGradAltEnable; } 
      bool GetTRLYFootprintStatus()  const { return fTRLYFootprintStatus; }  

      bool GetVaryTimeStatus(std::string dev="NONE",std::string type="NONE") const {
	 // determine if we will vary the times on event selection 
         if( dev.compare("tr")!=0 ) return false; // not the trolley
         if( dev.compare("tr")==0 && type.compare("db")==0 )   return fVaryDBTime_tr;  
         if( dev.compare("tr")==0 && type.compare("swap")==0 ) return fVarySwapTime_tr; 
	 return false;  // all other options failed, return false 
      }

      bool GetSystFitStatus(std::string type="NONE") const {
	 if( type.compare("shim")==0     ) return fVaryShimFit; 
	 if( type.compare("imp-grad")==0 ) return fVaryImpGradFit;
	 return false; 
      } 
     
      int GetTrolleyProbe()          const { return fTrolleyProbe;    } 
      int GetAxis()                  const { return fAxis;            }
      int GetFixedProbeListTag()     const { return fFXPRListTag;     }
      int GetBlindUnits()            const { return fBlindUnits;      }
      int GetRunPeriod()             const { return fRunPeriod;       } 
      int GetNumEventsToAvg()        const { return fNumEventsToAvg;  }  
      int GetNumEventsTimeWindow()   const { return fNumEventsTimeWindow; }  
      int GetImpGradFitDimension()   const { return fImpGradFitDim;   }  
      int GetImpGradFitOrder()       const { return fImpGradFitOrder; } 
      int GetSystDirNum()            const { return fSystDirNum;      }  
      int GetTRLYFootprint(double &x,double &e) { x = fTRLYFootprint; e = fTRLYFootprintErr; return 0; }  

      double GetBlindScale()         const { return fBlindScale;      }  
      double GetTrolleyAngle()       const { return fTrolleyAngle;    }  
      double GetDBZCurrent()         const { return fDBZCurrent;      } 
      double GetTempCor_pp()         const { return fTempCor_pp;      } 
      double GetTempCor_tr()         const { return fTempCor_tr;      } 

      double GetDeltaTime(std::string dev="NONE",std::string type="NONE"){
         // the step size over which to vary a time cut 
         if( dev.compare("tr")!=0 ) return 0; // not the trolley 
	 if( dev.compare("tr")==0 && type.compare("db")==0)   return fDBDeltaTime_tr; 
	 if( dev.compare("tr")==0 && type.compare("swap")==0) return fSwapDeltaTime_tr; 
	 return 0; // all other options failed, return 0
      } 

      std::string GetType()          const { return fType;            } 
      std::string GetDevice()        const { return fDevice;          } 
      std::string GetRunDate()       const { return fRunDate;         } 
      std::string GetAnaDate()       const { return fAnaDate;         } 
      std::string GetFitFunction()   const { return fFitFunc;         } 
      std::string GetBlindLabel()    const { return fBlindLabel;      } 
      std::string GetProductionTag() const { return fProdTag;         }
      std::string GetNMRANATag()     const { return fNMRANATag;       }
      std::string GetCutFile()       const { return fCutFile;         } 
      std::string GetPPID()          const { return fPPID;            } 
      std::string GetOscCorType()    const { return fOscCorType;      } 

      std::string GetValue(std::string key)                       const { return fParams[key]; }  
      std::string GetValue(std::string key,std::string subKey)    const { return fParams[key][subKey]; } 

      // identical to above, but more descriptive function name and templated 
      template <typename T> T GetValueFromKey(std::string key)                       const { return fParams[key]; }  
      template <typename T> T GetValueFromSubKey(std::string key,std::string subKey) const { return fParams[key][subKey]; } 

}; 

#endif 
