#include "../include/InputManager.h"
//______________________________________________________________________________
InputManager::InputManager(){
   Init(); 
}
//______________________________________________________________________________
InputManager::~InputManager(){
   ClearVectors(); 
}
//______________________________________________________________________________
int InputManager::Init(){
   fIsSimple        = false; 
   fIsFullAnalysis  = false; 
   fIsBlind         = false; 
   fUseP2PFit       = false; 
   fIsFinalLocation = false;
   fIsFreeProton    = false; 
   fUseAxis         = false;
   fLoadSwapTime    = false; 
   fLoadSCCTime     = false;
   fUseTimeWeight   = false; 
   fUseTempCor      = false;
   fUseTempCor_pp   = false;
   fUseOscCor       = false; 
   fRemoveFXPRDrift = false;
   fUseMisalignCor  = false; 
   fVaryDBTime_tr   = false;  
   fVarySwapTime_tr = false; 
   fVaryShimFit     = false; 
   fVaryImpGradFit  = false; 
   fSyst            = false;  
   fTempCor_pp      = 0;
   fNumEventsToAvg  = 0;  
   fNumEventsTimeWindow = 0; 
   fTrolleyAngle     = 0; 
   fDBZCurrent       = 0;
   fDBDeltaTime_tr   = 0; 
   fSwapDeltaTime_tr = 0; 
   fTrolleyProbe     = -1; 
   fAxis             = -1;
   fFXPRListTag      = -1;
   fRunPeriod        = -1;
   fBlindUnits       = -1;  
   fBlindScale       = -1;
   fImpGradFitDim    = -1;  
   fImpGradFitOrder  = -1;  
   fType             = "NONE";
   fDevice           = "NONE";  
   fRunDate          = "NONE";  
   fFitFunc          = "NONE";
   fProdTag          = "NONE"; 
   fNMRANATag        = "NONE";
   fPPID             = "NONE";  
   ClearVectors(); 
   return 0;
}
//______________________________________________________________________________
int InputManager::ClearVectors(){
   fRunList.clear();
   fFXPRList.clear();
   fRunLabel.clear();
   return 0;
}
//______________________________________________________________________________
bool InputManager::DoesKeyExist(std::string keyName){
   auto it_key = fParams.find(keyName);  // this is an iterator 
   if (it_key!=fParams.end() ){ 
      return true;  // not at the end of fParams -- found the key 
   }else{
      return false;  // at the end of fParams -- didn't find the key 
   }
}
//______________________________________________________________________________
int InputManager::GetRunList(std::vector<int> &v){
   const int N = fRunList.size();
   if(N==0){
      std::cout << "[InputManager::GetRunList]: No runs!" << std::endl;
      return 1;
   }
   for(int i=0;i<N;i++){
      v.push_back(fRunList[i]); 
   }
   return 0;
}
//______________________________________________________________________________
int InputManager::LoadFXPRList(){

   int N=0,ifp=0;
   char inpath[200];
   sprintf(inpath,"./input/probe-lists/fxpr-list_set-%d.csv",fFXPRListTag); 
   std::ifstream infile; 
   infile.open(inpath); 

   if( infile.fail() ){
      std::cout << "[InputManager::LoadFXPRList]: Cannot open the file: " << inpath << std::endl;
      return 1;
   }else{
      while( !infile.eof() ){
	 infile >> ifp;
	 fFXPRList.push_back(ifp);  
      }
      N = fFXPRList.size();
      if(N>1) fFXPRList.pop_back();
      infile.close();
   } 
   return 0;
}
//______________________________________________________________________________
int InputManager::GetFXPRList(std::vector<int> &v){
   const int N = fFXPRList.size();
   if(N==0){
      std::cout << "[InputManager::GetFXPRList]: No runs!" << std::endl;
      return 1;
   }
   for(int i=0;i<N;i++){
      v.push_back(fFXPRList[i]); 
   }
   return 0;
}
//______________________________________________________________________________
int InputManager::GetRunLabels(std::vector<std::string> &v){
   const int N = fRunList.size();
   if(N==0){
      std::cout << "[InputManager::GetRunLabels]: No runs!" << std::endl;
      return 1;
   }
   for(int i=0;i<N;i++){
      v.push_back(fRunLabel[i]); 
   }
   return 0;
}
//______________________________________________________________________________
int InputManager::Load(std::string inpath){
   int NRUNS=0; 
   int rc = gm2fieldUtil::Import::ImportJSON(inpath,fParams);
 
   bool typeKeyExist = DoesKeyExist("type"); 
   if(typeKeyExist){
      fType = fParams["type"]; 
      Parse();
   }else{
      std::cout << "[InputManager::Load]: Cannot find the key 'type'.  Exiting" << std::endl;
      return 1;
   }

   return rc; 
}
//______________________________________________________________________________
int InputManager::Parse(){
   int NRUNS=0; 

   // parse all the parameters we want
   bool axisStatus   = DoesKeyExist("axis"); 
   bool anaStatus    = DoesKeyExist("full-ana"); 
   bool locStatus    = DoesKeyExist("final-loc"); 
   bool p2pStatus    = DoesKeyExist("p2p-fit"); 
   bool runStatus    = DoesKeyExist("nruns"); 
   bool dateStatus   = DoesKeyExist("date"); 
   bool blindStatus  = DoesKeyExist("blinding"); 
   bool trlyStatus   = DoesKeyExist("trly-probe");
   bool fitStatus    = DoesKeyExist("fit"); 
   bool swapStatus   = DoesKeyExist("load-trly-swap-times");  
   bool sccStatus    = DoesKeyExist("load-trly-scc-times");  
   bool timeStatus   = DoesKeyExist("use-aba-time-weight");
   bool tempStatus   = DoesKeyExist("use-trly-temp-cor");
   bool oscStatus    = DoesKeyExist("osc-cor");
   bool runpStatus   = DoesKeyExist("run-period"); 
   bool prodStatus   = DoesKeyExist("prod-tag");  
   bool nmrAnaStatus = DoesKeyExist("nmr-ana-tag"); 
   bool cutStatus    = DoesKeyExist("cut-file"); 
   bool nevStatus    = DoesKeyExist("num-events-to-avg");  
   bool twStatus     = DoesKeyExist("num-events-time-window"); 
   bool angleStatus  = DoesKeyExist("trly_azi-angle"); 
   bool curStatus    = DoesKeyExist("dBz-current");
   bool mcorStatus   = DoesKeyExist("use-misalign-cor");  
   bool ppStatus     = DoesKeyExist("pp");
   bool impStatus    = DoesKeyExist("imp-grad"); 
   bool systStatus   = DoesKeyExist("syst");  

   // parameters common to all 
   std::string unitStr="";
   if(dateStatus)   fRunDate      = fParams["date"];
   if(prodStatus)   fProdTag      = fParams["prod-tag"]; 
   if(nmrAnaStatus) fNMRANATag    = fParams["nmr-ana-tag"]; 
   if(blindStatus){
      fIsBlind    = (bool)( (int)fParams["blinding"]["enable"] );
      fBlindLabel = fParams["blinding"]["label"];
      fBlindScale = (double)fParams["blinding"]["scale"]; 
      unitStr     = fParams["blinding"]["units"];
      if( unitStr.compare("Hz")==0  ) fBlindUnits = gm2fieldUtil::Constants::Hz;  
      if( unitStr.compare("ppm")==0 ) fBlindUnits = gm2fieldUtil::Constants::ppm;  
      if( unitStr.compare("ppb")==0 ) fBlindUnits = gm2fieldUtil::Constants::ppb;  
   } 

   if(runStatus){ 
      NRUNS = (int)fParams["nruns"];
      for(int i=0;i<NRUNS;i++){
	 fRunList.push_back( (int)fParams["run-list"][i] ); 
	 fRunLabel.push_back( fParams["run-label"][i] ); 
      }
   }

   if(mcorStatus)  fUseMisalignCor  = (bool)( (int)fParams["use-misalign-cor"] );  
  
   if(systStatus){
      fSyst             = (bool)( (int)fParams["syst"]["enable"]         ); 
      fVaryDBTime_tr    = (bool)( (int)fParams["syst"]["vary-db-time"]   ); 
      fVarySwapTime_tr  = (bool)( (int)fParams["syst"]["vary-swap-time"] ); 
      fVaryShimFit      = (bool)( (int)fParams["syst"]["vary-shim-fit"]  ); 
      fVaryImpGradFit   = (bool)( (int)fParams["syst"]["vary-igrad-fit"] );
      fDBDeltaTime_tr   = (double)fParams["syst"]["tr-db-delta"];  
      fSwapDeltaTime_tr = (double)fParams["syst"]["tr-swap-delta"];  
   }
 
   // calibration: production 
   if( fType.compare("calib-prod")==0 ){
      if(trlyStatus) fTrolleyProbe   = (int)fParams["trly-probe"]; 
      if(runpStatus) fRunPeriod      = fParams["run-period"];  
      if(axisStatus) fAxis           = (int)fParams["axis"];
      if(fitStatus)  fFitFunc        = fParams["fit"];
      if(swapStatus) fLoadSwapTime   = (bool)( (int)fParams["load-trly-swap-times"] ); 
      if(sccStatus)  fLoadSCCTime    = (bool)( (int)fParams["load-trly-scc-times"] ); 
      if(timeStatus) fUseTimeWeight  = (bool)( (int)fParams["use-aba-time-weight"] ); 
      if(tempStatus) fUseTempCor     = (bool)( (int)fParams["use-trly-temp-cor"] ); 
      if(oscStatus){
	 fUseOscCor       = (bool)( (int)fParams["osc-cor"]["enable"] );
         fOscCorType      = fParams["osc-cor"]["type"];  
	 fRemoveFXPRDrift = (bool)( (int)fParams["osc-cor"]["fxpr-remove-drift"] );
	 fFXPRListTag     = (int)fParams["osc-cor"]["fxpr-set"]; 
	 LoadFXPRList();
      } 
      if(cutStatus)  fCutFile        = fParams["cut-file"];
      if(nevStatus)  fNumEventsToAvg = (int)fParams["num-events-to-avg"];  
      if(twStatus )  fNumEventsTimeWindow = (int)fParams["num-events-time-window"]; 
      if(angleStatus) fTrolleyAngle  = (double)fParams["trly_azi-angle"]; 
      if(curStatus)   fDBZCurrent    = (double)fParams["dBz-current"];
      if(ppStatus){
	 fPPID           = fParams["pp"]["id"];  
	 fUseTempCor_pp  = (bool)( (int)fParams["pp"]["temp-cor-enable"] ); // 2 wire => 4 wire correction! 
	 fTempCor_pp     = (double)fParams["pp"]["temp-cor-value"];
         fIsFreeProton   = (bool)( (int)fParams["pp"]["free-proton-enable"] );  
      }
      if(impStatus){
	 fImpGradFitDim   = (int)fParams["imp-grad"]["fit-dimension"]; 
	 fImpGradFitOrder = (int)fParams["imp-grad"]["fit-order"]; 
      } 
   }

   // trolley Delta-B measurements
   if( fType.compare("trly-db")==0 ){
      if(trlyStatus) fTrolleyProbe  = (int)fParams["trly-probe"]; 
      if(runpStatus) fRunPeriod     = fParams["run-period"];  
      if(axisStatus) fAxis          = (int)fParams["axis"];
      if(sccStatus)  fLoadSCCTime   = (bool)( (int)fParams["load-trly-scc-times"] ); 
      if(timeStatus) fUseTimeWeight = (bool)( (int)fParams["use-aba-time-weight"] ); 
   }
 
   // simple input format (i.e., rapid swapping data from 6/1)  
   if( fType.compare("is-simple")==0 ){
      fIsSimple                      = true;
      if(trlyStatus) fTrolleyProbe   = (int)fParams["trly-probe"]; 
      if(axisStatus) fAxis           = (int)fParams["axis"]; 
      if(anaStatus) fIsFullAnalysis  = (bool)( (int)fParams["full-ana"]  );    
      if(fitStatus) fFitFunc         = fParams["fit"]; 
   }else if( fType.compare("is-complex")==0 ){
      fIsSimple                      = false;
      if(trlyStatus) fTrolleyProbe   = (int)fParams["trly-probe"]; 
      if(axisStatus) fAxis           = (int)fParams["axis"]; 
      if(anaStatus) fIsFullAnalysis  = (bool)( (int)fParams["full-ana"]  );    
      if(locStatus) fIsFinalLocation = (bool)( (int)fParams["final-loc"] );   
      if(p2pStatus) fUseP2PFit       = (bool)( (int)fParams["p2p-fit"]   ); 
      if(fitStatus) fFitFunc         = fParams["fit"]; 
   }

   return 0; 
}
//______________________________________________________________________________
int InputManager::Print(){
   char axis='N';
   if(fAxis==0) axis = 'x'; 
   if(fAxis==1) axis = 'y'; 
   if(fAxis==2) axis = 'z'; 
   std::cout << "------------------- Input Manager -------------------" << std::endl;
   std::cout << "Plunging Probe SN:  " << fPPID            << std::endl; 
   std::cout << "Trolley probe:      " << fTrolleyProbe    << std::endl;
   std::cout << "Is blinded:         " << fIsBlind         << std::endl;
   std::cout << "Blind label:        " << fBlindLabel      << std::endl;
   std::cout << "Production tag:     " << fProdTag         << std::endl; 
   std::cout << "NMR-ANA tag:        " << fNMRANATag       << std::endl; 
   std::cout << "Oscillation cor:    " << fUseOscCor       << std::endl;
   std::cout << "Misalignment cor:   " << fUseMisalignCor  << std::endl; 
   std::cout << "Trolley angle:      " << fTrolleyAngle    << std::endl; 
   std::cout << "dBz current:        " << fDBZCurrent      << std::endl; 
   if(fType.compare("calib-prod")==0){
         std::cout << "Run period:         " << fRunPeriod    << std::endl;
	 std::cout << "Axis:               " << fAxis << " (" << axis << ")" << std::endl;
         std::cout << "Load TRLY swap time " << fLoadSwapTime << std::endl;
   }else{
      if(fIsSimple){
	 std::cout << "FXPR set:           " << fFXPRListTag     << std::endl;
      }else if(!fIsSimple){
	 std::cout << "Run date:           " << fRunDate         << std::endl;
	 std::cout << "Is full analysis:   " << fIsFullAnalysis  << std::endl;
	 std::cout << "Is final location:  " << fIsFinalLocation << std::endl;
	 std::cout << "Use P2P fit:        " << fUseP2PFit       << std::endl;
	 std::cout << "Fit function:       " << fFitFunc         << std::endl;
	 std::cout << "Axis:               " << fAxis            << " (" << axis << ")" << std::endl;
      }
   }
   int N = fRunList.size(); 
   std::cout << "Number of runs:     " << N                << std::endl;
   for(int i=0;i<N;i++) std::cout << Form("run = %d, label = %s",fRunList[i],fRunLabel[i].c_str()) << std::endl;
   std::cout << "-----------------------------------------------------" << std::endl;
   return 0;
}
