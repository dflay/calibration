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
   fTrolleyProbe    = -1; 
   fAxis            = -1;
   fFXPRListTag     = -1;
   fType            = "NONE";
   fDevice          = "NONE";  
   fAnaDate         = "NONE";  
   fFitFunc         = "NONE";
   ClearVectors(); 
   return 0;
}
//______________________________________________________________________________
int InputManager::ClearVectors(){
   fRunList.clear();
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
   bool axisStatus  = DoesKeyExist("axis"); 
   bool anaStatus   = DoesKeyExist("full-ana"); 
   bool locStatus   = DoesKeyExist("final-loc"); 
   bool p2pStatus   = DoesKeyExist("p2p-fit"); 
   bool protStatus  = DoesKeyExist("free-proton-cor");
   bool runStatus   = DoesKeyExist("nruns"); 
   bool dateStatus  = DoesKeyExist("date"); 
   bool blindStatus = DoesKeyExist("blinding"); 
   bool trlyStatus  = DoesKeyExist("trly-probe");
   bool fxprStatus  = DoesKeyExist("fxpr-set"); 
   bool fitStatus   = DoesKeyExist("fit"); 
   bool swapStatus  = DoesKeyExist("load-trly-swap-times");  
   bool sccStatus   = DoesKeyExist("load-trly-scc-times");  
   bool timeStatus  = DoesKeyExist("use-aba-time-weight");

   // parameters common to all 
   if(dateStatus)  fAnaDate      = fParams["date"];
   if(blindStatus) fIsBlind      = (bool)( (int)fParams["blinding"]  ); 
   if(trlyStatus)  fTrolleyProbe = (int)fParams["trly-probe"]; 
   if(fxprStatus)  fFXPRListTag  = (int)fParams["fxpr-set"]; 

   if(runStatus){ 
      NRUNS = (int)fParams["nruns"];
      for(int i=0;i<NRUNS;i++){
	 fRunList.push_back( (int)fParams["run-list"][i] ); 
	 fRunLabel.push_back( fParams["run-label"][i] ); 
      }
   }
   
   // calibration: production 
   if( fType.compare("calib-prod")==0 ){
      if(axisStatus) fAxis          = (int)fParams["axis"];
      if(protStatus) fIsFreeProton  = (bool)( (int)fParams["free-proton-cor"] );  
      if(fitStatus)  fFitFunc       = fParams["fit"];
      if(swapStatus) fLoadSwapTime  = (bool)( (int)fParams["load-trly-swap-times"] ); 
      if(sccStatus)  fLoadSCCTime   = (bool)( (int)fParams["load-trly-scc-times"] ); 
      if(timeStatus) fUseTimeWeight = (bool)( (int)fParams["use-aba-time-weight"] ); 
   }
 
   // simple input format (i.e., rapid swapping data from 6/1)  
   if( fType.compare("is-simple")==0 ){
      fIsSimple                      = true;
      if(axisStatus) fAxis           = (int)fParams["axis"]; 
      if(anaStatus) fIsFullAnalysis  = (bool)( (int)fParams["full-ana"]  );    
      if(fitStatus) fFitFunc         = fParams["fit"]; 
   }else if( fType.compare("is-complex")==0 ){
      fIsSimple                      = false;
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
   std::cout << "Is blinded:         " << fIsBlind         << std::endl;
   std::cout << "Trolley probe:      " << fTrolleyProbe    << std::endl;
   if(fType.compare("calib-prod")==0){
	 std::cout << "Axis:               " << fAxis << " (" << axis << ")" << std::endl;
         std::cout << "Load TRLY swap time " << fLoadSwapTime << std::endl;
   }else{
      if(fIsSimple){
	 std::cout << "FXPR set:           " << fFXPRListTag     << std::endl;
      }else if(!fIsSimple){
	 std::cout << "Is full analysis:   " << fIsFullAnalysis  << std::endl;
	 std::cout << "Is final location:  " << fIsFinalLocation << std::endl;
	 std::cout << "Use P2P fit:        " << fUseP2PFit       << std::endl;
	 std::cout << "Analysis date:      " << fAnaDate         << std::endl;
	 std::cout << "Fit function:       " << fFitFunc         << std::endl;
	 std::cout << "Axis:               " << fAxis            << " (" << axis << ")" << std::endl;
      }
   }
   int N = fRunList.size(); 
   std::cout << "Number of runs:    " << N                << std::endl;
   for(int i=0;i<N;i++) std::cout << Form("run = %d, label = %s",fRunList[i],fRunLabel[i].c_str()) << std::endl;
   std::cout << "-----------------------------------------------------" << std::endl;
   return 0;
}
