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
   fUseAxis         = false; 
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
   json params;
   int NRUNS=0; 
   int rc = gm2fieldUtil::Import::ImportJSON(inpath,params);

   if(rc==0){
      fType = params["type"]; 
   }else{
      std::cout << "[InputManager::Load]: Cannot load data from " << inpath << std::endl;
      return 1; 
   }
   
   // parameters common to all 
   fAnaDate      = params["date"];
   fIsBlind      = (bool)( (int)params["blinding"]  ); 
   fTrolleyProbe = (int)params["trly-probe"]; 
   NRUNS = (int)params["nruns"];
   for(int i=0;i<NRUNS;i++){
      fRunList.push_back( (int)params["run-list"][i] ); 
      fRunLabel.push_back( params["run-label"][i] ); 
   }

   // calibration: production 
   if( fType.compare("calib-prod")==0 ){
      if(fUseAxis) fAxis = (int)params["axis"];
      fFXPRListTag = (int)params["fxpr-set"];   
   }

   // simple input format (i.e., rapid swapping data from 6/1)  
   if( fType.compare("is-simple")==0 ){
      fIsSimple        = true;
      fFXPRListTag     = (int)params["fxpr-set"];   
   }else if( fType.compare("is-complex")==0 ){
      fIsSimple        = false;
      fIsFullAnalysis  = (bool)( (int)params["full-ana"]  );    
      fIsFinalLocation = (bool)( (int)params["final-loc"] );   
      fUseP2PFit       = (bool)( (int)params["p2p-fit"]   ); 
      fFitFunc         = params["fit"]; 
      fAxis            = (int)params["axis"];
      fFXPRListTag     = (int)params["fxpr-set"];   
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
   std::cout << "Is blinded:        " << fIsBlind         << std::endl;
   std::cout << "Trolley probe:     " << fTrolleyProbe    << std::endl;
   if(fType.compare("calib-prod")==0){
	 if(fUseAxis) std::cout << "Axis:              " << fAxis << " (" << axis << ")" << std::endl;
   }else{
      if(fIsSimple){
	 std::cout << "FXPR set:          " << fFXPRListTag     << std::endl;
      }else if(!fIsSimple){
	 std::cout << "Is full analysis:  " << fIsFullAnalysis  << std::endl;
	 std::cout << "Is final location: " << fIsFinalLocation << std::endl;
	 std::cout << "Use P2P fit:       " << fUseP2PFit       << std::endl;
	 std::cout << "Analysis date:     " << fAnaDate         << std::endl;
	 std::cout << "Fit function:      " << fFitFunc         << std::endl;
	 std::cout << "Axis:              " << fAxis            << " (" << axis << ")" << std::endl;
      }
   }
   int N = fRunList.size(); 
   std::cout << "Number of runs:    " << N                << std::endl;
   for(int i=0;i<N;i++) std::cout << Form("run = %d, label = %s",fRunList[i],fRunLabel[i].c_str()) << std::endl;
   std::cout << "-----------------------------------------------------" << std::endl;
   return 0;
}
