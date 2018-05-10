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
   fIsFullAnalysis  = false; 
   fIsBlind         = false; 
   fUseP2PFit       = false; 
   fIsFinalLocation = false;
   fTrolleyProbe    = -1; 
   fAxis            = -1;
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
      std::cout << "[InputManager::GetLabelList]: No runs!" << std::endl;
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
      fIsFullAnalysis  = (bool)( (int)params["full-ana"]  );    
      fIsFinalLocation = (bool)( (int)params["final-loc"] );   
      fIsBlind         = (bool)( (int)params["blinding"]  ); 
      fUseP2PFit       = (bool)( (int)params["p2p-fit"]   ); 
      fAnaDate         = params["date"];
      fFitFunc         = params["fit"];  
      fTrolleyProbe    = (int)params["trly-probe"]; 
      fAxis            = (int)params["axis"];  
      NRUNS            = (int)params["nruns"];
      for(int i=0;i<NRUNS;i++){
	 fRunList.push_back( (int)params["run-list"][i] ); 
	 fRunLabel.push_back( params["run-label"][i] ); 
      }  
   }else{
      std::cout << "[InputManager::Load]: Cannot load data from " << inpath << std::endl;
      return 1;
   }
   return 0; 
}
//______________________________________________________________________________
int InputManager::Print(){
   char axis='N';
   if(fAxis==0) axis = 'x'; 
   if(fAxis==1) axis = 'y'; 
   if(fAxis==2) axis = 'z'; 
   int N = fRunList.size(); 
   std::cout << "------------------- Input Manager -------------------" << std::endl;
   std::cout << "Is full analysis:  " << fIsFullAnalysis  << std::endl;
   std::cout << "Is final location: " << fIsFinalLocation << std::endl;
   std::cout << "Is blinded:        " << fIsBlind         << std::endl;
   std::cout << "Use P2P fit:       " << fUseP2PFit       << std::endl;
   std::cout << "Analysis date:     " << fAnaDate         << std::endl;
   std::cout << "Fit function:      " << fFitFunc         << std::endl;
   std::cout << "Axis:              " << fAxis            << " (" << axis << ")" << std::endl;
   std::cout << "Trolley probe:     " << fTrolleyProbe    << std::endl;
   std::cout << "Number of runs:    " << N                << std::endl;
   for(int i=0;i<N;i++){
      std::cout << fRunList[i] << ", " << fRunLabel[i] << std::endl;
   }
   std::cout << "-----------------------------------------------------" << std::endl;
   return 0;
}
