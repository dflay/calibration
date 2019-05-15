#include "../include/Cut.h"
//______________________________________________________________________________
Cut::Cut(std::string filepath,int verbosity){
   fVerbosity = verbosity;
   if( filepath.compare("UNKNOWN.json")!=0 ){
      LoadData(filepath);
   }else{
      std::cout << "[Cut]: No cut data loaded.  Use the LoadData function." << std::endl;
   }
   if(fVerbosity>0) Print(); 
}
//______________________________________________________________________________
Cut::~Cut(){

}
//______________________________________________________________________________
int Cut::LoadData(std::string filepath){
   if(fVerbosity>0) std::cout << "[Cut]: Loading data from: " << filepath << std::endl; 
   int rc = gm2fieldUtil::Import::ImportJSON(filepath,fData);
   if(rc!=0) std::cout << "[Cut]: Cannot load cut data!" << std::endl;
   if(fVerbosity>0) std::cout << "[Cut]: --> Done." << std::endl;
   return rc; 
}
//______________________________________________________________________________
bool Cut::CheckEvent_nmrAna(int run,int trace,int nzc,double ampl,double freq){

   int nzc_min     = fData["nzc"]["min"]; 
   double ampl_min = fData["ampl"]["min"]; 
   double ampl_max = fData["ampl"]["max"];
   double freq_min = fData["freq"]["min"]; 

   bool goodEvent = false; 

   // check for analysis characteristics 
   if( (freq>freq_min) && (ampl>ampl_min && ampl<ampl_max) && (nzc>nzc_min) ){
      goodEvent = true;
   }else{
      goodEvent = false;
   } 

   // if it's a specific run or event, reject it 
   int NTR=0,theRun=0,theTrace=0;
   const int NRUNS = fData["nmr-ana"]["runs"].size();

   for(int i=0;i<NRUNS;i++){
      theRun = fData["nmr-ana"]["runs"][i]["number"]; 
      NTR    = fData["nmr-ana"]["runs"][i]["trace"].size();
      for(int j=0;j<NTR;j++){
	 theTrace = fData["nmr-ana"]["runs"][i]["trace"][j]; 
	 if(run==theRun && trace==theTrace) goodEvent = false;  
      }
   }

   if(!goodEvent){
      std::cout << "[Cut]: Run " << run << " trace " << trace << " rejected! " << std::endl;
      std::cout << "       ampl = " << ampl << " nzc = " << nzc << " freq = " << freq << std::endl;
   }

   return goodEvent;
}
//______________________________________________________________________________
bool Cut::CheckEvent_trly(){
   return true; 
}
//______________________________________________________________________________
int Cut::Print(){
   int nzc_min     = fData["nzc"]["min"]; 
   double ampl_min = fData["ampl"]["min"]; 
   double ampl_max = fData["ampl"]["max"];
   double freq_min = fData["freq"]["min"]; 

   std::cout << "------------------ Cut Data ------------------" << std::endl;
   std::cout << "nzc (min):      " << nzc_min  << std::endl;
   std::cout << "freq (min):     " << freq_min << std::endl;
   std::cout << "ampl (min,max): " << ampl_min << ", " << ampl_max << std::endl;

   int NTR=0,theRun=0,theTrace=0;
   const int NRUNS = fData["nmr-ana"]["runs"].size();
   for(int i=0;i<NRUNS;i++){
      theRun = fData["nmr-ana"]["runs"][i]["number"]; 
      NTR    = fData["nmr-ana"]["runs"][i]["trace"].size();
      for(int j=0;j<NTR;j++){
	 theTrace = fData["nmr-ana"]["runs"][i]["trace"][j]; 
	 std::cout << "run " << theRun << " trace " << theTrace << std::endl;  
      }
   }
   std::cout << "----------------------------------------------" << std::endl;
   return 0;
}
