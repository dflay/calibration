#include "../include/Cut.h"
//______________________________________________________________________________
Cut::Cut(std::string filepath,int verbosity){
   fVerbosity = verbosity;
   if( filepath.compare("UNKNOWN.json")!=0 ){
      LoadData(filepath);
   }else{
      std::cout << "[Cut]: No cut data loaded.  Use the Cut::LoadData function." << std::endl;
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

   // check analysis parameters  
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
      std::cout << "[Cut]: Run " << run << " trace " << trace << " rejected! ";
      std::cout << "ampl = " << ampl << " nzc = " << nzc << " freq = " << freq << std::endl;
   }

   return goodEvent;
}
//______________________________________________________________________________
bool Cut::CheckEvent_trly(){
   return true; 
}
//______________________________________________________________________________
int Cut::FilterFXPRData(int runPeriod,int probe,
                        std::vector<averageFixedProbeEvent_t> in,std::vector<averageFixedProbeEvent_t> &out,
                        std::string inpath){

   // FILTER PP DATA
   // we may have multiple runs with different data sets in them 
   // this function filters data based on an input file 
   // type = data type we are analyzing.  Can be either dB or shim 

   const int N = in.size();

   // load json data with cut info  
   json jData;
   int rc = gm2fieldUtil::Import::ImportJSON(inpath,jData);
   if(rc!=0) return 1;

   // create probe key 
   char pr[20];
   sprintf(pr,"probe-%02d",probe);
   std::string probeStr = pr;

   auto it_key = jData.find(probeStr);  // this is an iterator 
   if (it_key!=jData.end() ){
      // not at the end of jData -- found the key 
   }else{
      std::cout << Form("[Cut::FilterFXPRData]: No key named: %s.  Keeping all data and returning.",probeStr.c_str()) << std::endl;
      for(int i=0;i<N;i++) out.push_back(in[i]);
      return 0;  // at the end of jData -- didn't find the key 
   }

   // get values
   std::string timeStr   = jData[probeStr]["fxpr"]["time"];
   std::string state_str = jData[probeStr]["fxpr"]["state"];

   // determine if STD or DST  
   bool isDST = false;
   int M = timeStr.size();
   std::string substr = timeStr.substr(M-3,M-1);
   if( substr.compare("CST")==0 ){
      isDST = false; // standard time
   }else{
      isDST = true;  // daylight savings time 
   }
   unsigned long TIME = gm2fieldUtil::GetUTCTimeStampFromString(timeStr,isDST);  

   std::cout << Form("[Cut::FilterFXPRData]: Applying the cut to run %d, probe %02d data: time = %s (%lu), fxpr = cut %s",
	             runPeriod,probe,timeStr.c_str(),TIME,state_str.c_str()) << std::endl;

   unsigned long theTime=0; 
   for(int i=0;i<N;i++){
      theTime = in[i].time/1E+9; // convert to sec 
      if( state_str.compare("before")==0 ){
	 if(theTime<TIME) out.push_back(in[i]); 
      }else if( state_str.compare("after")==0 ){
	 if(theTime>TIME) out.push_back(in[i]); 
     }  
   }
   return 0;
}
//______________________________________________________________________________
int Cut::FilterPPData(int runPeriod,int probe,std::string type,std::string axis,
                      std::vector<plungingProbeAnaEvent_t> in,std::vector<plungingProbeAnaEvent_t> &out,
                      std::string inpath){
   // FILTER PP DATA
   // we may have multiple runs with different data sets in them 
   // this function filters data based on an input file 
   // type = data type we are analyzing.  Can be either dB or shim 
   
   const int N = in.size();

   if(type.compare("dB")!=0 && type.compare("shim")!=0 ){
      std::cout << "[FilterData]: Invalid key type " << type << std::endl;
      return 1;
   }

   // load json data with cut info  
   json jData;
   int rc = gm2fieldUtil::Import::ImportJSON(inpath,jData);
   if(rc!=0) return 1;

   // create probe key 
   char pr[20];
   sprintf(pr,"probe-%02d",probe);
   std::string probeStr = pr;

   auto it_key = jData.find(probeStr);  // this is an iterator 
   if (it_key!=jData.end() ){
      // not at the end of jData -- found the key 
   }else{
      std::cout << Form("[Cut::FilterPPData]: No key named: %s.  Keeping all data and returning.",probeStr.c_str()) << std::endl;
      for(int i=0;i<N;i++) out.push_back(in[i]);  
      return 0;  // at the end of jData -- didn't find the key 
   }

   // get values
   std::string timeStr   = jData[probeStr]["time"];
   std::string state_str = jData[probeStr][type]["state"]; 
   std::string axis_str  = jData[probeStr][type]["axis"]; 

   if(axis.compare(axis_str)!=0){
      // not the axis of interest, keep all data
      std::cout << Form("[Cut::FilterPPData]: Axis of input data: %s, cut data %s.  Keeping all data and returning."
                        ,axis.c_str(),axis_str.c_str()) << std::endl;
      for(int i=0;i<N;i++) out.push_back(in[i]);  
      return 0;
   }

   // now we will actually make the cut; start with the timestamp 

   // determine if STD or DST  
   bool isDST = false;
   int M = timeStr.size();
   std::string substr = timeStr.substr(M-3,M-1);
   // std::cout << st << "," << substr << std::endl; 
   if( substr.compare("CST")==0 ){
      isDST = false; // standard time
   }else{
      isDST = true;  // daylight savings time 
   }
   unsigned long TIME = gm2fieldUtil::GetUTCTimeStampFromString(timeStr,isDST);

   // if(fVerbosity>0){
      std::cout << Form("[Cut::FilterPPData]: Applying the cut to run %d, probe %02d data: time = %s (%lu), %s = cut %s",
	                runPeriod,probe,timeStr.c_str(),TIME,type.c_str(),state_str.c_str()) << std::endl;
   // }

   // apply the cut 
   unsigned long theTime=0;

   M = in[0].numTraces; 

   for(int i=0;i<N;i++){
      theTime = in[i].time[M-1]/1E+9; // input time is in ns, convert to sec.  use time of last event in a run  
      if( state_str.compare("before")==0 ){
	 if( theTime<TIME ) out.push_back(in[i]);
      }else if( state_str.compare("after")==0 ){
	 if( theTime>TIME ) out.push_back(in[i]);
      }
   }

   std::cout << "[Cut::FilterPPData]: Size of output vector = " << out.size() << std::endl;

   return 0;
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
