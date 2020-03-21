#include "../include/Cut.h"
//______________________________________________________________________________
Cut::Cut(std::string filepath,int verbosity){
   // kLOWER = 0;
   // kUPPER = 1;
   // kRANGE = 2;
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
   // type = data type we are analyzing.  Data type must be dB or shim 
   
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
   std::string timeStr   = jData[probeStr][type]["time"];
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

   // check if the time string is actually more than one timestamp (i.e., defines a cut range)
   int cutType=-1;
   bool isCutRange=false;
   std::vector<std::string> tsv;
   std::string::size_type n;
   n = timeStr.find(',');  // search for a comma character from beginning of string  
   if(n==std::string::npos){
      // no comma, this is a single timestamp
      tsv.push_back(timeStr);
   }else{
      // comma found, is a cut range 
      isCutRange = true;
      rc = SplitString(",",timeStr,tsv); // split to a vector 
   }

   // verify cut type  
   if( state_str.compare("lower")==0 ){
      cutType = kLOWER;
   }else if( state_str.compare("upper")==0){
      cutType = kUPPER;
   }else if( state_str.compare("range")==0){
      cutType = kRANGE;
   }

   // consistency check
   if(isCutRange==true && cutType!=2){
      std::cout << "[Cut::FilterPPData]: ERROR! Time vector size is > 1, state is not range!" << std::endl;
      std::cout << "                                   Keeping all data and returning." << std::endl;
      for(int i=0;i<N;i++) out.push_back(in[i]);  
      return 0;
   }

   // convert to UTC
   int M=0;
   bool isDST=false;
   unsigned long long TIME=0;
   const int NT = tsv.size();
   std::vector<unsigned long long> timeCut;
   std::string substr;
   for(int i=0;i<NT;i++){
      // determine if STD or DST  
      M = tsv[i].size();
      substr = tsv[i].substr(M-3,M-1);
      if( substr.compare("CST")==0 ){
         isDST = false; // standard time
      }else{
         isDST = true;  // daylight savings time 
      }
      TIME = 1E+9*gm2fieldUtil::GetUTCTimeStampFromString(tsv[i],isDST); // convert to ns
      timeCut.push_back(TIME);
   }

   std::cout << Form("[Cut::FilterPPData]: Applying the cut to run %d, probe %02d data: time = %s, %s cut = %s",
	             runPeriod,probe,timeStr.c_str(),type.c_str(),state_str.c_str()) << std::endl;

   // apply the cut 
   unsigned long theTime=0;

   plungingProbeAnaEvent_t d; 
   std::vector<plungingProbeAnaEvent_t> data; 

   // check every trace
   bool passCut=false;
   int k=0,l=0; 
   for(int i=0;i<N;i++){
      M = in[i].numTraces;
      for(int j=0;j<M;j++){
	 theTime            = in[i].time[j];
         d.run              = in[i].run;
         d.midasRun         = in[i].midasRun;
	 d.time[l]          = theTime;
	 d.t2Time[l]        = in[i].t2Time[j];
	 d.nzc[l]           = in[i].nzc[j];
	 d.traceNumber[l]   = in[i].traceNumber[j];
	 d.channelNumber[l] = in[i].channelNumber[j];
	 d.freq_LO[l]       = in[i].freq_LO[j];
	 d.freq_RF[l]       = in[i].freq_RF[j];
	 d.freq[l]          = in[i].freq[j];
	 d.freq_err[l]      = in[i].freq_err[j];
	 d.r[l]             = in[i].r[j];
	 d.y[l]             = in[i].y[j];
	 d.phi[l]           = in[i].phi[j];
	 d.temp[l]          = in[i].temp[j];
	 d.temp_err[l]      = in[i].temp_err[j]; 
	 if(cutType==kLOWER && theTime>timeCut[0]){
	    // advance num traces k, and trace index l 
	    passCut = true; 
	 } 
	 if(cutType==kUPPER && theTime<timeCut[0]){
	    // advance num traces k, and trace index l 
	    passCut = true; 
	 } 
	 if(cutType==kRANGE && theTime>timeCut[0] && theTime<timeCut[1]){
	    // advance num traces k, and trace index l
	    passCut = true; 
	 }
         if(passCut){
	    k++;
	    l++; 
	    d.numTraces = k; 
         }
	 // reset for next event
	 passCut = false; 
      }
      // push back on the data which is cut according to the above  
      if(k>0) data.push_back(d);  // push back only if we've accumulated events 
      // reset for next PP-DAQ run  
      k = 0;
      l = 0;
   }

   // fill output vector
   int ND = data.size();
   for(int i=0;i<ND;i++) out.push_back(data[i]); 

   // M = in[0].numTraces;
   // for(int i=0;i<N;i++){
   //    theTime = in[i].time[M-1];  // use last time stamp of PP-DAQ runs 
   //    if(cutType==kLOWER) if(theTime>timeCut[0]) out.push_back(in[i]);  
   //    if(cutType==kUPPER) if(theTime<timeCut[0]) out.push_back(in[i]);  
   //    if(cutType==kRANGE) if(theTime>timeCut[0]&&theTime<timeCut[1]) out.push_back(in[i]);  
   // }

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
