#include "../include/Consolidate.h"
//______________________________________________________________________________
int CopyTrolleyProbe(std::vector<trolleyAnaEvent_t> x,std::vector<trolleyAnaEvent_t> &y){
   // copy a trolley probe vector to a new one 
   trolleyAnaEvent_t data;
   const int N = x.size();
   if(N==0){
      std::cout << "[CopyTrolleyProbe]: No data!" << std::endl;
      return 1;
   }
   for(int i=0;i<N;i++){
      for(int j=0;j<NUM_TRLY;j++){
         data.time[j] = x[i].time[j];
         data.freq[j] = x[i].freq[j];
         data.r[j]    = x[i].r[j];
         data.y[j]    = x[i].y[j];
         data.phi[j]  = x[i].phi[j];
         data.temp[j] = x[i].temp[j];
      }
      y.push_back(data);
   }
   return 0;
}
//______________________________________________________________________________
int ApplyBlindingTRLY(double blind,std::vector<trolleyAnaEvent_t> &x){
   // apply a blinding value to the PP data  
   plungingProbeAnaEvent_t data;
   int M=0;
   const int N = x.size();
   if(N==0){
      std::cout << "[ApplyBlindingPP]: No data!" << std::endl;
      return 1;
   }
   for(int i=0;i<N;i++){
      for(int j=0;j<NUM_TRLY;j++){
	 x[i].freq[j] += blind;
      }
   }
   return 0;
}
//______________________________________________________________________________
int ConsolidateTrolleyData(int startIndex,int method,std::vector<gm2field::trolleyTimeStamp_t> trlyTime,
                           std::vector<gm2field::trolleyPosition_t> trlyPos,
                           std::vector<gm2field::trolleyMonitor_t> trlyMon,
                           std::vector<gm2field::trolleyProbeFrequency_t> trlyFreq, 
                           std::vector<trolleyAnaEvent_t> &trlyEvent){
   // read in trolley data into a more useful single data struct
   trolleyAnaEvent_t data;
   const int NT = trlyTime.size();
   const int NP = trlyPos.size();
   const int NM = trlyMon.size();
   const int NF = trlyFreq.size();
   if(NT<=1 || NP<=1 || NM<=1 || NF<=1) return 1;

   double arg_phi=0;
   bool validEvent = false; 
   for(int i=startIndex;i<NT;i++){  
      for(int j=0;j<NUM_TRLY;j++){
	 if(trlyFreq[i].Frequency[j][method]>0){
            arg_phi = trlyPos[i].Phi[j][0]*gm2fieldUtil::Constants::DEG_PER_RAD;
            if(arg_phi<0) arg_phi += 360.; 
            validEvent = true; 
	    data.time[j] = trlyTime[i].GpsCycleStart[j];
	    data.freq[j] = trlyFreq[i].Frequency[j][method]; 
	    data.r[j]    = trlyPos[i].X[j];  
	    data.y[j]    = trlyPos[i].Y[j];  
	    data.phi[j]  = arg_phi; 
	    data.temp[j] = trlyMon[i].TemperatureExt[j];
	 }
         // if(data.freq[j]<20E+3){
	 //    std::cout << gm2fieldUtil::GetStringTimeStampFromUTC( data.time[j]/1E+9 ) << " " << data.freq[j] << std::endl; 
         // }
      }
      if(validEvent) trlyEvent.push_back(data); 
      validEvent = false; 
   }

   return 0;
}
//______________________________________________________________________________
int GetTrolleyData(std::string date,int runNumber,int freqMethod,std::vector<trolleyAnaEvent_t> &trlyEvent){
   // combine data from various structs into a single struct
   std::vector<gm2field::trolleyProbeFrequency_t> trlyFreq;
   int rc = gm2fieldUtil::RootHelper::GetTrolleyFrequencies(runNumber,trlyFreq);
   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   std::vector<gm2field::trolleyTimeStamp_t> trlyTime;
   rc = gm2fieldUtil::RootHelper::GetTrolleyTime(runNumber,trlyTime);
   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   std::vector<gm2field::trolleyPosition_t> trlyPos;
   rc = gm2fieldUtil::RootHelper::GetTrolleyPosition(runNumber,trlyPos);
   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   std::vector<gm2field::trolleyMonitor_t> trlyMon;
   rc = gm2fieldUtil::RootHelper::GetTrolleyMonitor(runNumber,trlyMon);
   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   int startIndex = 10; 
   if( date.compare("")!=0 ){
      // only if we have a valid date do we go look for a start index 
      startIndex = FindStartIndexTRLY(date,runNumber);
   } 

   // FIXME: special case...
   if(runNumber==3030) startIndex = 10; 
   if(runNumber==4911) startIndex = 35; 

   rc = ConsolidateTrolleyData(startIndex,freqMethod,trlyTime,trlyPos,trlyMon,trlyFreq,trlyEvent);

   return rc;
}
//______________________________________________________________________________
int CopyPlungingProbe(plungingProbeAnaEvent_t x,plungingProbeAnaEvent_t &y){
   // copy a plunging probe vector to a new one 
   int M = x.numTraces;
   y.numTraces = M;
   y.run = x.run;
   for(int i=0;i<M;i++){
      y.time[i]     = x.time[i];
      y.freq_LO[i]  = x.freq_LO[i];
      y.freq_RF[i]  = x.freq_RF[i];
      y.freq[i]     = x.freq[i];
      y.freq_err[i] = x.freq_err[i];
      y.r[i]        = x.r[i];
      y.y[i]        = x.y[i];
      y.phi[i]      = x.phi[i];
      y.temp[i]     = x.temp[i];
      y.temp_err[i] = x.temp_err[i];
   }
   return 0;
}
//______________________________________________________________________________
int CopyPlungingProbe(std::vector<plungingProbeAnaEvent_t> x,std::vector<plungingProbeAnaEvent_t> &y){
   // copy a plunging probe vector to a new one 
   plungingProbeAnaEvent_t data;
   int M=0;
   const int N = x.size();
   if(N==0){
      std::cout << "[CopyPlungingProbe]: No data!" << std::endl;
      return 1;
   }
   for(int i=0;i<N;i++){
      M = x[i].numTraces;
      data.numTraces = M;
      data.run = x[i].run;
      for(int j=0;j<M;j++){
         data.time[j]     = x[i].time[j];
         data.freq_LO[j]  = x[i].freq_LO[j];
         data.freq_RF[j]  = x[i].freq_RF[j];
         data.freq[j]     = x[i].freq[j];
         data.freq_err[j] = x[i].freq_err[j];
         data.r[j]        = x[i].r[j];
         data.y[j]        = x[i].y[j];
         data.phi[j]      = x[i].phi[j];
         data.temp[j]     = x[i].temp[j];
         data.temp_err[j] = x[i].temp_err[j];
      }
      y.push_back(data);
   }
   return 0;
}
//______________________________________________________________________________
int ConsolidatePPData(int method,std::vector<gm2field::plungingProbeFrequency_t> ppData,
                      std::vector<gm2field::plungingProbeInfo_t> ppInfo,
                      std::vector<plungingProbeAnaEvent_t> &ppEvent){

   // collect the PP data into a single struct
   // includes a variable "drift" that is a measure of field drift during the PP run 

   gm2fieldUtil::TemperatureSensor *tempSensor = new gm2fieldUtil::TemperatureSensor("PT1000");
   tempSensor->SetDataPath(TEMP_DIR); 

   plungingProbeAnaEvent_t data;
 
   std::string timeStamp;
 
   double theTemp=0; 
   int lastRun = ppInfo[0].FlayRunNumber;

   int cntr=0,M=0;
   const int N = ppInfo.size();

   for(int i=0;i<N;i++){
      theTemp   = tempSensor->GetTemperature( ppInfo[i].Temperature );
      timeStamp = gm2fieldUtil::GetStringTimeStampFromUTC( ppInfo[i].TimeStamp/1E+9 ); 
      // std::cout << "Processing NMR-DAQ run " << ppInfo[i].FlayRunNumber << std::endl;
      if( ppInfo[i].FlayRunNumber==lastRun ){ 
         // gather frequencies and temperatures for each NMR-DAQ run to average over  
         data.run           = ppInfo[i].FlayRunNumber; 
         data.time[cntr]    = ppInfo[i].TimeStamp; 
         data.r[cntr]       = ppInfo[i].R;  
         data.y[cntr]       = ppInfo[i].Y;  
         data.phi[cntr]     = ppInfo[i].Phi;  
         data.temp[cntr]    = theTemp;  
         data.freq[cntr]    = ppData[i].Frequency[method];  
         data.freq_LO[cntr] = ppData[i].FLO; 
         data.freq_RF[cntr] = ppData[i].F0;
	 cntr++;
      }else{
	 // done gathering, push back on the vector  
         data.numTraces  = cntr; 
         ppEvent.push_back(data);
	 // print to screen 
	 // for(int j=0;j<cntr;j++){ 
	 //    std::cout << "run = "  << Form("%d",data.run)  << " " 
	 //              << "time = " << Form("%s",gm2fieldUtil::GetStringTimeStampFromUTC(data.time[j]/1E+9 ).c_str() )  << " " 
	 //              << "temp = " << Form("%.3lf deg C",data.temp[j])  << " " 
	 //              << "r = "    << Form("%.3lf mm",data.r[j]  )      << " " 
	 //              << "y = "    << Form("%.3lf mm",data.y[j]  )      << " " 
	 //              << "phi = "  << Form("%.3lf mm",data.phi[j])      << " "
	 //              << "freq = " << Form("%.3lf Hz",data.freq[j])     << std::endl; 
	 //              // << "drift = " << Form("%.3lf Hz",data.freq_drift)  << std::endl;
	 // }
	 // std::cout << "----------------" << std::endl;
         // reset the counter 
         cntr = 0;
         // fill the data struct for this "event" since it's good for the next run 
         data.run           = ppInfo[i].FlayRunNumber; 
         data.time[cntr]    = ppInfo[i].TimeStamp; 
         data.r[cntr]       = ppInfo[i].R;  
         data.y[cntr]       = ppInfo[i].Y;  
         data.phi[cntr]     = ppInfo[i].Phi;  
         data.temp[cntr]    = theTemp;  
         data.freq[cntr]    = ppData[i].Frequency[method];  
         data.freq_LO[cntr] = ppData[i].FLO; 
         data.freq_RF[cntr] = ppData[i].F0;
	 cntr++;
      }
      lastRun = ppInfo[i].FlayRunNumber;
   }

   // catch the last run if we end on the same run as the last 
   data.numTraces  = cntr; 
   ppEvent.push_back(data);

   const int NN = ppEvent.size();
   std::cout << "NMR-DAQ runs: " << endl;
   for(int i=0;i<NN;i++) std::cout << ppEvent[i].run << endl;

   delete tempSensor; 

   return 0;
}
//______________________________________________________________________________
int GetPlungingProbeData(int runNumber,int freqMethod,std::vector<plungingProbeAnaEvent_t> &ppEvent){
   // combine the PP data into a vector of a single struct 
   std::vector<gm2field::plungingProbeFrequency_t> ppFreq;
   int rc = gm2fieldUtil::RootHelper::GetPPFrequencies(runNumber,ppFreq);
   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }
   std::vector<gm2field::plungingProbeInfo_t> ppInfo;
   rc = gm2fieldUtil::RootHelper::GetPPInfo(runNumber,ppInfo);
   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   rc = ConsolidatePPData(freqMethod,ppFreq,ppInfo,ppEvent);

   return rc;
}
//______________________________________________________________________________
int ModifyPlungingProbeData(int method,plungingProbeAnaEvent_t &Data){
   // replace the frequency values with those calculated by the NMR-ANA framework  
   int runNumber = Data.run;
   std::vector<nmrAnaEvent_t> inData;
   // std::cout << "Trying NMR-DAQ run " << runNumber << std::endl; 
   char inpath[512];
   sprintf(inpath,"./input/NMR-ANA/run-%05d/results_pulse-stats.dat",runNumber);
   int rc = ImportNMRANAData(inpath,inData);
   if(rc!=0){
      std::cout << "[ModifyPlungingProbeData]: No data for NMR-DAQ run " << runNumber << "!" << std::endl;
      return 1;
   }
  
   const int N = inData.size();
   if(N!=Data.numTraces){
      std::cout << "[ModifyPlungingProbeData]: Inconsistent number of traces for MIDAS and NMR-ANA!" << std::endl;
      std::cout << "NMR-DAQ run: " << Data.run       << std::endl;
      std::cout << "MIDAS:       " << Data.numTraces << " traces" << std::endl;
      std::cout << "NMR-ANA:     " << N              << " traces" << std::endl;
      return 1;
   }

   Data.numTraces = N; 
   for(int i=0;i<N;i++) Data.freq[i] = inData[i].freq[method];

   return 0;
}
//______________________________________________________________________________
int FilterPlungingProbeData(std::vector<int> subRun,
                            std::vector<plungingProbeAnaEvent_t> x,
                            std::vector<plungingProbeAnaEvent_t> &y){
   // copy over events from vector x into vector y 
   // that correspond to the proper runs listed in subRun (NMR-DAQ runs) 
   int nmrDAQ_run=0; 
   const int N = x.size(); 
   const int M = subRun.size(); 
   plungingProbeAnaEvent_t data;
   for(int i=0;i<N;i++){
      nmrDAQ_run = x[i].run;
      for(int j=0;j<M;j++){
	 if(nmrDAQ_run==subRun[j]){
	    CopyPlungingProbe(x[i],data);
            y.push_back(data); 
         }
      }
   } 
   return 0; 
}
//______________________________________________________________________________
int ApplyBlindingPP(double blind,std::vector<plungingProbeAnaEvent_t> &x){
   // apply a blinding value to the PP data  
   plungingProbeAnaEvent_t data;
   int M=0;
   const int N = x.size();
   if(N==0){
      std::cout << "[ApplyBlindingPP]: No data!" << std::endl;
      return 1;
   }
   for(int i=0;i<N;i++){
      M = x[i].numTraces;
      for(int j=0;j<M;j++){
	 x[i].freq[j] += blind;
      }
   }
   return 0;
}
