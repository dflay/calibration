#include "../include/CustomImport.h"
#include "TRLYCoordinates.C"
#include "Cut.C"
#include "Logger.C"
//______________________________________________________________________________
int SetDataFileParameters(std::string version,std::string &fileName,std::string &dataPath){
   // determine file name prefix and data path based upon version tag 
   if( version.compare("run-1")==0 ){
      fileName = "default";
      dataPath = "default";
   }else if( version.compare("run-2")==0 ){
      fileName = "default";
      dataPath = "/mnt/nfs/g2field-server-2/newg2/DataProduction/Nearline/ArtTFSDir";  // use this for now... 
   }else if( version.compare("v09_04")==0 ){
      // grid production 
      fileName = "FieldPlainRootOutput_";
      dataPath = "/home/newg2/DataProduction/Offline/ArtTFSDir/v09_04";
   }else if( version.compare("calib")==0 ){
      // calibration (new production on g2field-server-2) 
      fileName = "default";
      dataPath = "/data2/newg2/DataProduction/Nearline/ArtTFSDir"; 
   }else if( version.compare("v9_21_01")==0 ){
      fileName = "FieldPlainRootOutput_";
      dataPath = "/data1/newg2/DataProduction/Offline/ArtTFSDir/v9_21_01";  
   }else if( version.compare("nearline")==0 ){
      fileName = "FieldGraphOut"; 
      dataPath = "/data1/newg2/DataProduction/Nearline/ArtTFSDir"; 
   }else if( version.compare("nearline-run-2")==0 ){
      fileName = "FieldGraphOut"; 
      dataPath = "/mnt/nfs/g2field-server-2/newg2/DataProduction/Nearline/ArtTFSDir"; 
   }else{
      std::cout << "No match for production tag: " << version << std::endl;
      return 1;
   }
   return 0;
}
//______________________________________________________________________________
int GetFixedProbeData(int run,int method,int probe,std::vector<fixedProbeEvent_t> &data,std::string version){
   // gather AVERAGED fixed probe data according to the probe list input

   std::vector<gm2field::fixedProbeFrequency_v1_t> fxpr; 
   std::vector<gm2field::fixedProbeFrequency_t> fxpr_new;

   std::vector<gm2field::fixedProbeHeader_v1_t> fxprHeader;
   std::vector<gm2field::fixedProbeHeader_t> fxprHeader_new;  

   int rc=0,N=0;
   if( version.compare("nearline")!=0 ){
      rc = GetFixedProbeFrequencies<gm2field::fixedProbeFrequency_t>(run,fxpr_new,version);
      rc = GetFixedProbeHeader<gm2field::fixedProbeHeader_t>(run,fxprHeader_new,version);
      N = fxpr_new.size(); 
   }else{
      rc = GetFixedProbeFrequencies<gm2field::fixedProbeFrequency_v1_t>(run,fxpr,version); 
      rc = GetFixedProbeHeader<gm2field::fixedProbeHeader_v1_t>(run,fxprHeader,version);
      N = fxpr.size(); 
   }

   fixedProbeEvent_t dataPt; 

   char arg_lay_id='t',arg_rad_id='t',arg_yoke_id='t';
   int arg_stat_id=0,arg_azi_id=0,arg_probe_id=0;
   unsigned long long arg_t=0,mean_t=0;
   double arg_f=0;

   for(int i=0;i<N;i++){
      if( version.compare("nearline")!=0 ){
	 arg_t       = fxpr_new[i].GpsTimeStamp[probe];
	 arg_f       = fxpr_new[i].Frequency[probe][method];
	 arg_azi_id  = fxprHeader_new[i].AziId[probe];
	 arg_lay_id  = fxprHeader_new[i].LayerId[probe];
	 arg_yoke_id = fxprHeader_new[i].YokeId[probe];
	 arg_rad_id  = fxprHeader_new[i].RadId[probe];
	 arg_stat_id = fxprHeader_new[i].StationId[probe];
      }else{
	 arg_t        = fxpr[i].GpsTimeStamp[probe];
	 arg_f        = fxpr[i].Frequency[probe][method];
	 arg_azi_id   = fxprHeader[i].AziId[probe];
	 arg_lay_id   = fxprHeader[i].LayerId[probe];
	 arg_yoke_id  = fxprHeader[i].YokeId[probe];
	 arg_rad_id   = fxprHeader[i].RadId[probe];
         arg_probe_id = fxprHeader[i].ProbeId[probe];  
	 arg_stat_id  = 0;
      }
      dataPt.time      = arg_t; 
      dataPt.freq      = arg_f;
      dataPt.yokeID    = arg_yoke_id; 
      dataPt.layerID   = arg_lay_id; 
      dataPt.radID     = arg_rad_id;
      dataPt.aziID     = arg_azi_id; 
      dataPt.probeID   = arg_probe_id;
      dataPt.stationID = arg_stat_id;  
      data.push_back(dataPt); 
   }

   return 0;
}
//______________________________________________________________________________
int GetFixedProbeData_avg(int run,int method,std::vector<int> probe,
                          std::vector<averageFixedProbeEvent_t> &data,std::string version,
                          bool subtractDrift,int period,unsigned long long t0){
   // gather AVERAGED fixed probe data according to the probe list input

   std::vector<gm2field::fixedProbeFrequency_v1_t> fxpr; 
   std::vector<gm2field::fixedProbeFrequency_t> fxpr_new; 

   int rc=0,N=0;
   if( version.compare("nearline")!=0 ){
      rc = GetFixedProbeFrequencies<gm2field::fixedProbeFrequency_t>(run,fxpr_new,version);
      N = fxpr_new.size(); 
   }else{
      rc = GetFixedProbeFrequencies<gm2field::fixedProbeFrequency_v1_t>(run,fxpr,version); 
      N = fxpr.size(); 
   }

   averageFixedProbeEvent_t dataPt; 
   const int NPR = probe.size();
   std::cout << "[GetFixedProbeData_avg]: Using " << NPR << " probes" << std::endl;

   unsigned long long arg_t=0;
   double t0d = t0/1E+9; 
   double arg_f=0,mean_t=0,mean_f=0,stdev_f=0,f0=0;
   std::vector<double> tt,ff,T,F; 

   int k=0;
   for(int i=0;i<N;i++){
      for(int j=0;j<NPR;j++){
         k = probe[j];
	 if( version.compare("nearline")!=0 ){
	    arg_t = fxpr_new[i].GpsTimeStamp[k];
	    arg_f = fxpr_new[i].Frequency[k][method];
	 }else{
	    arg_t = fxpr[i].GpsTimeStamp[k];
	    arg_f = fxpr[i].Frequency[k][method];
	 }
	 T.push_back(arg_t/1E+9);
	 F.push_back(arg_f);
      }
      // compute averages 
      mean_t  = gm2fieldUtil::Math::GetMean<double>(T);
      mean_f  = gm2fieldUtil::Math::GetMean<double>(F);
      stdev_f = gm2fieldUtil::Math::GetStandardDeviation<double>(F);
      // store results and clear vectors
      if(mean_t>=t0d){ 
	 tt.push_back(mean_t); 
	 ff.push_back(mean_f);
      }
      T.clear();
      F.clear();
   }

   // construct moving average over time 
   double x,y;
   std::vector<double> TT,FF;
   gm2fieldUtil::Algorithm::MovingAverage<double> *mvAvg1 = new gm2fieldUtil::Algorithm::MovingAverage<double>(period);
   gm2fieldUtil::Algorithm::MovingAverage<double> *mvAvg2 = new gm2fieldUtil::Algorithm::MovingAverage<double>(period);
   const int NN = tt.size();
   for(int i=0;i<NN;i++){
      x = mvAvg1->GetUpdatedAverage(tt[i]);
      y = mvAvg2->GetUpdatedAverage(ff[i]);
      TT.push_back(x);
      FF.push_back(y);
   }

   delete mvAvg1; 
   delete mvAvg2; 

   // determine slope of data 
   double slope=0,intercept=0,r=0;
   rc = gm2fieldUtil::Math::LeastSquaresFitting(TT,FF,intercept,slope,r);

   // subtract off slope
   double arg=0;
   for(int i=0;i<NN;i++){
      arg = FF[i] - (intercept + slope*TT[i]);
      if(subtractDrift){
	 dataPt.freq = arg;
      }else{
	 dataPt.freq = FF[i]; 
      } 
      dataPt.time    = (unsigned long long)(TT[i]*1E+9);  // this needs to be unsigned long long 
      dataPt.freqErr = 0;
      data.push_back(dataPt);  
   }

   // first data sample is the point to which we want to reference 
   // double F0 = data[0].freq; 
   // for(int i=0;i<NN;i++) data[i].freq -= F0;

   return 0;
} 
//______________________________________________________________________________
int GetSurfaceCoilData(int run,std::vector<surfaceCoilEvent_t> &data,std::string version){

   std::vector<gm2field::surfaceCoils_v1_t> scc; 
   std::vector<gm2field::surfaceCoils_t> scc_new; 

   int rc=0,N=0;

   if( version.compare("nearline")!=0 ){
      rc = GetSurfaceCoil<gm2field::surfaceCoils_t>(run,scc_new,version);
      N  = scc_new.size(); 
   }else{
      rc = GetSurfaceCoil<gm2field::surfaceCoils_v1_t>(run,scc,version); 
      N  = scc.size(); 
   }

   if(rc!=0){
      std::cout << "[GetSurfaceCoilData]: Cannot read in data!" << std::endl;
      return 1; 
   }
 
   surfaceCoilEvent_t dataPt;
   for(int i=0;i<N;i++){
      for(int j=0;j<NUM_SCC;j++){
	 if( version.compare("nearline")!=0 ){
	    dataPt.BotTime[j]     = scc_new[i].BotTime[j]; 
	    dataPt.TopTime[j]     = scc_new[i].TopTime[j]; 
	    dataPt.BotCurrents[j] = scc_new[i].BotCurrents[j]; 
	    dataPt.TopCurrents[j] = scc_new[i].TopCurrents[j]; 
	    dataPt.BotTemps[j]    = scc_new[i].BotTemps[j]; 
	    dataPt.TopTemps[j]    = scc_new[i].TopTemps[j]; 
	 }else{
	    dataPt.BotTime[j]     = scc[i].BotTime[j]; 
	    dataPt.TopTime[j]     = scc[i].TopTime[j]; 
	    dataPt.BotCurrents[j] = scc[i].BotCurrents[j]; 
	    dataPt.TopCurrents[j] = scc[i].TopCurrents[j]; 
	    dataPt.BotTemps[j]    = scc[i].BotTemps[j]; 
	    dataPt.TopTemps[j]    = scc[i].TopTemps[j]; 
	 }
      }
      // azi coils
      for(int j=0;j<4;j++){
	 if( version.compare("nearline")!=0 ){
	    dataPt.AzCurrents[j]  = scc_new[i].AzCurrents[j]; 
	    dataPt.AzTemps[j]     = scc_new[i].AzTemps[j]; 
	 }else{
	    dataPt.AzCurrents[j]  = scc[i].AzCurrents[j]; 
	    dataPt.AzTemps[j]     = scc[i].AzTemps[j]; 
	 }
      }
      data.push_back(dataPt); 
   }

   return 0;
}
//______________________________________________________________________________
int GetTrolleyData(int run,int method,std::vector<trolleyAnaEvent_t> &trlyEvent,std::string version){

   int rc=0;

   std::vector<gm2field::trolleyProbeFrequency_v1_t> trlyFreq;
   std::vector<gm2field::trolleyProbeFrequency_t> trlyFreq_new;

   std::vector<gm2field::trolleyPosition_v1_t> trlyPos;
   std::vector<gm2field::trolleyPosition_t> trlyPos_new;

   std::vector<gm2field::trolleyTimeStamp_t> trlyTime; // only one version of trolley time 
   std::vector<gm2field::trolleyMonitor_t> trlyMon;    // only one version of trolley monitor 

   int NT=0,NP=0,NM=0,NF=0;

   // these items don't change based upon data production    
   rc = GetTrolleyTimeStamps<gm2field::trolleyTimeStamp_t>(run,trlyTime,version);
   rc = GetTrolleyMonitor<gm2field::trolleyMonitor_t>(run,trlyMon,version);
   NT = trlyTime.size();
   NM = trlyMon.size();

   if(version.compare("nearline")==0){
      rc = GetTrolleyFrequencies<gm2field::trolleyProbeFrequency_v1_t>(run,trlyFreq,version);
      rc = GetTrolleyPosition<gm2field::trolleyPosition_v1_t>(run,trlyPos,version);
      NP = trlyPos.size();
      NF = trlyFreq.size();
   }else{
      rc = GetTrolleyFrequencies<gm2field::trolleyProbeFrequency_t>(run,trlyFreq_new,version);
      rc = GetTrolleyPosition<gm2field::trolleyPosition_t>(run,trlyPos_new,version);
      NP = trlyPos_new.size();
      NF = trlyFreq_new.size();
   }

   if(NT<=1 || NP<=1 || NM<=1 || NF<=1) return 1;

   trolleyAnaEvent_t data;

   int startIndex=22;  // for continuous mode 
   // if( date.compare("")!=0 ){
   //    // only if we have a valid date do we go look for a start index 
   //    startIndex = FindStartIndexTRLY(date,runNumber);
   // }

   // special case...
   if(run==3030) startIndex = 10;
   if(run==4911) startIndex = 35;

   double arg_phi=0,arg_freq=0,arg_time=0,arg_temp=0;
   bool validEvent = false;
   for(int i=startIndex;i<NT;i++){
      for(int j=0;j<NUM_TRLY;j++){
         if(version.compare("nearline")==0){
            arg_phi  = trlyPos[i].Phi[j][0]*gm2fieldUtil::Constants::DEG_PER_RAD;
            arg_freq = trlyFreq[i].Frequency[j][method];
         }else{
            // arg_phi  = trlyPos_new[i].Phi[j][0]*gm2fieldUtil::Constants::DEG_PER_RAD;
            arg_phi  = trlyPos_new[i].Phi[j][0];
            arg_freq = trlyFreq_new[i].Frequency[j][method];
         }
         arg_time = trlyTime[i].GpsCycleStart[j];
         arg_temp = trlyMon[i].TemperatureExt[j];
         if(arg_phi<0) arg_phi += 360.;
         if(arg_freq>0){
            validEvent   = true;
            data.time[j] = arg_time;
            data.freq[j] = arg_freq;
            data.r[j]    = GetTrolleyProbeTransverseCoordinate(j,"r"); 
            data.y[j]    = GetTrolleyProbeTransverseCoordinate(j,"y"); 
            data.phi[j]  = arg_phi;
            data.temp[j] = arg_temp;
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
int GetPlungingProbeData(int run,int prMethod,int ppMethod,std::vector<plungingProbeAnaEvent_t> &data,
                         std::string version,std::string nmrAnaVersion,std::string cutFile,bool useNMRANA){

   // prMethod = production method (from art analysis) 
   // ppMethod = UMass analysis method 

   int rc=0;

   std::vector<gm2field::plungingProbeFrequency_v1_t> ppData;
   std::vector<gm2field::plungingProbeFrequency_t> ppData_new;
   std::vector<gm2field::plungingProbeInfo_v1_t> ppInfo;
   std::vector<gm2field::plungingProbeInfo_t> ppInfo_new;

   int N=0;
   int lastRun=0;
   if( version.compare("nearline")!=0 ){
      rc      = GetPlungingProbeFrequencies<gm2field::plungingProbeFrequency_t>(run,ppData_new,version);
      rc      = GetPlungingProbeInfo<gm2field::plungingProbeInfo_t>(run,ppInfo_new,version);
      lastRun = ppInfo_new[0].FlayRunNumber;
      N       = ppInfo_new.size();
   }else{
      rc      = GetPlungingProbeFrequencies<gm2field::plungingProbeFrequency_v1_t>(run,ppData,version);
      rc      = GetPlungingProbeInfo<gm2field::plungingProbeInfo_v1_t>(run,ppInfo,version);
      lastRun = ppInfo[0].FlayRunNumber;
      N       = ppInfo.size();
   }
   if(rc!=0) return 1;

   if( TMath::Abs(lastRun) > 99999 || lastRun<=0 ){
      std::cout << "[GetPlungingProbeData]: Invalid starting run " << lastRun << std::endl;
      return 1;
   }

   // collect the PP data into a single struct
   // includes a variable "drift" that is a measure of field drift during the PP run 

   gm2fieldUtil::Temperature::Sensor *tempSensor = new gm2fieldUtil::Temperature::Sensor("PT1000");
   tempSensor->SetDataPath(TEMP_DIR);

   plungingProbeAnaEvent_t subRun;

   std::string timeStamp;

   int theRun=0,NZero=0,theCh=0;
   double theAmpl=0,theTemp=0,theR=0,theY=0,thePhi=0,theFreq=0,theLO=0,theRF=0;
   unsigned long long theTime=0;

   int cntr=0,M=0;

   for(int i=0;i<N;i++){
      if( version.compare("nearline")!=0 ){
         theTemp = tempSensor->GetTemperature( ppInfo_new[i].Temperature );
         theAmpl = ppInfo_new[i].Amplitude;
         NZero   = ppInfo_new[i].NZero;   
         theFreq = ppData_new[i].Frequency[prMethod]; 
         theLO   = ppData_new[i].FLO; 
         theRF   = ppData_new[i].F0; 
         theTime = ppInfo_new[i].TimeStamp;
         theRun  = ppInfo_new[i].FlayRunNumber;
         theR    = ppInfo_new[i].R;
         theY    = ppInfo_new[i].Y;
         thePhi  = ppInfo_new[i].Phi;
      }else{
         theTemp = tempSensor->GetTemperature( ppInfo[i].Temperature );
         theAmpl = ppInfo[i].Amplitude;
         NZero   = 0; 
         theFreq = ppData[i].Frequency[prMethod]; 
         theLO   = ppData[i].FLO; 
         theRF   = ppData[i].F0; 
         theTime = ppInfo[i].TimeStamp;
         theRun  = ppInfo[i].FlayRunNumber;
         theR    = ppInfo[i].R;
         theY    = ppInfo[i].Y;
         thePhi  = ppInfo[i].Phi;
      }
      timeStamp = gm2fieldUtil::GetStringTimeStampFromUTC(theTime/1E+9);
      // std::cout << "EVENT " << i << " TIME " << timeStamp << " RUN " << theRun << std::endl;
      if(theRun==lastRun){
         // gather frequencies and temperatures for each NMR-DAQ run to average over 
         subRun.run                 = theRun;
         subRun.midasRun            = run; 
	 subRun.traceNumber[cntr]   = cntr+1;  
	 subRun.channelNumber[cntr] = theCh;  
         subRun.time[cntr]          = theTime;
         subRun.r[cntr]             = theR;
         subRun.y[cntr]             = theY;
         subRun.phi[cntr]           = thePhi;
         subRun.temp[cntr]          = theTemp;
         subRun.freq[cntr]          = theFreq;
         subRun.freq_LO[cntr]       = theLO;
         subRun.freq_RF[cntr]       = theRF;
         subRun.nzc[cntr]           = NZero; 
         subRun.ampl[cntr]          = theAmpl; 
         cntr++;
      }else{
         // done gathering, push back on the vector  
         subRun.numTraces = cntr;
         data.push_back(subRun);
         // print to screen 
         // for(int j=0;j<cntr;j++){  
         //    std::cout << "run = "  << Form("%d",subRun.run)  << " " 
         //              << "time = " << Form("%s",gm2fieldUtil::GetStringTimeStampFromUTC(subRun.time[j]/1E+9 ).c_str() )  << " "
         //              << "temp = " << Form("%.3lf deg C",subRun.temp[j])  << " "
         //              << "r = "    << Form("%.3lf mm",subRun.r[j]  )      << " "
         //              << "y = "    << Form("%.3lf mm",subRun.y[j]  )      << " "
         //              << "phi = "  << Form("%.3lf mm",subRun.phi[j])      << " "
         //              << "freq = " << Form("%.3lf Hz",subRun.freq[j])     << std::endl;
         // }
         // std::cout << "----------------" << std::endl;
         // reset the counter 
         cntr = 0;
         // fill the data struct for this "event" since it's good for the next run 
         subRun.run                 = theRun; 
         subRun.midasRun            = run; 
	 subRun.traceNumber[cntr]   = cntr+1;  
	 subRun.channelNumber[cntr] = theCh;  
         subRun.time[cntr]          = theTime;
         subRun.r[cntr]             = theR;
         subRun.y[cntr]             = theY;  
         subRun.phi[cntr]           = thePhi; 
         subRun.temp[cntr]          = theTemp;
	 subRun.freq[cntr]          = theFreq;
	 subRun.freq_LO[cntr]       = theLO; 
	 subRun.freq_RF[cntr]       = theRF;
         subRun.nzc[cntr]           = NZero; 
         subRun.ampl[cntr]          = theAmpl; 
         cntr++;
      }
      lastRun = theRun;
   }

   // catch the last run if we end on the same run as the last 
   subRun.numTraces = cntr;
   data.push_back(subRun);

   const int NN = data.size();
   if(useNMRANA){
      std::cout << "[GetPlungingProbeData]: Using NMR-ANA results... " << endl;
      for(int i=0;i<NN;i++){
	 rc = ModifyPlungingProbeData(ppMethod,data[i],nmrAnaVersion,cutFile);
      }
      std::cout << "[GetPlungingProbeData]: --> Done." << std::endl;
   }

   delete tempSensor;

   return 0;
}
//______________________________________________________________________________
int ModifyPlungingProbeData(int method,plungingProbeAnaEvent_t &data,std::string nmrAnaVersion,std::string cutFile){
   // replace the frequency values with those calculated by the NMR-ANA framework 
   char msg[400];
   std::string timeStamp,MSG;

   int runNumber = data.run;
   std::vector<nmrAnaEvent_t> inData;
   // std::cout << "Trying NMR-DAQ run " << runNumber << std::endl; 
   char inpath[512];
   sprintf(inpath,"./input/NMR-ANA/%s/run-%05d/results_pulse-stats.dat",nmrAnaVersion.c_str(),runNumber);
   int rc = ImportNMRANAData(inpath,inData,cutFile);
   if(rc!=0){
      sprintf(msg,"[ModifyPlungingProbeData]: No data for NMR-DAQ run %d!",runNumber);
      MSG = msg; 
      rc = Logger::PrintMessage(Logger::kINFO,"default",MSG,'a'); 
      return 1;
   }

   bool isDiff = false;
   const int N = inData.size();
   if(N!=data.numTraces){
      sprintf(msg,"[ModifyPlungingProbeData]: Inconsistent number of traces for MIDAS and NMR-ANA!");
      MSG = msg; 
      rc = Logger::PrintMessage(Logger::kINFO,"default",MSG,'a'); 
      sprintf(msg,"NMR-DAQ run: %d",data.run);
      MSG = msg; 
      rc = Logger::PrintMessage(Logger::kINFO,"default",MSG,'a'); 
      sprintf(msg,"MIDAS:       %d traces",data.numTraces);
      MSG = msg; 
      rc = Logger::PrintMessage(Logger::kINFO,"default",MSG,'a'); 
      sprintf(msg,"NMR-ANA:     %d traces",N             );
      MSG = msg; 
      rc = Logger::PrintMessage(Logger::kINFO,"default",MSG,'a'); 
      // return 1;
      isDiff = true;
   }

   int M = data.numTraces; 

   if(isDiff){ 
      MSG = "----------- INCOMING PP DATA ------------";
      rc = Logger::PrintMessage(Logger::kINFO,"default",MSG,'a'); 

      for(int i=0;i<M;i++){
	 timeStamp = gm2fieldUtil::GetStringTimeStampFromUTC( data.time[i]/1E+9 ); 
	 sprintf(msg,"run %05d, trace %02d: time = %s, ampl = %.3lf, freq = %.3lf, temp = %.3lf",
	       data.run,data.traceNumber[i],timeStamp.c_str(),data.ampl[i],data.freq[i],data.temp[i]);
	 MSG = msg;
	 rc = Logger::PrintMessage(Logger::kINFO,"default",MSG,'a'); 
      } 
   }

   // now we may have excluded traces here, so we have to be careful 
   int k=0;
   for(int i=0;i<M;i++){       // loop over MIDAS data 
      for(int j=0;j<N;j++){    // loop over NMR-ANA data
	 if(data.traceNumber[i]==inData[j].pulse){
	    // we have a match -- update the struct
            // data.run            = data.run;  
	    data.traceNumber[k] = data.traceNumber[i];  
	    data.time[k]        = data.time[i];
	    data.r[k]           = data.r[i];
	    data.y[k]           = data.y[i];  
	    data.phi[k]         = data.phi[i]; 
	    data.temp[k]        = data.temp[i];
	    data.freq[k]        = inData[j].freq[method];
	    data.freq_LO[k]     = data.freq_LO[i]; 
	    data.freq_RF[k]     = data.freq_RF[i];
	    data.nzc[k]         = data.nzc[i]; 
	    data.ampl[k]        = data.ampl[i]; 
	    k++; // move to next index of modified data set 
	 }
      }
   }
   
   if(N!=k){
      sprintf(msg,"[ModifyPlungingProbeData]: Error! N (NMR-ANA) = %d, k matches = %d",N,k);
      MSG = msg;
      rc = Logger::PrintMessage(Logger::kERR,"default",MSG,'a'); 
      return 1;
   }
   data.numTraces = k;

   // print for confirmation
   if(isDiff){
      MSG = "----------- UPDATED PP DATA ------------";
      rc = Logger::PrintMessage(Logger::kINFO,"default",MSG,'a'); 
      for(int i=0;i<k;i++){
	 timeStamp = gm2fieldUtil::GetStringTimeStampFromUTC( data.time[i]/1E+9 ); 
	 sprintf(msg,"run %05d, trace %02d: time = %s, ampl = %.3lf, freq = %.3lf, temp = %.3lf",
	       data.run,data.traceNumber[i],timeStamp.c_str(),data.ampl[i],data.freq[i],data.temp[i]);
         MSG = msg; 
	 rc = Logger::PrintMessage(Logger::kINFO,"default",MSG,'a'); 
      }
   } 
   
   return 0;
}
//______________________________________________________________________________
int LoadImageResults(std::string inpath,std::vector<imageResult_t> &data){

   std::string stype,strial,sheight,simg,simgErr,simgc,simgcErr; 

   imageResult_t res;
   std::ifstream infile;
   infile.open(inpath.c_str());
   if( infile.fail() ){
      std::cout << "Cannot open the file: " << inpath << std::endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,stype   ,',');
         std::getline(infile,strial  ,',');
         std::getline(infile,sheight ,',');
         std::getline(infile,simg    ,',');
         std::getline(infile,simgErr ,',');
         std::getline(infile,simgc   ,',');
         std::getline(infile,simgcErr);
         res.type        = stype; 
         res.trial       = std::atoi( strial.c_str()   ); 
         res.height      = std::atof( sheight.c_str()  ); 
         res.img         = std::atof( simg.c_str()     ); 
         res.img_err     = std::atof( simgErr.c_str()  ); 
         res.img_cor     = std::atof( simgc.c_str()    ); 
         res.img_cor_err = std::atof( simgcErr.c_str() ); 
	 data.push_back(res); 
      }
      infile.close(); 
      data.pop_back(); 
   }
   return 0;
}
//______________________________________________________________________________
int LoadImageParameters(std::string inpath,std::string type,std::vector<imageParameter_t> &data){
   // load relevant parameters for image measurements 

   json input;
   int rc = gm2fieldUtil::Import::ImportJSON(inpath,input);

   imageParameter_t in;
   in.numEvents = (int)input["num-events"];

   std::string index;
   const int N = input[type].size();
   for(int i=0;i<N;i++){
      in.trial     = (int)input[type][i]["trial"];
      in.midasRun  = (int)input[type][i]["midas-run"];
      in.nmrDAQRun = (int)input[type][i]["nmr-daq-run"];
      in.height    = (double)input[type][i]["fxpr-height"];
      // std::cout << Form("trial: %d, MIDAS: %d, NMR-DAQ: %d, height = %.3lf mm",in.trial,in.midasRun,in.nmrDAQRun,in.height) << std::endl;
      data.push_back(in);
   }
   return rc;
}
//______________________________________________________________________________
int LoadScanTimes(int runNumber,int runPeriod,std::string prodVersion,std::vector<double> &time){
   // load in the by-hand determined swap times 
   std::vector<std::string> stime;
   unsigned long int aTime;
   std::string st;
   char inpath[200];
   sprintf(inpath,"./input/scan-times/run-%d/%s/run-%d.txt",runPeriod,prodVersion.c_str(),runNumber);
   std::ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      std::cout << "Cannot open the file: " << inpath << std::endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,st);
         aTime = gm2fieldUtil::GetUTCTimeStampFromString(st,true);
         stime.push_back(st);
         time.push_back(aTime);
      }
      time.pop_back();
      infile.close();
   }

   std::cout << "Time stamps from run " << runNumber << ":" << std::endl;
   const int N = time.size();
   for(int i=0;i<N;i++) std::cout << Form("%s (%.0lf)",stime[i].c_str(),time[i]) << std::endl;

   return 0;
}
//______________________________________________________________________________
int LoadTimes(int probe,int runPeriod,std::string prodVersion,std::string type,std::string dev,std::vector<double> &time){
   // load in the by-hand determined swap times 
   std::vector<std::string> stime;
   unsigned long int aTime;
   std::string st;
   char inpath[200];
   sprintf(inpath,"./input/%s-times/run-%d/%s/%s-%02d.txt",type.c_str(),runPeriod,prodVersion.c_str(),dev.c_str(),probe);
   std::ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      std::cout << "Cannot open the file: " << inpath << std::endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,st);
         aTime = gm2fieldUtil::GetUTCTimeStampFromString(st,true);
         stime.push_back(st);
         time.push_back(aTime);
      }
      time.pop_back();
      infile.close();
   }

   const int N = time.size();
   for(int i=0;i<N;i++) std::cout << Form("%s (%.0lf)",stime[i].c_str(),time[i]) << std::endl;

   return 0;
}
//______________________________________________________________________________
int LoadSCCTimes(int probe,int runPeriod,std::string prodVersion,std::string dev, 
                 std::vector<double> &sccOff,std::vector<double> &sccOn){
   // load in the by-hand determined SCC times 
   std::vector<double> time;
   int rc = LoadTimes(probe,runPeriod,prodVersion,"scc",dev,time); 
   // now fill the vectors appropriately 
   int N = time.size();
   for(int i=0;i<N;i++){
      // odd -- SCC is on 
      if(i%2!=0) sccOn.push_back(time[i]);
      // even -- SCC is off 
      if(i%2==0) sccOff.push_back(time[i]);
   }
   return 0;
}
//______________________________________________________________________________
int LoadIMGTimes(std::string type,int trial,std::vector<double> &time){
   // load in the by-hand determined swap times 
   std::vector<std::string> stime;
   unsigned long int aTime;
   std::string st;
   char inpath[200];
   sprintf(inpath,"./input/img-times/%s_t%02d.txt",type.c_str(),trial);
   std::ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      std::cout << "Cannot open the file: " << inpath << std::endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,st);
         aTime = gm2fieldUtil::GetUTCTimeStampFromString(st,true);
         stime.push_back(st);
         time.push_back(aTime);
      }
      time.pop_back();
      infile.close();
   }

   // std::cout << "[LoadIMGTimes]: Read in the values: " << std:: endl;
   // const int N = time.size();
   // for(int i=0;i<N;i++) std::cout << Form("%s (%.0lf)",stime[i].c_str(),time[i]) << std::endl;

   return 0;
}

//______________________________________________________________________________
int LoadResultsProdFinalData(const char *inpath,result_prod_t &data){

   const int NLines = 1;
   const int SIZE = 2048;  
   char buf[SIZE]; 

   int cntr=0;
   std::string stype,sx,sShot,sMisalign,sFree;

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      for(int i=0;i<NLines;i++) infile.getline(buf,SIZE); 
      while( !infile.eof() ){
         std::getline(infile,stype    ,',');
         std::getline(infile,sx       ,','); 
         std::getline(infile,sShot    ,',');
         std::getline(infile,sMisalign,',');
         std::getline(infile,sFree);
	 cntr++;
	 if(cntr==1){
	    data.diff    = std::atof( sx.c_str() ); 
	    data.diffErr = std::atof( sShot.c_str() );
	    data.mErr    = std::atof( sMisalign.c_str() );
	    data.pErr    = std::atof( sFree.c_str() );
         }else if(cntr==2){
	    data.diff_aba    = std::atof( sx.c_str() ); 
	    data.diffErr_aba = std::atof( sShot.c_str() );
	    data.mErr_aba    = std::atof( sMisalign.c_str() );
	    data.pErr_aba    = std::atof( sFree.c_str() );
         }else if(cntr==3){
	    data.diff_opt    = std::atof( sx.c_str() ); 
	    data.diffErr_opt = std::atof( sShot.c_str() );
	    data.mErr_opt    = std::atof( sMisalign.c_str() );
	    data.pErr_opt    = std::atof( sFree.c_str() );
         } 
      }
      infile.close();
   }

   return 0;
}
//______________________________________________________________________________
int LoadImagesData(const char *inpath,int probe,double &image,double &image_err){

   int ip=0;
   std::string sp,sx,sdx;

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,sp,',');
         std::getline(infile,sx,',');
         std::getline(infile,sdx);
         ip = std::atof( sp.c_str() ); 
         if(probe==ip){
	    image     = std::atof( sx.c_str() ); 
	    image_err = std::atof( sdx.c_str() );
         } 
      }
      infile.close();
   }

   return 0;
}
//______________________________________________________________________________
int LoadResultsProdData(const char *inpath,result_prod_t &data){

   const int NLines = 1;
   const int SIZE = 2048;  
   char buf[SIZE]; 

   int cntr=0;
   std::string sx,sdx;

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      for(int i=0;i<NLines;i++) infile.getline(buf,SIZE); 
      while( !infile.eof() ){
         std::getline(infile,sx,',');
         std::getline(infile,sdx);
	 cntr++;
	 if(cntr==1){
	    data.diff    = std::atof( sx.c_str() ); 
	    data.diffErr = std::atof( sdx.c_str() );
         }else if(cntr==2){
	    data.diff_aba    = std::atof( sx.c_str() ); 
	    data.diffErr_aba = std::atof( sdx.c_str() );
         }else if(cntr==3){
	    data.diff_opt    = std::atof( sx.c_str() ); 
	    data.diffErr_opt = std::atof( sdx.c_str() );
         } 
      }
      infile.close();
   }

   // std::cout << data.diff     << " +/- " << data.diffErr     << std::endl; 
   // std::cout << data.diff_aba << " +/- " << data.diffErr_aba << std::endl; 

   return 0;
}
//______________________________________________________________________________
int LoadMisalignmentData(const char *inpath,misalignment_t &data){

   std::string sname,sdq,sdB_q,sdq_aba,sdB_q_aba,sdq_opt,sdB_q_opt;

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,sname    ,',');
         std::getline(infile,sdq      ,',');
         std::getline(infile,sdB_q    ,',');
         std::getline(infile,sdq_aba  ,',');
         std::getline(infile,sdB_q_aba,',');
         std::getline(infile,sdq_opt  ,',');
         std::getline(infile,sdB_q_opt);
         if( sname.compare("r")==0 ){
	    data.dx       = std::atof( sdq.c_str()       );
	    data.dB_x     = std::atof( sdB_q.c_str()     );
	    data.dx_aba   = std::atof( sdq_aba.c_str()   );
	    data.dB_x_aba = std::atof( sdB_q_aba.c_str() );
	    data.dx_opt   = std::atof( sdq_opt.c_str()   );
	    data.dB_x_opt = std::atof( sdB_q_opt.c_str() );
         }else if( sname.compare("y")==0 ) {
	    data.dy       = std::atof( sdq.c_str()       );
	    data.dB_y     = std::atof( sdB_q.c_str()     );
	    data.dy_aba   = std::atof( sdq_aba.c_str()   );
	    data.dB_y_aba = std::atof( sdB_q_aba.c_str() );
	    data.dy_opt   = std::atof( sdq_opt.c_str()   );
	    data.dB_y_opt = std::atof( sdB_q_opt.c_str() );
	 }else if( sname.compare("z")==0 ){
	    data.dz       = std::atof( sdq.c_str()       );
	    data.dB_z     = std::atof( sdB_q.c_str()     );
	    data.dz_aba   = std::atof( sdq_aba.c_str()   );
	    data.dB_z_aba = std::atof( sdB_q_aba.c_str() );
	    data.dz_opt   = std::atof( sdq_opt.c_str()   );
	    data.dB_z_opt = std::atof( sdB_q_opt.c_str() );
	 }
      }
      infile.close();
   }

   return 0;
}
//______________________________________________________________________________
int LoadMisalignmentCorData(const char *inpath,std::vector<misalignCor_t> &data){

   misalignCor_t dataPt; 
   std::string sname,sv,se;

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,sname,',');
         std::getline(infile,sv   ,',');
         std::getline(infile,se);
	 dataPt.label = sname;
         dataPt.val   = std::atof( sv.c_str() ); 
         dataPt.err   = std::atof( se.c_str() ); 
	 data.push_back(dataPt); 
      }
      data.pop_back();
      infile.close();
   }

   return 0;
}
//______________________________________________________________________________
int LoadCalibSwapData(const char *inpath,std::vector<calibSwap_t> &data){

   int i=0;
   std::string stime,sf,sfe,st,ste,sx,sxe,sy,sye,sz,sze;

   calibSwap_t dataPt;

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,stime,',');
         std::getline(infile,sf   ,',');
         std::getline(infile,sfe  ,',');
         std::getline(infile,st   ,',');
         std::getline(infile,ste  ,',');
         std::getline(infile,sx   ,',');
         std::getline(infile,sxe  ,',');
         std::getline(infile,sy   ,',');
         std::getline(infile,sye  ,',');
         std::getline(infile,sz   ,',');
         std::getline(infile,sze);
	 dataPt.time    = std::atof( stime.c_str() );
	 dataPt.freq    = std::atof( sf.c_str()    );
	 dataPt.freqErr = std::atof( sfe.c_str()   );
	 dataPt.temp    = std::atof( st.c_str()    );
	 dataPt.tempErr = std::atof( ste.c_str()   );
	 dataPt.r       = std::atof( sx.c_str()    );
	 dataPt.rErr    = std::atof( sxe.c_str()   );
	 dataPt.y       = std::atof( sy.c_str()    );
	 dataPt.yErr    = std::atof( sye.c_str()   );
	 dataPt.phi     = std::atof( sz.c_str()    );
	 dataPt.phiErr  = std::atof( sze.c_str()   );
	 data.push_back(dataPt); 
      }
      infile.close();
      data.pop_back();
   }

   return 0;
}
//______________________________________________________________________________
int LoadImposedGradientData(const char *inpath,imposed_gradient_t &data){

   int i=0;
   std::string sp,sgrad,sgrad_err;

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,sp,',');
         std::getline(infile,sgrad,',');
         std::getline(infile,sgrad_err);
	 data.pos      = std::atof( sp.c_str()        );
	 data.grad     = std::atof( sgrad.c_str()     );
	 data.grad_err = std::atof( sgrad_err.c_str() );
      }
      infile.close();
   }

   return 0;
}
//______________________________________________________________________________
int LoadImposedGradientData(const char *inpath,std::vector<imposed_gradient_t> &data){

   int i=0;
   std::string sp,sgrad,sgrad_err;

   imposed_gradient_t dataPt;

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,sp,',');
         std::getline(infile,sgrad,',');
         std::getline(infile,sgrad_err);
	 dataPt.pos      = std::atof( sp.c_str()        );
	 dataPt.grad     = std::atof( sgrad.c_str()     );
	 dataPt.grad_err = std::atof( sgrad_err.c_str() );
	 data.push_back(dataPt); 
      }
      infile.close();
      data.pop_back();
   }

   return 0;
}
//______________________________________________________________________________
int LoadImposedAziGradData_bak(const char *inpath,int probe,double &dBdz){

   int ipr=0;
   std::string stp,sbg,sig_82,sg,sg_per_A,ssc;  // bare grad, grad at 0.82 A, grad, grad/Amp,shim current (A) 

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,stp     ,',');
         std::getline(infile,sbg     ,',');
         std::getline(infile,sig_82  ,',');
         std::getline(infile,sg      ,',');
         std::getline(infile,sg_per_A,',');
         std::getline(infile,ssc);
         ipr = std::atoi( stp.c_str() ); 
         if(ipr==probe) dBdz = std::atof( sg.c_str() ); 
      }
      infile.close();
   }

   return 0;
}
//______________________________________________________________________________
int LoadImposedAziGradData_old(const char *inpath,int probe,double &dBdz,double &dBdz_err){

   int ipr=0;
   std::string stp,sgr,sgr_err;  // probe, imposed gradient @ 0.82 A (Hz/mm)  

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,stp,',');
         std::getline(infile,sgr,',');
         std::getline(infile,sgr_err);
         ipr = std::atoi( stp.c_str() ); 
         if(ipr==probe){
	    dBdz     = std::atof( sgr.c_str() ); 
            dBdz_err = std::atof( sgr_err.c_str() );
         } 
      }
      infile.close();
   }

   return 0;
}
//______________________________________________________________________________
int LoadImposedAziGradData(const char *inpath,double &dBdz,double &dBdz_err){

   const int NLines = 1;
   const int SIZE = 2048; 
   char buf[SIZE]; 

   int cntr=0;
   std::string sgr,sgr_err;  // imposed gradient @ 0.82 A (Hz/mm)  

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      for(int i=0;i<NLines;i++) infile.getline(buf,SIZE); 
      while( !infile.eof() ){
         std::getline(infile,sgr,',');
	 std::getline(infile,sgr_err);
	 cntr++;
	 if(cntr==1){
	    dBdz     = std::atof( sgr.c_str() ); 
	    dBdz_err = std::atof( sgr_err.c_str() );
	 } 
      }
      infile.close();
   }

   std::cout << "[LoadImposedAziGradData]: dBdz = " << dBdz << " dBdz_err = " << dBdz_err << std::endl;

   return 0;
}
//______________________________________________________________________________
int LoadFitPars(const char *inpath,std::vector<double> &x,std::vector<double> &dx){

   double ix,idx;
   const int NLines = 1;
   const int SIZE = 2048; 
   char buf[SIZE];  
   std::string sx,sdx;    

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      for(int i=0;i<NLines;i++) infile.getline(buf,SIZE); 
      while( !infile.eof() ){
         std::getline(infile,sx,',');
         std::getline(infile,sdx);
	 ix  = std::atof( sx.c_str() ); 
	 idx = std::atof( sdx.c_str() ); 
	 x.push_back(ix); 
	 dx.push_back(idx); 
      }
      x.pop_back(); 
      dx.pop_back(); 
      infile.close();
   }
   return 0;
}
//______________________________________________________________________________
int LoadTrolleyDeltaBData(const char *inpath,trolleyDeltaB_t &data){

   int i=0;
   std::string sp,snq,snqE,ssq,ssqE;

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,sp  ,',');
         std::getline(infile,snq ,',');
         std::getline(infile,snqE,',');
         std::getline(infile,ssq ,',');
         std::getline(infile,ssqE);
         data.probeID[i]      = i+1; 
	 data.normQuad[i]     = std::atof( snq.c_str()  );
	 data.normQuad_err[i] = std::atof( snqE.c_str() );
	 data.skewQuad[i]     = std::atof( ssq.c_str()  );
	 data.skewQuad_err[i] = std::atof( ssqE.c_str() );
	 i++;
      }
      infile.close();
   }

   return 0;
}
//______________________________________________________________________________
int LoadTrolleyPositionData(const char *inpath,trolleyProbePosition_t &data){

   int i=0;
   std::string sp,sr,sdr,sy,sdy,sphi,sdphi;

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,sp,',');
         std::getline(infile,sr,',');
         std::getline(infile,sdr,',');
         std::getline(infile,sy,',');
         std::getline(infile,sdy,',');
         std::getline(infile,sphi,',');
         std::getline(infile,sdphi);
	 data.r[i]    = std::atof( sr.c_str()    );
	 data.dr[i]   = std::atof( sdr.c_str()   );
	 data.y[i]    = std::atof( sy.c_str()    );
	 data.dy[i]   = std::atof( sdy.c_str()   );
	 data.phi[i]  = std::atof( sphi.c_str()  );
	 data.dphi[i] = std::atof( sdphi.c_str() );
	 i++;
      }
      infile.close();
   }

   return 0;
}
//______________________________________________________________________________
int ImportNMRANAData(const char *inpath,std::vector<nmrAnaEvent_t> &Data,std::string cutFile){
   // load cuts 
   Cut *myCuts = new Cut(cutFile); 

   // load data from the NMR-ANA framework 
   nmrAnaEvent_t inData;

   // allowed amplitude in volts 
   // double VMIN = 0.3; 
   // double VMAX = 1.5; 

   int cntr=0;
   int N=0,k=0;
   int irun,ipulse,izc,NUM_PULSES=0;
   double inoise,inc,ifa,ifb,ifc,ifa_ph,ifb_ph,ifc_ph,iampl;

   // const int MAX = 1000;
   // char buf[MAX];

   ifstream infile;
   infile.open(inpath,ios::in);

   bool goodEvent = false;

   if(infile.fail()){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      // cout << "Opening the file: " << inpath << endl;
      // for(int i=0;i<1;i++) infile.getline(buf,MAX);
      while( !infile.eof() ){
         infile >> irun >> ipulse >> iampl >> inoise >> izc >> inc >> ifa >> ifb >> ifc >> ifa_ph >> ifb_ph >> ifc_ph;
         goodEvent = myCuts->CheckEvent_nmrAna(irun,ipulse,izc,iampl,ifc_ph); 
         if(goodEvent){
         // if( izc>10&&ifc_ph>0 && (iampl>VMIN&&iampl<VMAX) ){
	    inData.run     = irun;
	    inData.ampl    = iampl;
	    inData.noise   = inoise;
	    inData.pulse   = ipulse; 
	    inData.zc      = izc; 
	    inData.nc      = inc; 
	    inData.freq[0] = ifa;
	    inData.freq[1] = ifb;
	    inData.freq[2] = ifc;
	    inData.freq[3] = ifa_ph;
	    inData.freq[4] = ifb_ph;
	    inData.freq[5] = ifc_ph;
	    Data.push_back(inData);
	    cntr++;
         }else{
	    // skip the event 
	    // std::cout << Form("[ImportData]: WARNING: rejected run %d pulse %d",irun,ipulse) << std::endl;
	    // std::cout << Form("              ampl = %.3lf V, zc = %d, freq = %.3lf",iampl,izc,ifc_ph) << std::endl;
         }
      }
      infile.close();
      if(cntr>0) Data.pop_back();
   }
   
   delete myCuts; 
   
   return 0;
}
//______________________________________________________________________________
int ImportNMRANAData_new(int run,std::string version,std::vector<nmrAnaEvent_t> &data,std::string cutFile){
  // a new way to load NMR-ANA results
  // need the temperature sensor 
  // ppRunSummary_t myRunSummary;
  gm2fieldUtil::Temperature::Sensor *tempSensor = new gm2fieldUtil::Temperature::Sensor("PT1000");
  tempSensor->SetDataPath(TEMP_DIR);
  // cut data 
  Cut *myCut = new Cut(cutFile);

  nmrAnaEvent_t inData;

  bool goodEvent = false;
  int tr=0,zc=0;
  double ampl=0,freq=0,res=0;
  std::string sp,sch,stime,szc,snc,sampl,snoise,st2,stemp,slo,srf;
  std::string sfmid,sflin,sflsq,sfmid_ph,sflin_ph,sflsq_ph;

  char inpath[512];
  sprintf(inpath,"./input/NMR-ANA/%s/run-%05d/results.csv",version.c_str(),run);
  std::ifstream infile;
  infile.open(inpath);
  if( infile.fail() ){
    std::cout << "Cannot open the file: " << inpath << std::endl;
    return 1;
  }else{
    while( !infile.eof() ){
      std::getline(infile,sp      ,',');
      std::getline(infile,sch     ,',');
      std::getline(infile,stime   ,',');
      std::getline(infile,szc     ,',');
      std::getline(infile,snc     ,',');
      std::getline(infile,sampl   ,',');
      std::getline(infile,snoise  ,',');
      std::getline(infile,st2     ,',');
      std::getline(infile,stemp   ,',');
      std::getline(infile,slo     ,',');
      std::getline(infile,srf     ,',');
      std::getline(infile,sfmid   ,',');
      std::getline(infile,sflin   ,',');
      std::getline(infile,sflsq   ,',');
      std::getline(infile,sfmid_ph,',');
      std::getline(infile,sflin_ph,',');
      std::getline(infile,sflsq_ph);
      // check the event
      try{
        tr   = std::atoi( sp.c_str()  );
        zc   = std::atoi( szc.c_str() );
        ampl = std::atof( sampl.c_str() );
        freq = std::atof( sflsq_ph.c_str() );
        res  = std::atof( stemp.c_str() );
      }catch(const std::invalid_argument &ia){
        std::cout << "[LoadNMRANAData]: Error: " << ia.what() << ".  Maybe EOF?" << std::endl;
      }
      goodEvent = myCut->CheckEvent_nmrAna(run,tr,zc,ampl,freq);
      if(goodEvent){
        inData.time     = std::stoull( stime.c_str() );
        inData.temp     = tempSensor->GetTemperature(res);
        inData.freq_LO  = std::atof( slo.c_str() );
        inData.freq_pi2 = std::atof( srf.c_str() );
        inData.ampl     = ampl;
        inData.noise    = std::atof( snoise.c_str() );
        inData.t2       = std::atof( st2.c_str() );
        inData.freq[0]  = std::atof( sfmid.c_str() );
        inData.freq[1]  = std::atof( sflin.c_str() );
        inData.freq[2]  = std::atof( sflsq.c_str() );
        inData.freq[3]  = std::atof( sfmid_ph.c_str() );
        inData.freq[4]  = std::atof( sflin_ph.c_str() );
        inData.freq[5]  = std::atof( sflsq_ph.c_str() );
        inData.nc       = std::atof( snc.c_str() );
        inData.pulse    = tr;
        inData.run      = run;
        data.push_back(inData);
      }
    }
    infile.close();
    data.pop_back();
  }
  delete tempSensor;
  return 0;
}
//______________________________________________________________________________
int LoadRunSummaryData(const char *inpath,runSummary_t &x){

   std::string tag,value;
   std::string runNo        = "run_number";
   std::string date         = "date";
   std::string startTime    = "start_time";
   std::string endTime      = "end_time";
   std::string numPulses    = "num_pulses";
   std::string numSamples   = "num_samples_per_pulse";
   std::string adcID        = "adc_id";
   std::string adcCh        = "adc_ch_number";
   std::string clockFreq    = "adc_clock_frequency";
   std::string expectedFreq = "expected_IF_frequency";
   std::string loFreq       = "LO_frequency";
   std::string pi2Freq      = "pi2_frequency_1";
   std::string pi2Power     = "pi2_power_1";
   std::string pi2Voltage   = "pi2_voltage_1";
   std::string bncVoltage   = "bnc_voltage";
   std::string nTypeVoltage = "ntype_voltage";
   std::string tempSensor   = "temp_sensor";
   std::ifstream infile;

   infile.open(inpath);
   if( infile.fail() ){
      std::cout << "Cannot open the file: " << inpath << std::endl;
      return 1;
   }else{
      while( !infile.eof() ){
         infile >> tag >> value;
         if( tag.compare(runNo)==0 )        x.runNumber    = std::atoi( value.c_str() );
         if( tag.compare(date)==0 )         x.date         = value;
         if( tag.compare(startTime)==0 )    x.startTime    = value;
         if( tag.compare(endTime)==0 )      x.endTime      = value;
         if( tag.compare(numPulses)==0 )    x.numPulses    = std::atoi( value.c_str() );
         if( tag.compare(numSamples)==0 )   x.numSamples   = std::atoi( value.c_str() );
         if( tag.compare(adcID)==0 )        x.adcID        = std::atoi( value.c_str() );
         if( tag.compare(adcCh)==0 )        x.adcChannel   = std::atoi( value.c_str() );
         if( tag.compare(clockFreq)==0 )    x.adcClockFreq = std::atof( value.c_str() );
         if( tag.compare(expectedFreq)==0 ) x.expFreq      = std::atof( value.c_str() );
         if( tag.compare(loFreq)==0 )       x.loFreq       = std::atof( value.c_str() );
         if( tag.compare(bncVoltage)==0 )   x.bncVoltage   = std::atof( value.c_str() );
         if( tag.compare(nTypeVoltage)==0 ) x.nTypeVoltage = std::atof( value.c_str() );
         if( tag.compare(pi2Freq)==0 )      x.pi2Freq      = std::atof( value.c_str() );
         if( tag.compare(pi2Power)==0 )     x.pi2Power     = std::atof( value.c_str() );
         if( tag.compare(pi2Voltage)==0 )   x.pi2Voltage   = std::atof( value.c_str() );
         if( tag.compare(tempSensor)==0 )   x.tempSensor   = value;
      }
      infile.close();
   }
   return 0;
}

//______________________________________________________________________________
int LoadNMRDAQEventData(const char *inpath,std::vector<NMRDAQEvent_t> &event){
   NMRDAQEvent_t data;  
   std::string sPulseNum,sChNum,sTime,sTemp,sx,sy,sz;

   std::ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      std::cout << "Cannot open the file: " << inpath << std::endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,sPulseNum,',');
         std::getline(infile,sChNum   ,',');
         std::getline(infile,sTime    ,',');
         std::getline(infile,sTemp    ,',');
         std::getline(infile,sx       ,',');
         std::getline(infile,sy       ,',');
         std::getline(infile,sz);
	 data.pulseNum    = std::atoi( sPulseNum.c_str() );
	 data.chNum       = std::atoi( sChNum.c_str() );
	 data.timestamp   = std::stoull(sTime.c_str());
	 data.temperature = std::atof( sTemp.c_str() );
	 data.x           = std::atof( sx.c_str() );
	 data.y           = std::atof( sy.c_str() );
	 data.z           = std::atof( sz.c_str() );
	 event.push_back(data); 
      }
      event.pop_back();
   }

   return 0;
}
//______________________________________________________________________________
int ImportPPEncoderCalib_csv(const char *inpath,std::vector<double> &y){

   int N=0,ix1=0,ix2=0,ix3=0;
   string sx1,sx2,sx3;

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,sx1,',');
         std::getline(infile,sx2,',');
         std::getline(infile,sx3);
         ix1 = std::atof( sx1.c_str() );
         ix2 = std::atof( sx2.c_str() );
         ix3 = std::atof( sx3.c_str() );
         y.push_back(ix1);
         y.push_back(ix2);
         y.push_back(ix3);
      }
      infile.close();
      N = y.size();
      if(N>3) y.pop_back();
   }

   return 0;
}
//______________________________________________________________________________
int ImportDeltaBFileList_csv(const char *inpath,
                    std::vector<int> &x1,std::vector<std::string> &x2,
                    std::vector<double> &x3){

   int ix1=0;
   double ix3;
   string sx1,sx2,sx3;

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,sx1,',');
         std::getline(infile,sx2,',');
         std::getline(infile,sx3);
         ix1 = atoi( sx1.c_str() );
         ix3 = atof( sx3.c_str() );
         x1.push_back(ix1);
         x2.push_back(sx2);
         x3.push_back(ix3);
      }
      infile.close();
      x1.pop_back();
      x2.pop_back();
      x3.pop_back();
   }

   return 0;
}
//______________________________________________________________________________
int LoadPerturbationData_json(const char *inpath,perturbation_t &pert){
   json obj; 
   std::string inpath_str = inpath; 
   int rc = gm2fieldUtil::Import::ImportJSON(inpath_str,obj);
   if(rc!=0){
      std::cout << "Cannot open the file: " << inpath_str << std::endl; 
      return 1;
   }else{
      pert.sigma         = obj["sigma"];  
      pert.sigma_err     = obj["sigma-err"];  
      pert.chi           = obj["chi"];  
      pert.chi_err       = obj["chi-err"];  
      pert.eps           = obj["eps"];  
      pert.eps_err       = obj["eps-err"];  
      pert.delta_m       = obj["delta-mat"];  
      pert.delta_m_err   = obj["delta-mat-err"];  
      pert.delta_eps     = obj["delta-eps"];  
      pert.delta_eps_err = obj["delta-eps-err"];  
      pert.delta_mag     = obj["delta-mag"];  
      pert.delta_mag_err = obj["delta-mag-err"];
      pert.delta_t       = obj["delta-t"];  
      pert.delta_t_err   = obj["delta-t-err"];  
   } 
   return 0;
}
//______________________________________________________________________________
int LoadPerturbationData(const char *inpath,int probe,perturbation_t &pert){

   int cntr=0;
   string s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12;

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,s1,',');
         std::getline(infile,s2,',');
         std::getline(infile,s3,',');
         std::getline(infile,s4,',');
         std::getline(infile,s5,',');
         std::getline(infile,s6,',');
         std::getline(infile,s7,',');
         std::getline(infile,s8,',');
         std::getline(infile,s9,',');
         std::getline(infile,s10,',');
         std::getline(infile,s11,',');
         std::getline(infile,s12);
	 cntr++; 
	 if(cntr==1){
	    pert.sigma         = std::atof( s1.c_str() );
	    pert.sigma_err     = std::atof( s2.c_str() );
	    pert.chi           = std::atof( s3.c_str() );
	    pert.chi_err       = std::atof( s4.c_str() );
	    pert.eps           = std::atof( s5.c_str() );
	    pert.eps_err       = std::atof( s6.c_str() );
	    pert.delta_m       = std::atof( s7.c_str() );
	    pert.delta_m_err   = std::atof( s8.c_str() );
	    pert.delta_eps     = std::atof( s9.c_str() );
	    pert.delta_eps_err = std::atof( s10.c_str() );
	    pert.delta_mag     = std::atof( s11.c_str() );
	    pert.delta_mag_err = std::atof( s12.c_str() );
	 }
      }
      infile.close();
   }
 
   // // get the correct value of the image 
   // double image=0,image_err=0;
   // char inpath_image[200];
   // sprintf(inpath_image,"./input/perturbation/images.csv"); 
   // LoadImagesData(inpath_image,probe,image,image_err);
   // pert.delta_mag     = image; 
   // pert.delta_mag_err = image_err;   
 
   return 0;
}

//______________________________________________________________________________
int LoadFieldData(const char *inpath,nmr_meas_t &data){

   int cntr=0;
   std::string sL,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12;

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,sL,',');
         std::getline(infile,s1,',');
         std::getline(infile,s2,',');
         std::getline(infile,s3,',');
         std::getline(infile,s4,',');
         std::getline(infile,s5,',');
         std::getline(infile,s6,',');
         std::getline(infile,s7,',');
	 std::getline(infile,s8,',');
	 std::getline(infile,s9,',');
	 std::getline(infile,s10,',');
	 std::getline(infile,s11,',');
	 std::getline(infile,s12);
	 cntr++;
	 if(cntr==1){
	    data.name          = sL;
	    data.freq          = std::atof( s1.c_str() );
	    data.freq_err      = std::atof( s2.c_str() );
	    data.freq_fxpr     = std::atof( s3.c_str() );
	    data.freq_fxpr_err = std::atof( s4.c_str() );
	    data.freq_trly     = std::atof( s5.c_str() );
	    data.freq_trly_err = std::atof( s6.c_str() );
            data.p2p_err       = std::atof( s8.c_str() ); 
            data.r2r_err       = 0.;  // no run-to-run uncertainty here  
	 }
      }
      infile.close();
   }

   std::cout << data.name << " " 
             << Form("raw     = %.3lf +/- %.3lf",data.freq,data.freq_err)           << " " 
             << Form("fxpr    = %.3lf +/- %.3lf",data.freq_fxpr,data.freq_fxpr_err) << " "  
             << Form("trly    = %.3lf +/- %.3lf",data.freq_trly,data.freq_trly_err) << " " 
             << Form("p2p_err = %.3lf",data.p2p_err)                                << std::endl; 

   return 0;
}
//______________________________________________________________________________
int LoadGradientData(const char *inpath,grad_meas_t &x){

   grad_meas_t data;

   std::string sL,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12;

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,sL,',');
         std::getline(infile,s1,',');
         std::getline(infile,s2,',');
         std::getline(infile,s3,',');
         std::getline(infile,s4,',');
         std::getline(infile,s5,',');
         std::getline(infile,s6,',');
         std::getline(infile,s7,',');
         std::getline(infile,s8,',');
         std::getline(infile,s9,',');
         std::getline(infile,s10,',');
         std::getline(infile,s11,',');
         std::getline(infile,s12);
	 x.name           = sL;
         x.grad           = std::atof( s1.c_str() );
         x.grad_err       = std::atof( s2.c_str() );
         x.grad_fxpr      = std::atof( s3.c_str() );
         x.grad_fxpr_err  = std::atof( s4.c_str() );
         x.grad_trly      = std::atof( s5.c_str() );
	 x.grad_trly_err  = std::atof( s6.c_str() );
         x.drift_fxpr     = std::atof( s9.c_str() );
         x.drift_fxpr_err = std::atof( s10.c_str() );
         x.drift_trly     = std::atof( s11.c_str() );
	 x.drift_trly_err = std::atof( s12.c_str() );
      }
      infile.close();
   }

   // const int N = x.size();
   // for(int i=0;i<N;i++){
   //    std::cout << x[i].name << " " 
   //              << x[i].grad      << " +/- " << x[i].grad_err      << " " 
   //              << x[i].grad_fxpr << " +/- " << x[i].grad_fxpr_err << " " 
   //              << x[i].grad_trly << " +/- " << x[i].grad_trly_err << std::endl;
   // }
   // std::cout << "--------------" << std::endl; 

   return 0;
}
//______________________________________________________________________________
int LoadGradientData(const char *inpath,std::vector<grad_meas_t> &x){

   grad_meas_t data;

   std::string sL,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12;

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,sL,',');
         std::getline(infile,s1,',');
         std::getline(infile,s2,',');
         std::getline(infile,s3,',');
         std::getline(infile,s4,',');
         std::getline(infile,s5,',');
         std::getline(infile,s6,',');
         std::getline(infile,s7,',');
         std::getline(infile,s8,',');
         std::getline(infile,s9,',');
         std::getline(infile,s10,',');
         std::getline(infile,s11,',');
         std::getline(infile,s12);
	 data.name           = sL;
         data.grad           = std::atof( s1.c_str() );
         data.grad_err       = std::atof( s2.c_str() );
         data.grad_fxpr      = std::atof( s3.c_str() );
         data.grad_fxpr_err  = std::atof( s4.c_str() );
         data.grad_trly      = std::atof( s5.c_str() );
	 data.grad_trly_err  = std::atof( s6.c_str() );
         data.drift_fxpr     = std::atof( s9.c_str() );
         data.drift_fxpr_err = std::atof( s10.c_str() );
         data.drift_trly     = std::atof( s11.c_str() );
	 data.drift_trly_err = std::atof( s12.c_str() );
	 x.push_back(data);
      }
      infile.close();
      x.pop_back();
   }

   // const int N = x.size();
   // for(int i=0;i<N;i++){
   //    std::cout << x[i].name << " " 
   //              << x[i].grad      << " +/- " << x[i].grad_err      << " " 
   //              << x[i].grad_fxpr << " +/- " << x[i].grad_fxpr_err << " " 
   //              << x[i].grad_trly << " +/- " << x[i].grad_trly_err << std::endl;
   // }
   // std::cout << "--------------" << std::endl; 

   return 0;
}
//______________________________________________________________________________
int LoadDeltaBData_trlyXYZ(const char *inpath,int probe,std::vector<deltab_t> &x){

   int ipr=0;
   deltab_t data;
   std::string stp,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12;

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,stp,',');
         std::getline(infile,s1 ,',');
         std::getline(infile,s2 ,',');
         std::getline(infile,s3 ,',');
         std::getline(infile,s4 ,',');
         std::getline(infile,s5 ,',');
         std::getline(infile,s6 ,',');
         std::getline(infile,s7 ,',');
         std::getline(infile,s8 ,',');
         std::getline(infile,s9 ,',');
         std::getline(infile,s10,',');
         std::getline(infile,s11,',');
         std::getline(infile,s12);
         ipr = std::atoi( stp.c_str() );
	 // std::cout << "PROBE " << ipr << std::endl;
         if(ipr==probe){
	    // std::cout << "--> MATCH" << std::endl;
	    // radial 
	    data.name = "rad-grad";
	    data.dB             = std::atof( s1.c_str()  );
	    data.dB_err         = std::atof( s2.c_str()  );
	    data.dB_fxpr        = std::atof( s3.c_str()  );
	    data.dB_fxpr_err    = std::atof( s4.c_str()  );
	    data.dB_trly        = 0;
	    data.dB_trly_err    = 0;
	    x.push_back(data);
	    // vertical 
	    data.name = "vert-grad";
	    data.dB             = std::atof( s5.c_str()  );
	    data.dB_err         = std::atof( s6.c_str()  );
	    data.dB_fxpr        = std::atof( s7.c_str()  );
	    data.dB_fxpr_err    = std::atof( s8.c_str()  );
	    data.dB_trly        = 0;
	    data.dB_trly_err    = 0;
	    x.push_back(data);
	    // azimuthal 
	    data.name = "azi";
	    data.dB             = std::atof( s9.c_str()  );
	    data.dB_err         = std::atof( s10.c_str()  );
	    data.dB_fxpr        = std::atof( s11.c_str()  );
	    data.dB_fxpr_err    = std::atof( s12.c_str()  );
	    data.dB_trly        = 0;
	    data.dB_trly_err    = 0;
	    x.push_back(data);
	 }
      }
      infile.close();
      // x.pop_back();
   }

   // const int N = x.size();
   // std::cout << N << " entries found" << std::endl;
   // for(int i=0;i<N;i++){
   //    std::cout << x[i].name << " " 
   //              << x[i].dB         << " +/- " << x[i].dB_err         << " " 
   //              << x[i].dB_fxpr    << " +/- " << x[i].dB_fxpr_err    << " " 
   //              << x[i].dB_trly    << " +/- " << x[i].dB_trly_err    << " " 
   //              << x[i].drift_fxpr << " +/- " << x[i].drift_fxpr_err << " " 
   //              << x[i].drift_trly << " +/- " << x[i].drift_trly_err << std::endl;
   // }
   // std::cout << "--------------" << std::endl; 

   return 0;
}
//______________________________________________________________________________
int LoadDeltaBData_trlyXY(const char *inpath,int probe,std::vector<deltab_t> &x){

   int ipr=0;
   deltab_t data;
   std::string stp,s1,s2,s3,s4,s5,s6,s7,s8,s9;

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,stp,',');
         std::getline(infile,s1 ,',');
         std::getline(infile,s2 ,',');
         std::getline(infile,s3 ,',');
         std::getline(infile,s4 ,',');
         std::getline(infile,s5 ,',');
         std::getline(infile,s6 ,',');
         std::getline(infile,s7 ,',');
         std::getline(infile,s8);
         ipr = std::atoi( stp.c_str() );
	 // std::cout << "PROBE " << ipr << std::endl;
         if(ipr==probe){
	    // std::cout << "--> MATCH" << std::endl;
	    // radial 
	    data.name = "rad-grad";
	    data.dB             = std::atof( s1.c_str()  );
	    data.dB_err         = std::atof( s2.c_str()  );
	    data.dB_fxpr        = std::atof( s3.c_str()  );
	    data.dB_fxpr_err    = std::atof( s4.c_str()  );
	    data.dB_trly        = 0;
	    data.dB_trly_err    = 0;
	    x.push_back(data);
	    // vertical 
	    data.name = "vert-grad";
	    data.dB             = std::atof( s5.c_str()  );
	    data.dB_err         = std::atof( s6.c_str()  );
	    data.dB_fxpr        = std::atof( s7.c_str()  );
	    data.dB_fxpr_err    = std::atof( s8.c_str()  );
	    data.dB_trly        = 0;
	    data.dB_trly_err    = 0;
	    x.push_back(data);
	 }
      }
      infile.close();
      // x.pop_back();
   }

   // const int N = x.size();
   // std::cout << N << " entries found" << std::endl;
   // for(int i=0;i<N;i++){
   //    std::cout << x[i].name << " " 
   //              << x[i].dB         << " +/- " << x[i].dB_err         << " " 
   //              << x[i].dB_fxpr    << " +/- " << x[i].dB_fxpr_err    << " " 
   //              << x[i].dB_trly    << " +/- " << x[i].dB_trly_err    << " " 
   //              << x[i].drift_fxpr << " +/- " << x[i].drift_fxpr_err << " " 
   //              << x[i].drift_trly << " +/- " << x[i].drift_trly_err << std::endl;
   // }
   // std::cout << "--------------" << std::endl; 

   return 0;
}
//______________________________________________________________________________
int LoadDeltaBData(const char *inpath,std::vector<deltab_t> &x){

   deltab_t data;
   std::string sL,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12;

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,sL,',');
         std::getline(infile,s1,',');
         std::getline(infile,s2,',');
         std::getline(infile,s3,',');
         std::getline(infile,s4,',');
         std::getline(infile,s5,',');
         std::getline(infile,s6,','); 
         std::getline(infile,s7,',');
         std::getline(infile,s8,',');
         std::getline(infile,s9,',');
         std::getline(infile,s10,',');
         std::getline(infile,s11,',');
         std::getline(infile,s12);
	 data.name           = sL;
         data.dB             = std::atof( s1.c_str()  );
         data.dB_err         = std::atof( s2.c_str()  );
         data.dB_fxpr        = std::atof( s3.c_str()  );
         data.dB_fxpr_err    = std::atof( s4.c_str()  );
         data.dB_trly        = std::atof( s5.c_str()  );
	 data.dB_trly_err    = std::atof( s6.c_str()  );
         data.drift_fxpr     = std::atof( s9.c_str()  );   // WARNING: note that we're starting at the 9th entry! 
         data.drift_fxpr_err = std::atof( s10.c_str() );
         data.drift_trly     = std::atof( s11.c_str() );
	 data.drift_trly_err = std::atof( s12.c_str() );
	 x.push_back(data);
      }
      infile.close();
      x.pop_back();
   }

   // const int N = x.size();
   // for(int i=0;i<N;i++){
   //    std::cout << x[i].name << " " 
   //              << x[i].dB         << " +/- " << x[i].dB_err         << " " 
   //              << x[i].dB_fxpr    << " +/- " << x[i].dB_fxpr_err    << " " 
   //              << x[i].dB_trly    << " +/- " << x[i].dB_trly_err    << " " 
   //              << x[i].drift_fxpr << " +/- " << x[i].drift_fxpr_err << " " 
   //              << x[i].drift_trly << " +/- " << x[i].drift_trly_err << std::endl;
   // }
   // std::cout << "--------------" << std::endl; 

   return 0;
}
//______________________________________________________________________________
int FindStartIndexTRLY(std::string date,int runNumber){
   // find the true start index to keep events for the trolley data because of the 
   // issue with interactive vs continuous mode
   int index=22,theRun=0; 
   std::string sr,sL;

   char inpath[500]; 
   sprintf(inpath,"./input/runlists/%s/trly-mode.csv",date.c_str() );

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,sr,',');
         std::getline(infile,sL);
         theRun = std::atoi(sr.c_str()); 
         if(runNumber==theRun){
	    std::cout << "Found run " << runNumber << std::endl;
	    if( sL.compare("interactive")==0 ) index = 0;
	    if( sL.compare("continuous")==0  ) index = 22;
         }
      }
      infile.close();
   }

   std::cout << "--> For run " << runNumber << ", the start index is " << index << std::endl; 

   return index;
  
}
//______________________________________________________________________________
int ImportBlinding(blind_t &data){
   // load the blinding data 
   int cntr=0;
   std::string sv,sv2,serr_x,serr_y,serr_z,serr_f1,serr_f2;

   std::string inpath = "./misc/blind/values.csv"; 

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,sv,',');
         std::getline(infile,sv2,',');
         std::getline(infile,serr_x,',');
         std::getline(infile,serr_y,',');
         std::getline(infile,serr_z,',');
         std::getline(infile,serr_f1,',');
         std::getline(infile,serr_f2);
         cntr++;
	 if(cntr==1){ // for some reason we need this...
	    data.value_pp = std::atof( sv.c_str() ); 
	    data.value_tr = std::atof( sv2.c_str() ); 
	    data.err_x    = std::atof( serr_x.c_str() ); 
	    data.err_y    = std::atof( serr_y.c_str() ); 
	    data.err_z    = std::atof( serr_z.c_str() ); 
	    data.err_pp   = std::atof( serr_f1.c_str() ); 
	    data.err_tr   = std::atof( serr_f2.c_str() );
	 } 
      }
      infile.close();
   }

   return 0;
}
//______________________________________________________________________________
int SortRuns(std::vector<std::string> label,std::vector<int> allRuns,
             std::vector<int> &run,std::vector<int> &driftRun,std::vector<int> &index){

   const int NRUN = allRuns.size(); 
   int bareRun=0,bareIndex=0;
   int gradRun=0,gradIndex=0;
   int bare2Run=0,bare2Index=0;
   for(int i=0;i<NRUN;i++){
      if(label[i].compare("bare")==0){
         bareRun   = allRuns[i];
         bareIndex = i;
         run.push_back(allRuns[i]);
      }else if(label[i].compare("bare-2")==0){
         bare2Run   = allRuns[i];
         bare2Index = i;
      }else{
         gradRun   = allRuns[i];
         gradIndex = i;
         run.push_back(allRuns[i]);
      }
   }

   // drift runs 
   driftRun.push_back( allRuns[bareIndex]  ); 
   driftRun.push_back( allRuns[bare2Index] );

   // save indices 
   index.push_back( bareIndex  );  
   index.push_back( gradIndex  );  
   index.push_back( bare2Index );  

   return 0;

}
//______________________________________________________________________________
int SortRunsAlt(std::vector<std::string> label,std::vector<int> allRuns,
             std::vector<int> &run,std::vector<int> &index){

   const int NRUN = allRuns.size(); 
   int bareRun=0,bareIndex=0;
   int gradRun=0,gradIndex=0;
   int bare2Run=0,bare2Index=0;
   for(int i=0;i<NRUN;i++){
      if(label[i].compare("bare")==0){
         bareIndex = i;
      }else if(label[i].compare("bare-2")==0){
         bare2Index = i;
      }else{
         gradIndex = i;
      }
      run.push_back(allRuns[i]);
   }

   // save indices 
   index.push_back( bareIndex  );  
   index.push_back( gradIndex  );  
   index.push_back( bare2Index );  

   return 0;

}
//______________________________________________________________________________
int ImportResults(std::string inpath,result_t &data){
   // load the results data 
   int cntr=0;
   std::string sl,s1,s1err,s2,s2err,s3,s3err;

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,sl,',');
         std::getline(infile,s1,',');
         std::getline(infile,s1err,',');
         std::getline(infile,s2,',');
         std::getline(infile,s2err,',');
         std::getline(infile,s3,',');
         std::getline(infile,s3err);
         if( sl.compare("pp-raw")==0){
	    data.ppRaw          = std::atof( s1.c_str() ); 
	    data.ppRaw_err      = std::atof( s1err.c_str() );
	    data.ppRaw_fxpr     = std::atof( s2.c_str() ); 
	    data.ppRaw_fxpr_err = std::atof( s2err.c_str() ); 
	    data.ppRaw_trly     = std::atof( s3.c_str() ); 
	    data.ppRaw_trly_err = std::atof( s3err.c_str() ); 
	 }else if(sl.compare("pp-free")==0){
	    data.ppFree          = std::atof( s1.c_str() ); 
	    data.ppFree_err      = std::atof( s1err.c_str() );
	    data.ppFree_fxpr     = std::atof( s2.c_str() ); 
	    data.ppFree_fxpr_err = std::atof( s2err.c_str() ); 
	    data.ppFree_trly     = std::atof( s3.c_str() ); 
	    data.ppFree_trly_err = std::atof( s3err.c_str() ); 
         }else if(sl.compare("trly")==0){
	    data.trly          = std::atof( s1.c_str() ); 
	    data.trly_err      = std::atof( s1err.c_str() );
	    data.trly_fxpr     = std::atof( s2.c_str() ); 
	    data.trly_fxpr_err = std::atof( s2err.c_str() ); 
	    data.trly_trly     = std::atof( s3.c_str() ); 
	    data.trly_trly_err = std::atof( s3err.c_str() ); 
         }else if(sl.compare("drift-shim")==0){
	    data.driftShim          = std::atof( s1.c_str() ); 
	    data.driftShim_err      = std::atof( s1err.c_str() );
	    data.driftShim_fxpr     = std::atof( s2.c_str() ); 
	    data.driftShim_fxpr_err = std::atof( s2err.c_str() ); 
	    data.driftShim_trly     = std::atof( s3.c_str() ); 
	    data.driftShim_trly_err = std::atof( s3err.c_str() ); 
         }else if(sl.compare("diff")==0){
	    data.diff          = std::atof( s1.c_str() ); 
	    data.diff_err      = std::atof( s1err.c_str() );
	    data.diff_fxpr     = std::atof( s2.c_str() ); 
	    data.diff_fxpr_err = std::atof( s2err.c_str() ); 
	    data.diff_trly     = std::atof( s3.c_str() ); 
	    data.diff_trly_err = std::atof( s3err.c_str() ); 
         }
      }
      infile.close();
   }

   return 0;
}
