#include "../include/FXPRFuncs.h"
//______________________________________________________________________________
int GetAverageFXPRVectorsAlt(int method,
                             unsigned long long runEndPoint,unsigned long long runStartPoint,unsigned long long tStep,
                             std::vector<int> probeList,std::vector<gm2field::fixedProbeFrequency_t> fxprData,
                             std::vector<unsigned long long> &time,std::vector<double> &freq){
   // fill vectors with all fxpr data based on time slices defined by run start/stop points 
   // (non-continuous data in time) 

   const int NPTS = fxprData.size();
   const int NFP  = probeList.size(); 
   int startIndex = probeList[0]; 
   int stopIndex  = probeList[NFP-1]; 

   // grab data for first run (signified by runEndPoint) 
   unsigned long long t0     = 0;  
   unsigned long long tStart = fxprData[0].GpsTimeStamp[startIndex];
   unsigned long long tStop  = runEndPoint;
   GetAverageFXPRVectors(method,t0,tStart,tStop,tStep,probeList,fxprData,time,freq);

   // grab data for last run (signified by runStartPoint)  
   double mean=0,stdev=0,mean_last=0;
   tStart = runStartPoint;
   tStop  = fxprData[NPTS-1].GpsTimeStamp[stopIndex];
   GetAverageFXPRVectors(method,t0,tStart,tStop,tStep,probeList,fxprData,time,freq);
 
   // for(int i=0;i<NPTS;i++){
   //    theTime = tStart + i*tStep;
   //    GetAverageFXPR(method,theTime,probeList,fxprData,mean,stdev);
   //    if( fabs(mean-mean_last)>60 && i!=0 ){
   //       std::cout << Form("FXPR Avg is suspicious! mean = %.3lf Hz, last mean = %.3lf Hz",mean,mean_last) << std::endl;
   //       mean = mean_last;
   //       continue;   // skip this event 
   //    }
   //    arg = theTime - t0; 
   //    time.push_back(arg);
   //    freq.push_back(mean);
   //    mean_last = mean;
   // }
   return 0;
}
//______________________________________________________________________________
int GetAverageFXPRVectorsNew(int method,
                             unsigned long long t0,unsigned long long tStart,unsigned long long tStop,unsigned long long tStep, 
                             std::vector<int> probeList,std::vector<gm2field::fixedProbeFrequency_t> fxprData,
                             std::vector<fixedProbeEvent_t> &fxprDataAvg){
   // fill vectors with fxpr avg at every event 
   // just take straight average of the probes we care about 
   // FIXME: change the spacing of the data? 
   fixedProbeEvent_t data;  
   TStopwatch *watch = new TStopwatch();
   const int NEvents = fxprData.size(); 
   const int NPTS    = (tStop-tStart)/tStep;
   const int NPR     = probeList.size(); 
   unsigned long long theTime=0,arg=0;
   double mean_t=0,mean_f=0,arg_t=0,arg_f=0;
   std::vector<double> T,F;
   watch->Start(); 
   int probeIndex=0; 
   for(int i=0;i<NEvents;i++){
      for(int j=0;j<NPR;j++){
	 probeIndex = probeList[j]; 
	 arg_t = fxprData[i].GpsTimeStamp[probeIndex] - t0;
         arg_f = fxprData[i].Frequency[probeIndex][method]; 
	 T.push_back(arg_t/1E+9);  
	 F.push_back(arg_f);  
      }
      // compute averages 
      mean_t = gm2fieldUtil::Math::GetMean<double>(T); 
      mean_f = gm2fieldUtil::Math::GetMean<double>(F);      
      // store results and clear vectors 
      data.time    = mean_t; 
      data.freq    = mean_f;
      data.freqErr = 0; 
      fxprDataAvg.push_back(data);  
      T.clear(); 
      F.clear(); 
      // show time profile 
      watch->Stop();
      // std::cout << Form("Event %d complete (duration: %.3E sec)",i+1,watch->RealTime()) << std::endl;
      watch->Start();
   }
   return 0;
}
//______________________________________________________________________________
int GetAverageFXPRVectorsNew(int method,
                             unsigned long long t0,unsigned long long tStart,unsigned long long tStop,unsigned long long tStep, 
                             std::vector<int> probeList,std::vector<gm2field::fixedProbeFrequency_t> fxprData,
                             std::vector<unsigned long long> &time,std::vector<double> &freq){
   // fill vectors with fxpr avg at every event 
   // just take straight average of the probes we care about  
   TStopwatch *watch = new TStopwatch();
   const int NEvents = fxprData.size(); 
   const int NPTS    = (tStop-tStart)/tStep;
   const int NPR     = probeList.size(); 
   unsigned long long theTime=0,arg=0;
   double mean=0,stdev=0,mean_last=0;
   double mean_t=0,mean_f=0,arg_t=0,arg_f=0;
   std::vector<unsigned long long> T; 
   std::vector<double> F;
   // std::cout << "Processing " << NPTS << " events..." << std::endl;
   watch->Start(); 
   int probeIndex=0; 
   for(int i=0;i<NEvents;i++){
      for(int j=0;j<NPR;j++){
	 probeIndex = probeList[j]; 
	 arg_t = fxprData[i].GpsTimeStamp[probeIndex] - t0;
         arg_f = fxprData[i].Frequency[probeIndex][method]; 
	 T.push_back(arg_t);  
	 F.push_back(arg_f);  
      }
      // compute averages 
      mean_t = gm2fieldUtil::Math::GetMean<unsigned long long>(T); 
      mean_f = gm2fieldUtil::Math::GetMean<double>(F);      
      // store results and clear vectors  
      time.push_back(mean_t);
      freq.push_back(mean_f);
      T.clear(); 
      F.clear(); 
      watch->Stop();
      // std::cout << Form("Event %d complete (duration: %.3E sec)",i+1,watch->RealTime()) << std::endl;
      watch->Start();
   }
   return 0;
}
//______________________________________________________________________________
int GetAverageFXPRVectors(int method,
                          unsigned long long t0,unsigned long long tStart,unsigned long long tStop,unsigned long long tStep, 
                          std::vector<int> probeList,std::vector<gm2field::fixedProbeFrequency_t> fxprData,
                          std::vector<unsigned long long> &time,std::vector<double> &freq){
   // fill vectors with fxpr avg at specified intervals in time 
   TStopwatch *watch = new TStopwatch();
   const int NPTS = (tStop-tStart)/tStep;
   unsigned long long theTime=0,arg=0;
   double mean=0,stdev=0,mean_last=0;
   // std::cout << "Processing " << NPTS << " events..." << std::endl;
   watch->Start();  
   for(int i=0;i<NPTS;i++){
      theTime = tStart + i*tStep;
      GetAverageFXPR(method,theTime,probeList,fxprData,mean,stdev);
      if( fabs(mean-mean_last)>60 && i!=0 ){
         std::cout << Form("FXPR Avg is suspicious! mean = %.3lf Hz, last mean = %.3lf Hz",mean,mean_last) << std::endl;
         mean = mean_last;
         continue;   // skip this event 
      }
      arg = theTime - t0; 
      time.push_back(arg);
      freq.push_back(mean);
      mean_last = mean;
      watch->Stop();
      // std::cout << Form("Event %d complete (duration: %.3lf sec)",i+1,watch->RealTime()) << std::endl;
      watch->Start();
   }
   return 0;
}
//______________________________________________________________________________
int GetAverageFXPR(int method,unsigned long long time,std::vector<int> probe,
                   std::vector<gm2field::fixedProbeFrequency_t> fxprData,
                   int NN,double *TIME,double *FREQ,double &mean,double &stdev){
   // using the list of fixed probes, find the average fixed probe field at a specific time
   // make sure we have data
   int eventNum=0;
   const int N = probe.size();
   if(N==0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   // first check if the time is within the bounds we have 
   unsigned long long theTime = 0;  
   unsigned long long tStart = 9E+10;
   unsigned long long tStop  = 0; 
   for(int i=0;i<N;i++){
      // check first event 
      theTime = fxprData[0].GpsTimeStamp[i];
      if(theTime<tStart) tStart = theTime;
      if(theTime>tStop)  tStop  = theTime;
      // check last event 
      theTime = fxprData[NN-1].GpsTimeStamp[i];
      if(theTime<tStart) tStart = theTime;
      if(theTime>tStop)  tStop  = theTime;
   }  

   // now check bracketing range 
   bool inRange=false,isBefore=false,isAfter=false;
 
   if(time<tStart){
      isBefore = true;
   }else if( time>=tStart && time<=tStop ){
      inRange = true;
   }else if(time>tStop){
      isAfter = true;
   }

   // std::cout << Form("start = %.3lf, time = %.3lf, stop = %.3lf",tStart/1E+9,time/1E+9,tStop/1E+9) << std::endl; 
  
   double t_star=0;
   double t_k=0,t_km1=0;
   double y_k=0,y_km1=0;
   double r=0;

   // gather the events we care about 
   double theFreq=0;
   std::vector<double> freq;
   if(inRange){
      // need linear interpolation  
      for(int i=0;i<N;i++){
	 // eventNum = gm2fieldUtil::Algorithm::FindFPEvent(time,fxprData,probe[i]); 
         theFreq = GetInterpolatedFXPRFreq(method,probe[i],time,fxprData,NN,TIME,FREQ);
	 freq.push_back(theFreq); 
	 // std::cout << Form("probe %d, event %d: %.3lf Hz",probe[i],eventNum,fxprData[eventNum].Frequency[probe[i]][method]) << std::endl;
      }
      // calculate the mean and standard deviation 
      mean  = gm2fieldUtil::Math::GetMean<double>(freq); 
      stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(freq); 
   }else if(isBefore){
      // need linear extrapolation 
      for(int i=0;i<N;i++) freq.push_back(fxprData[0].Frequency[probe[i]][method]); 
      y_k      = gm2fieldUtil::Math::GetMean<double>(freq); 
      freq.clear(); 
      for(int i=0;i<N;i++) freq.push_back(fxprData[1].Frequency[probe[i]][method]); 
      y_km1    = gm2fieldUtil::Math::GetMean<double>(freq); 
      stdev    = gm2fieldUtil::Math::GetStandardDeviation<double>(freq); 
      t_star   = time/1E+9; 
      t_k      = fxprData[0].GpsTimeStamp[probe[0]]/1E+9;
      t_km1    = fxprData[1].GpsTimeStamp[probe[0]]/1E+9;
      r        = (t_star-t_km1)/(t_k-t_km1);
      mean     = y_km1 + r*(y_k-y_km1); 
      stdev    = 0;
   }else if(isAfter){
      // need linear extrapolation 
      for(int i=0;i<N;i++) freq.push_back(fxprData[NN-1].Frequency[probe[i]][method]); 
      y_k      = gm2fieldUtil::Math::GetMean<double>(freq); 
      freq.clear(); 
      for(int i=0;i<N;i++) freq.push_back(fxprData[NN-2].Frequency[probe[i]][method]); 
      y_km1    = gm2fieldUtil::Math::GetMean<double>(freq); 
      stdev    = gm2fieldUtil::Math::GetStandardDeviation<double>(freq); 
      t_star   = time/1E+9; 
      t_k      = fxprData[NN-1].GpsTimeStamp[probe[N-1]]/1E+9;
      t_km1    = fxprData[NN-2].GpsTimeStamp[probe[N-1]]/1E+9;
      r        = (t_star-t_km1)/(t_k-t_km1);
      mean     = y_km1 + r*(y_k-y_km1); 
      stdev    = 0;
      // std::cout << "LINEAR EXTRAPOLATION (AFTER)" << std::endl;
      // std::cout << Form("t_k = %.3lf, y_k = %.3lf,t_k-1 = %.3lf, y_k-1 = %.3lf,t_* = %.3lf, y_* = %.3lf",t_k,y_k,t_km1,y_km1,t_star,mean) << std::endl; 
   }

   return 0; 
}
//______________________________________________________________________________
int GetAverageFXPR(double time,std::vector<fixedProbeEvent_t> fxprData,double &mean,double &stdev){
   // using the average fixed probes, find the average fixed probe field at a specific time

   // gather times and frequencies
   const int NEvents = fxprData.size();  
   std::vector<double> T,F; 
   for(int i=0;i<NEvents;i++){
      T.push_back( fxprData[i].time ); 
      F.push_back( fxprData[i].freq ); 
   }

   // check bracketing range 
   bool inRange=false,isBefore=false,isAfter=false;

   double tStart = fxprData[0].time; 
   double tStop  = fxprData[NEvents-1].time;
 
   if(time<tStart){
      isBefore = true;
   }else if( time>=tStart && time<=tStop ){
      inRange = true;
   }else if(time>tStop){
      isAfter = true;
   }
  
   double t_star=0;
   double t_k=0,t_km1=0;
   double y_k=0,y_km1=0;
   double r=0;

   int N = 2;  
   int lo=0,hi=0,start=0,stop=0;

   // gather the events we care about 
   double theFreq=0;
   std::vector<double> freq;
   if(inRange){
      // get average over neighboring +/- N events
      // first find the index we need 
      gm2fieldUtil::Algorithm::BinarySearch<double>(T,time,lo,hi);   
      start = lo - N; 
      stop  = hi + N;
      if(start<0) start = 0;
      if(stop>NEvents) stop = NEvents; 
      for(int i=start;i<stop;i++) freq.push_back(F[i]); 
      // calculate the mean and standard deviation 
      mean  = gm2fieldUtil::Math::GetMean<double>(freq); 
      stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(freq); 
   }else if(isBefore){
      // need linear extrapolation 
      y_k      = F[1]; 
      y_km1    = F[0]; 
      t_star   = time; 
      t_k      = T[1];
      t_km1    = T[0];
      r        = (t_star-t_km1)/(t_k-t_km1);
      mean     = y_km1 + r*(y_k-y_km1); 
      stdev    = 0;
   }else if(isAfter){
      // need linear extrapolation 
      y_k      = F[NEvents-1]; 
      y_km1    = F[NEvents-2]; 
      t_star   = time; 
      t_k      = T[NEvents-1];
      t_km1    = T[NEvents-2];
      r        = (t_star-t_km1)/(t_k-t_km1);
      mean     = y_km1 + r*(y_k-y_km1); 
      stdev    = 0;
   }

   return 0; 
}
//______________________________________________________________________________
int GetAverageFXPR(int method,unsigned long long time,std::vector<int> probe,
                   std::vector<gm2field::fixedProbeFrequency_t> fxprData,double &mean,double &stdev){
   // using the list of fixed probes, find the average fixed probe field at a specific time
   // make sure we have data
   int eventNum=0;
   const int N = probe.size();
   if(N==0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   // first check if the time is within the bounds we have 
   int NEvents = fxprData.size();
   unsigned long long theTime = 0;  
   unsigned long long tStart = 9E+10;
   unsigned long long tStop  = 0; 
   for(int i=0;i<N;i++){
      // check first event 
      theTime = fxprData[0].GpsTimeStamp[i];
      if(theTime<tStart) tStart = theTime;
      if(theTime>tStop)  tStop  = theTime;
      // check last event 
      theTime = fxprData[NEvents-1].GpsTimeStamp[i];
      if(theTime<tStart) tStart = theTime;
      if(theTime>tStop)  tStop  = theTime;
   }  

   // now check bracketing range 
   bool inRange=false,isBefore=false,isAfter=false;
 
   if(time<tStart){
      isBefore = true;
   }else if( time>=tStart && time<=tStop ){
      inRange = true;
   }else if(time>tStop){
      isAfter = true;
   }

   // std::cout << Form("start = %.3lf, time = %.3lf, stop = %.3lf",tStart/1E+9,time/1E+9,tStop/1E+9) << std::endl; 
  
   double t_star=0;
   double t_k=0,t_km1=0;
   double y_k=0,y_km1=0;
   double r=0;

   // gather the events we care about 
   double theFreq=0;
   std::vector<double> freq;
   if(inRange){
      // need linear interpolation  
      for(int i=0;i<N;i++){
	 // eventNum = gm2fieldUtil::Algorithm::FindFPEvent(time,fxprData,probe[i]); 
         theFreq = GetInterpolatedFXPRFreq(method,probe[i],time,fxprData);
	 freq.push_back(theFreq); 
	 // std::cout << Form("probe %d, event %d: %.3lf Hz",probe[i],eventNum,fxprData[eventNum].Frequency[probe[i]][method]) << std::endl;
      }
      // calculate the mean and standard deviation 
      mean  = gm2fieldUtil::Math::GetMean<double>(freq); 
      stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(freq); 
   }else if(isBefore){
      // need linear extrapolation 
      for(int i=0;i<N;i++) freq.push_back(fxprData[0].Frequency[probe[i]][method]); 
      y_k      = gm2fieldUtil::Math::GetMean<double>(freq); 
      freq.clear(); 
      for(int i=0;i<N;i++) freq.push_back(fxprData[1].Frequency[probe[i]][method]); 
      y_km1    = gm2fieldUtil::Math::GetMean<double>(freq); 
      stdev    = gm2fieldUtil::Math::GetStandardDeviation<double>(freq); 
      t_star   = time/1E+9; 
      t_k      = fxprData[0].GpsTimeStamp[probe[0]]/1E+9;
      t_km1    = fxprData[1].GpsTimeStamp[probe[0]]/1E+9;
      r        = (t_star-t_km1)/(t_k-t_km1);
      mean     = y_km1 + r*(y_k-y_km1); 
      stdev    = 0;
   }else if(isAfter){
      // need linear extrapolation 
      for(int i=0;i<N;i++) freq.push_back(fxprData[NEvents-1].Frequency[probe[i]][method]); 
      y_k      = gm2fieldUtil::Math::GetMean<double>(freq); 
      freq.clear(); 
      for(int i=0;i<N;i++) freq.push_back(fxprData[NEvents-2].Frequency[probe[i]][method]); 
      y_km1    = gm2fieldUtil::Math::GetMean<double>(freq); 
      stdev    = gm2fieldUtil::Math::GetStandardDeviation<double>(freq); 
      t_star   = time/1E+9; 
      t_k      = fxprData[NEvents-1].GpsTimeStamp[probe[N-1]]/1E+9;
      t_km1    = fxprData[NEvents-2].GpsTimeStamp[probe[N-1]]/1E+9;
      r        = (t_star-t_km1)/(t_k-t_km1);
      mean     = y_km1 + r*(y_k-y_km1); 
      stdev    = 0;
      // std::cout << "LINEAR EXTRAPOLATION (AFTER)" << std::endl;
      // std::cout << Form("t_k = %.3lf, y_k = %.3lf,t_k-1 = %.3lf, y_k-1 = %.3lf,t_* = %.3lf, y_* = %.3lf",t_k,y_k,t_km1,y_km1,t_star,mean) << std::endl; 
   }

   return 0; 
}
//______________________________________________________________________________
double GetInterpolatedFXPRFreq(int method,int probe,unsigned long long time,
                               std::vector<gm2field::fixedProbeFrequency_t> fxprData,
                               int NN,double *TIME,double *FREQ){
   // linear interpolated frequency for a given fixed probe
   // gather all data for the probe
   for(int i=0;i<NN;i++){
      TIME[i] = fxprData[i].GpsTimeStamp[probe]/1E+9; 
      FREQ[i] = fxprData[i].Frequency[probe][method]; 
   }
   // now find the boundaries in time 
   int lo=0,hi=0;
   gm2fieldUtil::Algorithm::BinarySearch<double>(NN,TIME,time/1E+9,lo,hi);
   // now get the mean frequency from a linear interpolation 
   double freq = gm2fieldUtil::Math::LinearInterpolation(time/1E+9,TIME[lo],FREQ[lo],TIME[hi],FREQ[hi]); 
   return freq;  
}
//______________________________________________________________________________
double GetInterpolatedFXPRFreq(int method,int probe,unsigned long long time,std::vector<gm2field::fixedProbeFrequency_t> fxprData){
   // linear interpolated frequency for a given fixed probe
   // gather all data for the probe
   const int N = fxprData.size();
   std::vector<double> TIME,FREQ; 
   for(int i=0;i<N;i++){
      TIME.push_back( fxprData[i].GpsTimeStamp[probe]/1E+9 ); 
      FREQ.push_back( fxprData[i].Frequency[probe][method] ); 
   }
   // now find the boundaries in time 
   int lo=0,hi=0;
   gm2fieldUtil::Algorithm::BinarySearch<double>(TIME,time/1E+9,lo,hi);
   // now get the mean frequency from a linear interpolation 
   double freq = gm2fieldUtil::Math::LinearInterpolation(time/1E+9,TIME[lo],FREQ[lo],TIME[hi],FREQ[hi]); 
   return freq;  
}
//______________________________________________________________________________
int FitFXPR(int lo,int hi,std::vector<fixedProbeEvent_t> fxprData,double *par){
   // fit the Fixed proe data to a line based on start and stop indices  
   std::vector<double> T,F;
   const int N = fxprData.size();
   int start = lo;
   int stop  = hi;
   if(hi>N) stop  = N;
   if(lo<0) start = 0;
   for(int i=start;i<stop;i++){
      T.push_back( fxprData[i].time ); 
      F.push_back( fxprData[i].freq ); 
   }

   double intercept=0,slope=0,r=0;
   gm2fieldUtil::Math::LeastSquaresFitting(T,F,intercept,slope,r);

   par[0] = intercept;
   par[1] = slope;

   return 0;
}
//______________________________________________________________________________
int FilterSingle(int nev,double T,std::vector<fixedProbeEvent_t> in,std::vector<fixedProbeEvent_t> &out){

   std::vector<double> ff,tt;
   const int N = in.size();
   for(int i=0;i<N;i++){
      tt.push_back(in[i].time); 
      ff.push_back(in[i].freq); 
   }

   int lo=0,hi=0;
   gm2fieldUtil::Algorithm::BinarySearch(tt,T,lo,hi); 
   int end   = lo-10;
   int start = end-nev;
   for(int i=start;i<end;i++){
      out.push_back(in[i]);
   } 

   return 0;

}
// //______________________________________________________________________________
// int Filter(int type,std::vector<double> T,
//            std::vector<fixedProbeEvent_t> in,std::vector<fixedProbeEvent_t> &out){
// 
//    const int N = in.size();
//    const int M = T.size();  
// 
//    std::vector<double> ff,tt;
//    for(int i=0;i<N;i++){
//       tt.push_back(in[i].time); 
//       ff.push_back(in[i].freq); 
//    }
//    double mean = gm2fieldUtil::Math::GetMean<double>(x); 
// 
//    int lo=0,hi=0;
//    int eventNum=0;
//    double aTime=0; 
//    for(int i=0;i<M;i++){     // transition loop
//       gm2fieldUtil::Algorithm::BinarySearch(tt,T[i],lo,hi); 
//       end = lo;
//       start = lo-30;
//       for(int j=start;j<end;j++){
// 	 out.push_back(in[j]);
//       } 
//  
//    }
// 
//    return 0;
// 
// }
//______________________________________________________________________________
int FindTransitionTimes(double thr,std::vector<fixedProbeEvent_t> fxprData,
                        std::vector<double> &timeLo,std::vector<double> &timeHi){

   std::string theTime;
   double aTime=0,aFreq=0,diff=0,freq=0;
   double freq_last = fxprData[0].freq;

   // find the mean 
   std::vector<double> x;
   const int N = fxprData.size(); 
   for(int i=0;i<N;i++) x.push_back( fxprData[i].freq ); 
   double mean = gm2fieldUtil::Math::GetMean<double>(x); 

   // split into a low list and a high list 
   std::vector<double> tHi,fHi,tLo,fLo; 
   for(int i=0;i<N;i++){
      aTime = fxprData[i].time; 
      aFreq = fxprData[i].freq; 
      if(aFreq>mean){
	 tHi.push_back(aTime);  
	 fHi.push_back(aFreq);
      }else{
	 tLo.push_back(aTime);  
	 fLo.push_back(aFreq);
      }  
   }

   int cntr=0;

   // find peaks 
   int NH = tHi.size();
   int i=0;
   do{
      diff = TMath::Abs(fHi[i]-fHi[i-1]);
      if(diff<thr){
	 // probably good data; keep moving 
	 i++;
      }else{
	 // probably an endpoint! mark it 
	 cntr++;
	 if(cntr>1) timeHi.push_back(tHi[i]);
	 i += 30;
      }
   }while(i<NH+1); 

   cntr=0; 
   // find peaks
   i=0;
   int NL = tLo.size();
   do{
      diff = TMath::Abs(fLo[i]-fLo[i-1]);
      if(diff<thr){
	 // probably good data; keep moving 
	 i++;
      }else{
	 // probably an endpoint! mark it 
	 cntr++;
	 if(cntr>1) timeLo.push_back(tLo[i]);
	 i += 30;
      }
   }while(i<NL+1); 

   return 0;
}
