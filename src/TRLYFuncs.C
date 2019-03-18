#include "../include/TRLYFuncs.h"
//______________________________________________________________________________
int FindTransitionTimes(int step,double thr,std::vector<trolleyAnaEvent_t> Data,
                        std::vector< std::vector<double> > &timeLo,std::vector< std::vector<double> > &timeHi){
   // find transition times for all trolley probes 
   int rc=0;
   std::vector<double> tLo,tHi; 
   const int N = 17;
   for(int i=0;i<N;i++){
      rc = FindTransitionTimes(i,step,thr,Data,tLo,tHi);
      timeLo.push_back(tLo); 
      timeHi.push_back(tHi); 
   } 
   return 0;
}
//______________________________________________________________________________
int FindTransitionTimes(int probe,int step,double thr,std::vector<trolleyAnaEvent_t> Data,
                        std::vector<double> &timeLo,std::vector<double> &timeHi){

   // thr  = threshold for change (Hz) 
   // step = how many events to skip after finding a transition 

   double aTime=0,aFreq=0,diff,freq=0;

   // find the mean 
   std::vector<double> x;
   const int N = Data.size();
   for(int i=0;i<N;i++) x.push_back( Data[i].freq[probe] );
   double mean = gm2fieldUtil::Math::GetMean<double>(x);

   // split into a low list and a high list 
   std::vector<double> tHi,fHi,tLo,fLo;
   for(int i=0;i<N;i++){
      aTime = Data[i].time[probe]/1E+9;
      aFreq = Data[i].freq[probe];
      if(aFreq>=mean){
         tHi.push_back(aTime);
         fHi.push_back(aFreq);
      }else{
         tLo.push_back(aTime);
         fLo.push_back(aFreq);
      }
   }

   int cntr=0;
   // find high peaks 
   int NH = tHi.size();
   int i=1;
   do{
      diff = TMath::Abs(fHi[i]-fHi[i-1]);
      if(diff<thr){
         // probably good data; keep moving 
         i++;
      }else{
         // probably an endpoint! mark it 
         cntr++;
         // if(cntr%2!=0) timeHi.push_back(tHi[i]); 
	 if(cntr>0){
	    timeHi.push_back(tHi[i]);
	 }
	 // if(cntr==1){
	 //    std::cout << "First high transition: " << gm2fieldUtil::GetStringTimeStampFromUTC(tHi[i]) << std::endl;
         // }
         i += step;
      }
   }while(i<NH+1);

   // find low peaks
   cntr=0;
   i=1;
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
         // if(cntr%2==0) timeLo.push_back(tLo[i]);
	 i += step;
      }
   }while(i<NL+1);

   // there's one more to add on at the end
   // timeLo.push_back(tLo[NL-1]); 

   return 0;
}
//______________________________________________________________________________
int FilterSingle(std::string var,int probe,int nev,double T,std::vector<trolleyAnaEvent_t> in,
                 std::vector<double> &x){

   std::vector<double> tt;
   const int N = in.size();
   for(int i=0;i<N;i++) tt.push_back(in[i].time[probe]/1E+9);

   int lo=0,hi=0;
   gm2fieldUtil::Algorithm::BinarySearch(tt,T,lo,hi);
   int end   = lo;
   int start = end-nev;
   for(int i=start;i<end;i++){
      if( var.compare("time")==0 ) x.push_back(in[i].time[probe]/1E+9); 
      if( var.compare("freq")==0 ) x.push_back(in[i].freq[probe]     );
      if( var.compare("temp")==0 ) x.push_back(in[i].temp[probe]     );
      if( var.compare("r")==0    ) x.push_back(in[i].r[probe]        );
      if( var.compare("y")==0    ) x.push_back(in[i].y[probe]        );
      if( var.compare("phi")==0  ) x.push_back(in[i].phi[probe]      );
   }

   return 0;
}
//______________________________________________________________________________
int GetTRLYStatsAtTime(bool UseTempCor,int probe,int nev,double fLO,std::vector<double> time,
                       std::vector<trolleyAnaEvent_t> Data,std::vector<calibSwap_t> &Event){

   // find the mean field at the times specified in the time vector 

   calibSwap_t theEvent; 

   const int N = time.size();
   int M=0,rc=0;
   double delta=0;
   double mean_freq=0,stdev_freq=0;
   double mean_temp=0,stdev_temp=0;
   double mean_x=0,stdev_x=0;
   double mean_y=0,stdev_y=0;
   double mean_z=0,stdev_z=0;
   std::vector<double> freq,temp,x,y,z;
   for(int i=0;i<N;i++){
      // std::cout << "Looking for time " << gm2fieldUtil::GetStringTimeStampFromUTC(time[i]) << std::endl;
      // find events 
      // rc = FilterSingle("time",probe,nev,time[i],Data,time);
      rc = FilterSingle("freq",probe,nev,time[i],Data,freq);
      rc = FilterSingle("temp",probe,nev,time[i],Data,temp);
      rc = FilterSingle("r"   ,probe,nev,time[i],Data,x);
      rc = FilterSingle("y"   ,probe,nev,time[i],Data,y);
      rc = FilterSingle("phi" ,probe,nev,time[i],Data,z);
      // now add in the LO 
      M = freq.size();
      for(int j=0;j<M;j++){
	 if(UseTempCor){
	    // apply a temperature correction if necessary
	    // accounts for the trolley being at a temperature other than 25 deg c  
	    delta = (10.36E-9)*(temp[j]-25.0);
	 }
	 freq[j] = (freq[j]+fLO)/(1. - delta); 
      }
      // now get mean of events 
      mean_freq  = gm2fieldUtil::Math::GetMean<double>(freq);
      stdev_freq = gm2fieldUtil::Math::GetStandardDeviation<double>(freq);
      mean_temp  = gm2fieldUtil::Math::GetMean<double>(temp);
      stdev_temp = gm2fieldUtil::Math::GetStandardDeviation<double>(temp);
      mean_x     = gm2fieldUtil::Math::GetMean<double>(x);
      stdev_x    = gm2fieldUtil::Math::GetStandardDeviation<double>(x);
      mean_y     = gm2fieldUtil::Math::GetMean<double>(y);
      stdev_y    = gm2fieldUtil::Math::GetStandardDeviation<double>(y);
      mean_z     = gm2fieldUtil::Math::GetMean<double>(z);
      stdev_z    = gm2fieldUtil::Math::GetStandardDeviation<double>(z);
      // store result
      theEvent.time    = time[i]; 
      theEvent.freq    = mean_freq;
      theEvent.freqErr = stdev_freq; 
      theEvent.temp    = mean_temp; 
      theEvent.tempErr = stdev_temp; 
      theEvent.r       = mean_x;       
      theEvent.rErr    = stdev_x;       
      theEvent.y       = mean_y;       
      theEvent.yErr    = stdev_y;       
      theEvent.phi     = mean_z;       
      theEvent.phiErr  = stdev_z;       
      Event.push_back(theEvent);  
      // set up for next time 
      freq.clear();
      temp.clear();
      x.clear();
      y.clear();
      z.clear();
   }

   return 0;
}
//______________________________________________________________________________
int GetTRLYStats_sccToggle(int probe,int nev,std::vector<double> time,std::vector<trolleyAnaEvent_t> Data,
                           std::vector<double> &TIME,std::vector<double> &MEAN,std::vector<double> &STDEV){

   // find the mean field at the times specified in the time vector 

   const int N = time.size();
   int M=0,rc=0;
   double mean=0,stdev=0;
   std::vector<double> freq;
   for(int i=0;i<N;i++){
      // std::cout << "Looking for time " << gm2fieldUtil::GetStringTimeStampFromUTC(time[i]) << std::endl;
      // find events 
      rc = FilterSingle("freq",probe,nev,time[i],Data,freq);
      // now get mean of events 
      mean  = gm2fieldUtil::Math::GetMean<double>(freq);
      stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(freq);
      // store result
      TIME.push_back(time[i]); 
      MEAN.push_back(mean);
      STDEV.push_back(stdev);
      // set up for next time 
      freq.clear();
   }

   return 0;
}
