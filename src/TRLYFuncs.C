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
   
   // std::cout << "[FilterSingle]: Looking for variable " << var << std::endl;

   int lo=0,hi=0,start=0,end=0;
   if(nev<0){
      start = 0;
      end   = N;
   }else{
      gm2fieldUtil::Algorithm::BinarySearch(tt,T,lo,hi);
      end   = lo;
      start = end-nev;
   }

   // have to adjust entries if we pick up the start of the data set  
   if(lo==0){
      start = 0;
      end   = nev; 
   }

   // if the number of events is LARGER than the lowest index found 
   if(lo>0 && lo<nev){
      start = 0; // just lose the extra events we'd want  
      end   = lo; 
   }

   // if( var.compare("time")==0 ) { 
   //    std::cout << "[FilterSingle]: The KEY is " << gm2fieldUtil::GetStringTimeStampFromUTC(T) << std::endl;
   //    std::cout << "[FilterSingle]: lo index = " << lo << " hi index = " << hi << std::endl;
   //    std::cout << "[FilterSingle]: start index = " << start << " end index = " << end << std::endl;
   // }

   for(int i=start;i<end;i++){
      if( var.compare("time")==0 ) x.push_back(in[i].time[probe]/1E+9); 
      if( var.compare("freq")==0 ) x.push_back(in[i].freq[probe]     );
      if( var.compare("temp")==0 ) x.push_back(in[i].temp[probe]     );
      if( var.compare("r")==0    ) x.push_back(in[i].r[probe]        );
      if( var.compare("y")==0    ) x.push_back(in[i].y[probe]        );
      if( var.compare("phi")==0  ) x.push_back(in[i].phi[probe]      );
   }

   // const int NX = x.size();
   // if( var.compare("time")==0 ) { 
   //    std::cout << "[FilterSingle]: Accumulated events: " << NX << std::endl;
   //    for(int i=0;i<NX;i++){
   //       std::cout << Form("--> event %03d, %s",i,gm2fieldUtil::GetStringTimeStampFromUTC(x[i]).c_str()) << std::endl;
   //    }
   // }

   return 0;
}
//______________________________________________________________________________
int GetTRLYStatsAtTime_old(bool UseTempCor,bool UseOscCor,int probe,int nev,double fLO,
                           std::vector<double> time,std::vector<averageFixedProbeEvent_t> fxpr,
                           std::vector<trolleyAnaEvent_t> Data,std::vector<calibSwap_t> &Event){

   // find the mean field at the times specified in the time vector 

   calibSwap_t theEvent; 

   // first gather the FXPR data 
   std::vector<double> tt,ff;
   const int NFP = fxpr.size();
   for(int i=0;i<NFP;i++){
      tt.push_back( fxpr[i].time/1E+9 );
      ff.push_back( fxpr[i].freq );
   }

   const int N = time.size();
   int M=0,rc=0,lo=0,hi=0;
   double f0,theFreq,stdev; 
   double delta_t=0,delta_osc=0;
   double mean_freq=0,stdev_freq=0;
   double mean_temp=0,stdev_temp=0;
   double mean_x=0,stdev_x=0;
   double mean_y=0,stdev_y=0;
   double mean_z=0,stdev_z=0;
   std::vector<double> trt,freq,temp,x,y,z;
   for(int i=0;i<N;i++){
      // std::cout << "Looking for time " << gm2fieldUtil::GetStringTimeStampFromUTC(time[i]) << std::endl;
      // find events 
      // rc = FilterSingle("time",probe,nev,time[i],Data,time);
      rc = FilterSingle("time",probe,nev,time[i],Data,trt);
      rc = FilterSingle("freq",probe,nev,time[i],Data,freq);
      rc = FilterSingle("temp",probe,nev,time[i],Data,temp);
      rc = FilterSingle("r"   ,probe,nev,time[i],Data,x);
      rc = FilterSingle("y"   ,probe,nev,time[i],Data,y);
      rc = FilterSingle("phi" ,probe,nev,time[i],Data,z);
      // now add in the LO 
      M = freq.size();
      for(int j=0;j<M;j++){
	 if(UseTempCor){
	    // FIXME: apply a temperature correction if necessary
	    // accounts for the trolley being at a temperature other than 25 deg c  
	    delta_t = (4.0E-9)*(temp[j]-25.0);
	 }
         if(UseOscCor){
	    gm2fieldUtil::Algorithm::BinarySearch<double>(tt,trt[j],lo,hi);
	    // take average over +/- 5 events in time 
	    rc = GetStats_vec(lo-5,hi+5,ff,theFreq,stdev);
	    if(j==0) f0 = theFreq;
	    delta_osc = theFreq - f0; 
         }
	 freq[j] = (freq[j]+fLO-delta_osc)/(1. - delta_t); 
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
      trt.clear();
      freq.clear();
      temp.clear();
      x.clear();
      y.clear();
      z.clear();
   }

   return 0;
}
//______________________________________________________________________________
int GetTRLYStatsAtTime(bool UseTempCor,bool UseOscCor,int probe,int nev,double fLO,
                       std::vector<double> time,std::vector<averageFixedProbeEvent_t> fxpr,
                       std::vector<trolleyAnaEvent_t> Data,std::vector<calibSwap_t> &Event){

   // find the mean field at the times specified in the time vector 
   calibSwap_t theEvent; 

   // do oscillation correction and obtain ALL data associated with toggle times in time vector 
   std::vector<double> trTime,trFreq,trFreq_cor; 
   int rc = CorrectOscillation_trly(probe,nev,time,fxpr,Data,trTime,trFreq,trFreq_cor);
   
   if(UseOscCor) std::cout << "[GetTRLYStatsAtTime]: USING OSCILLATION CORRECTION" << std::endl;

   // now need to average over each toggle 

   int n=0;
   int M = trTime.size();
   const int NT = time.size();
   double lastTime=0;
   double arg_freq=0,delta_t=0;
   double mean_freq=0,stdev_freq=0;
   double mean_temp=0,stdev_temp=0;
   double mean_x=0,stdev_x=0;
   double mean_y=0,stdev_y=0;
   double mean_z=0,stdev_z=0;
   std::vector<double> tt,freq,temp,x,y,z;
   for(int i=0;i<NT;i++){
      // find events that satisfy the timestamp for frequency and apply corrections 
      // WARNING: be careful to not accumulate events from previous swaps!  
      for(int j=0;j<M;j++){
	 if(trTime[j]>lastTime && trTime[j]<time[i]){
	    tt.push_back(trTime[j]); 
	    if(UseTempCor){
	       // FIXME: apply a temperature correction if necessary
	       // accounts for the trolley being at a temperature other than 25 deg c  
	       delta_t = (4.0E-9)*(temp[j]-25.0);
	    }
	    if(UseOscCor){
	       arg_freq = (fLO + trFreq_cor[j])/(1.-delta_t);  
	    }else{
	       arg_freq = (fLO + trFreq[j])/(1.-delta_t);  
	    }
	    freq.push_back(arg_freq); 
	 }
      }
      // std::cout << "SWAP " << i+1 << ": GATHERED " << freq.size() << " EVENTS" << std::endl;
      // gather all other variables
      rc = FilterSingle("temp",probe,nev,time[i],Data,temp);
      rc = FilterSingle("r"   ,probe,nev,time[i],Data,x);
      rc = FilterSingle("y"   ,probe,nev,time[i],Data,y);
      rc = FilterSingle("phi" ,probe,nev,time[i],Data,z);
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
      n = tt.size();
      lastTime = tt[n-1];
      tt.clear(); 
      freq.clear();
      temp.clear();
      x.clear();
      y.clear();
      z.clear();
   } 

   return 0;
}
//______________________________________________________________________________
int GetTRLYStats_sccToggle_old(bool useOscCor,int probe,int nev,std::vector<double> time,
                               std::vector<averageFixedProbeEvent_t> fxpr,std::vector<trolleyAnaEvent_t> Data,
                               std::vector<double> &TIME,std::vector<double> &MEAN,std::vector<double> &STDEV){

   // find the mean field at the times specified in the time vector 

   // first gather the FXPR data 
   std::vector<double> tt,ff;
   const int NFP = fxpr.size();
   for(int i=0;i<NFP;i++){
      tt.push_back( fxpr[i].time/1E+9 );
      ff.push_back( fxpr[i].freq );
   }

   const int N = time.size();
   int M=0,rc=0,lo=0,hi=0;
   double mean=0,stdev=0,f0=0,theFreq=0,delta_osc=0;
   std::vector<double> trt,freq;
   for(int i=0;i<N;i++){
      // std::cout << "Looking for time " << gm2fieldUtil::GetStringTimeStampFromUTC(time[i]) << std::endl;
      // find events 
      rc = FilterSingle("time",probe,nev,time[i],Data,trt);
      rc = FilterSingle("freq",probe,nev,time[i],Data,freq);
      if(useOscCor){
	 for(int j=0;j<M;j++){
	    gm2fieldUtil::Algorithm::BinarySearch<double>(tt,trt[j],lo,hi);
	    // take average over +/- 5 events in time 
	    rc = GetStats_vec(lo-5,hi+5,ff,theFreq,stdev);
	    if(j==0) f0 = theFreq;
	    delta_osc = theFreq - f0; 
	    freq[j] -= delta_osc; 
	 }
      }
      // now get mean of events 
      mean  = gm2fieldUtil::Math::GetMean<double>(freq);
      stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(freq);
      // store result
      TIME.push_back(time[i]); 
      MEAN.push_back(mean);
      STDEV.push_back(stdev);
      // set up for next time
      trt.clear(); 
      freq.clear();
   }

   return 0;
}
//______________________________________________________________________________
int GetTRLYStats_sccToggle(bool useOscCor,int probe,int nev,std::vector<double> time,
                           std::vector<averageFixedProbeEvent_t> fxpr,std::vector<trolleyAnaEvent_t> Data,
                           std::vector<double> &TIME,std::vector<double> &MEAN,std::vector<double> &STDEV){

   // find the mean field at the times specified in the time vector 

   // do oscillation correction and obtain ALL data associated with toggle times in time vector 
   std::vector<double> trTime,trFreq,trFreq_cor; 
   int rc = CorrectOscillation_trly(probe,nev,time,fxpr,Data,trTime,trFreq,trFreq_cor);

   const int NT = time.size();
   int M = trTime.size();
   if(useOscCor){ 
      std::cout << "[GetTRLYStats_sccToggle]: USING OSCILLATION CORRECTION" << std::endl;
      // for(int i=0;i<M;i++) std::cout << Form("Event %d, change = %.3lf Hz",i+1,trFreq[i]-trFreq_cor[i]) << std::endl;
   }   
       
   // now need to average over each toggle 

   int n=0;
   double lastTime=0,mean=0,stdev=0;
   std::vector<double> tt,ff; 
   for(int i=0;i<NT;i++){
      // find events that satisfy the timestamp
      // WARNING: must be careful to not accumulate events from previous toggles!  
      for(int j=0;j<M;j++){
	 if(trTime[j]>lastTime && trTime[j]<time[i]){
	    tt.push_back(trTime[j]);
	    if(useOscCor){
	       ff.push_back(trFreq_cor[j]);
            }else{
	       ff.push_back(trFreq[j]);
            }
	 }
      }
      // get stats
      mean  = gm2fieldUtil::Math::GetMean<double>(ff); 
      stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(ff);
      // store results  
      TIME.push_back(time[i]); 
      MEAN.push_back(mean); 
      STDEV.push_back(stdev);
      // set up for next time 
      n = tt.size();
      lastTime = tt[n-1]; 
      tt.clear();
      ff.clear();
   } 

   return 0;
}
