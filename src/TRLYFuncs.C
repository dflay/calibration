#include "../include/TRLYFuncs.h"
//______________________________________________________________________________
int GetTrolleyProbePositions(trolleyProbePosition_t &data){
    // this is in mm!
 
    double r1 = 17.5;
    double r2 = 35.;

    for(int i=0;i<17;i++){
       data.phi[i]  = 0; 
       data.dphi[i] = 0;
       data.dr[i]   = 0; 
       data.dy[i]   = 0;
    } 

    data.r[0] = 0.;
    data.y[0] = 0.;

    data.r[1] = 0.;
    data.y[1] = -r1;

    data.r[2] = r1;
    data.y[2] = 0.;

    data.r[3] = 0.;
    data.y[3] = r1;

    data.r[4] = -r1;
    data.y[4] = 0.;

    data.r[5] = r2*TMath::Cos(9*TMath::Pi()/6);
    data.y[5] = r2*TMath::Sin(9*TMath::Pi()/6);

    data.r[6] = r2*TMath::Cos(10*TMath::Pi()/6);
    data.y[6] = r2*TMath::Sin(10*TMath::Pi()/6);

    data.r[7] = r2*TMath::Cos(11*TMath::Pi()/6);
    data.y[7] = r2*TMath::Sin(11*TMath::Pi()/6);

    data.r[8] = r2*TMath::Cos(12*TMath::Pi()/6);
    data.y[8] = r2*TMath::Sin(12*TMath::Pi()/6);

    data.r[9] = r2*TMath::Cos(13*TMath::Pi()/6);
    data.y[9] = r2*TMath::Sin(13*TMath::Pi()/6);

    data.r[10] = r2*TMath::Cos(14*TMath::Pi()/6);
    data.y[10] = r2*TMath::Sin(14*TMath::Pi()/6);

    data.r[11] = r2*TMath::Cos(15*TMath::Pi()/6);
    data.y[11] = r2*TMath::Sin(15*TMath::Pi()/6);

    data.r[12] = r2*TMath::Cos(16*TMath::Pi()/6);
    data.y[12] = r2*TMath::Sin(16*TMath::Pi()/6);

    data.r[13] = r2*TMath::Cos(17*TMath::Pi()/6);
    data.y[13] = r2*TMath::Sin(17*TMath::Pi()/6);

    data.r[14] = r2*TMath::Cos(18*TMath::Pi()/6);
    data.y[14] = r2*TMath::Sin(18*TMath::Pi()/6);

    data.r[15] = r2*TMath::Cos(19*TMath::Pi()/6);
    data.y[15] = r2*TMath::Sin(19*TMath::Pi()/6);

    data.r[16] = r2*TMath::Cos(20*TMath::Pi()/6);
    data.y[16] = r2*TMath::Sin(20*TMath::Pi()/6);

    return 0;

}
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
int FilterSingle(int probe,int nev,double T,std::vector<trolleyAnaEvent_t> in,std::vector<double> &out){

   std::vector<double> tt;
   const int N = in.size();
   for(int i=0;i<N;i++) tt.push_back(in[i].time[probe]/1E+9);

   int lo=0,hi=0;
   gm2fieldUtil::Algorithm::BinarySearch(tt,T,lo,hi);
   int end   = lo-10;
   int start = end-nev;
   for(int i=start;i<end;i++){
      // std::cout << Form("%s: %.3lf",gm2fieldUtil::GetStringTimeStampFromUTC(in[i].time[probe]/1E+9).c_str(),in[i].freq[probe]) << std::endl; 
      out.push_back(in[i].freq[probe]);
   }
   // std::cout << "--------------------------" << std::endl;

   return 0;

}
//______________________________________________________________________________
int GetTRLYStatsAtTime(int probe,int nev,double fLO,std::vector<double> time,std::vector<trolleyAnaEvent_t> Data,
                       std::vector<double> &MEAN,std::vector<double> &STDEV){

   // find the mean field at the times specified in the time vector 

   const int N = time.size();
   int M=0,rc=0;
   double mean=0,stdev=0;
   std::vector<double> freq;
   for(int i=0;i<N;i++){
      // std::cout << "Looking for time " << gm2fieldUtil::GetStringTimeStampFromUTC(time[i]) << std::endl;
      // find events 
      rc = FilterSingle(probe,nev,time[i],Data,freq);
      // now get mean of events 
      mean  = gm2fieldUtil::Math::GetMean<double>(freq);
      stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(freq);
      // store result
      MEAN.push_back(mean+fLO);
      STDEV.push_back(stdev);
      // set up for next time 
      freq.clear();
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
      rc = FilterSingle(probe,nev,time[i],Data,freq);
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
