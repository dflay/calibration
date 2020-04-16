#include "../include/CustomAlgorithms.h"
//______________________________________________________________________________
int GetAvgSCCMagnitude(std::string type,std::vector<surfaceCoilEvent_t> data,double &mean,double &err){
   // Get surface coil maximum magnitude (excluding zeroes) 
   std::vector<double> x;
   double mag=0;
   double thr = 0.100; // cut threshold
   int M = 100;
   if(type.compare("azi")==0) M = 4;
   const int N = data.size(); 
   for(int i=0;i<N;i++){
      for(int j=0;j<M;j++){
	 if(type.compare("top")==0) mag = TMath::Abs(data[i].TopCurrents[j] );  
	 if(type.compare("bot")==0) mag = TMath::Abs(data[i].BotCurrents[j] ); 
	 if(type.compare("azi")==0) mag = TMath::Abs(data[i].AzCurrents[j]  );
	 if(mag>thr) x.push_back(mag); // cut out zeroes 
      }
   }
   mean = gm2fieldUtil::Math::GetMean<double>(x); 
   err  = gm2fieldUtil::Math::GetStandardDeviation<double>(x);
   return 0; 
}
//______________________________________________________________________________
int CheckShimGrad(int probe,std::vector<grad_meas_t> &data){
   // make sure we have non-zero entries for all values
   // zeroes are usually because we don't have a measurement 
   // first find the maximum
   char msg[200];
   const int N = data.size();
   double max          = data[0].grad;
   double max_err      = data[0].grad_err;
   double max_fxpr     = data[0].grad_fxpr;
   double max_fxpr_err = data[0].grad_fxpr_err;
   for(int i=0;i<N;i++){
      if( fabs(max)<=fabs(data[i].grad) ){
         max          = data[i].grad;
         max_err      = data[i].grad_err;
         max_fxpr     = data[i].grad_fxpr;
         max_fxpr_err = data[i].grad_fxpr_err;
      }
   }

   bool isBad = gm2fieldUtil::Math::IsInfOrNaN<double>(max);
   if(isBad){
      sprintf(msg,"[CustomAlgorithms::CheckShimGrad]: Warning for maximum gradients (probe %02d)! inf or NaN.  Using 1000 Hz...",probe);
      Logger::PrintMessage(Logger::kINFO,"default",msg,'a');
      max          = 1E+3;
      max_err      = 1E+3;
      max_fxpr     = 1E+3;
      max_fxpr_err = 1E+3;
   }

   // now check to see if we need to update a value 
   for(int i=0;i<N;i++){
      if( data[i].grad!=0 ){
         // things are good, do nothing
      }else{
         // invalid gradient.  Make an update using the worst gradient we have 
         sprintf(msg,"[CustomAlgorithms::CheckShimGrad]: Warning for probe %02d, axis %d! Using worst gradient of %.3lf +/- %.3lf Hz (ABA %.3lf +/- %.3lf Hz)",
                           probe,i,max,max,max_fxpr,max_fxpr);
         Logger::PrintMessage(Logger::kINFO,"default",msg,'a');
	 // use 100% uncertainty 
         data[i].grad          = max;
         data[i].grad_err      = TMath::Abs(max);
         data[i].grad_fxpr     = max_fxpr;
         data[i].grad_fxpr_err = TMath::Abs(max_fxpr);
      }
   }
   return 0;
}
//______________________________________________________________________________
int CopyTrolleyProbe(std::vector<trolleyAnaEvent_t> x,std::vector<trolleyAnaEvent_t> &y){
   // copy a trolley probe vector to a new one 
   trolleyAnaEvent_t data;
   const int N = x.size();
   if(N==0){
      std::cout << "[CustomAlgorithms::CopyTrolleyProbe]: No data!" << std::endl;
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
int CopyPlungingProbe(plungingProbeAnaEvent_t data,plungingProbeAnaEvent_t &y){
   // copy a plunging probe struct to a new one 
   int M = data.numTraces;
   y.numTraces = M;
   y.run       = data.run;
   y.midasRun  = data.midasRun;
   for(int i=0;i<M;i++){
      y.time[i]          = data.time[i];
      y.t2Time[i]        = data.t2Time[i];
      y.nzc[i]           = data.nzc[i];
      y.traceNumber[i]   = data.traceNumber[i];
      y.channelNumber[i] = data.channelNumber[i];
      y.freq_LO[i]       = data.freq_LO[i];
      y.freq_RF[i]       = data.freq_RF[i];
      y.freq[i]          = data.freq[i];
      y.freq_err[i]      = data.freq_err[i];
      y.r[i]             = data.r[i];
      y.y[i]             = data.y[i];
      y.phi[i]           = data.phi[i];
      y.temp[i]          = data.temp[i];
      y.temp_err[i]      = data.temp_err[i];
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
      CopyPlungingProbe(x[i],data);
      y.push_back(data);
   }

   return 0;
}
// //______________________________________________________________________________
// int CopyPlungingProbe_nmrAna_to_ppAna_single(int method,nmrAnaEvent_t data,plungingProbeAnaEvent_t &y){
//    // copy a plunging probe nmrAnaEvent_t vector to a plungingProbeAnaEvent_t struct
//    // note that the vector does NOT have a MIDAS run attached to it necessarily, 
//    // and that the index moves from the that of the VECTOR to that of the data MEMBERS
//    int M       = data.numTraces; 
//    y.run       = data.run;
//    y.numTraces = M;
//    y.midasRun  = 0;
//    for(int i=0;i<M;i++){
//       y.time[i]          = data[i].time;
//       y.t2Time[i]        = data[i].t2;
//       y.nzc[i]           = data[i].zc;
//       y.traceNumber[i]   = data[i].pulse;
//       y.channelNumber[i] = data[i].ch;
//       y.freq_LO[i]       = data[i].freq_LO;
//       y.freq_RF[i]       = data[i].freq_pi2;
//       y.freq[i]          = data[i].freq[method];
//       y.freq_err[i]      = 0;
//       y.r[i]             = 0;
//       y.y[i]             = 0;
//       y.phi[i]           = 0;
//       y.temp[i]          = data[i].temp;
//       y.temp_err[i]      = 0;
//    }
//    return 0;
// }
//______________________________________________________________________________
int CopyPlungingProbe_nmrAna_to_ppAna(int method,std::vector<nmrAnaEvent_t> data,plungingProbeAnaEvent_t &y){
   // copy a plunging probe nmrAnaEvent_t vector to a plungingProbeAnaEvent_t struct
   // note that the vector does NOT have a MIDAS run attached to it necessarily, 
   // and that the index moves from the that of the VECTOR to that of the data MEMBERS 
   int M       = data.size();
   y.run       = data[0].run;
   y.numTraces = M;
   y.midasRun  = 0;
   for(int i=0;i<M;i++){
      y.time[i]          = data[i].time;
      y.t2Time[i]        = data[i].t2;
      y.nzc[i]           = data[i].zc;
      y.traceNumber[i]   = data[i].pulse;
      y.channelNumber[i] = data[i].ch;
      y.freq_LO[i]       = data[i].freq_LO;
      y.freq_RF[i]       = data[i].freq_pi2;
      y.freq[i]          = data[i].freq[method];
      y.freq_err[i]      = 0;
      y.r[i]             = 0;
      y.y[i]             = 0;
      y.phi[i]           = 0;
      y.temp[i]          = data[i].temp;
      y.temp_err[i]      = 0;
   }
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
   int L=6;
   plungingProbeAnaEvent_t data;
   for(int i=0;i<N;i++){ // PP-DAQ run loop
      nmrDAQ_run = x[i].run; 
      // std::cout << "Looking for run " << nmrDAQ_run << std::endl;
      for(int j=0;j<M;j++){ // sub run loop 
         if(nmrDAQ_run==subRun[j]){
            CopyPlungingProbe(x[i],data);
            y.push_back(data);
	    // for(int k=0;k<L;k++){
	    //    std::cout << Form("[FilterPlungingProbeData]: FOUND: run %d, trace %d, freq = %.3lf Hz",
            //                   data.run,data.traceNumber[k],data.freq[k]) << std::endl;
	    // }
         }
      }
   }
   return 0;
}
//______________________________________________________________________________
int GetWeightedAverageStats(std::vector<double> x,std::vector<double> dx,double &mean,double &err,double &stdev){
   // get weighted mean, stdev of a data set x +/- dx
   double arg=0;
   const int N = x.size();
   std::vector<double> weight;
   for(int i=0;i<N;i++){
      if(dx[i]!=0){
	 arg = 1./( dx[i]*dx[i] );
      }else{
	 arg = 1.;
      }
      weight.push_back(arg);
   }

   int rc = gm2fieldUtil::Math::GetWeightedMean<double>(x,weight,mean,err);
   stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(x);

   return 0;
}
//______________________________________________________________________________
int FindGalilEvent(int probe,unsigned long long trlyTime,
                   std::vector<trolleyAnaEvent_t> trly,
                   std::vector<gm2field::galilTrolley_t> galil){
   // need to apply an offset to the galil data set 
   double sf = 2597./2658.; // scale factor for galil data
   double dt = trly[0].time[probe]/1E+9 - sf*galil[0].TimeStamp/1E+3;

   // now search through Galil data to find the right event that's 
   // within a few seconds on the trolley time
   int theEvent=-1; 
   double thr = 2.;
   double arg_tr=trlyTime/1E+9,arg_ga=0;
   const int N = galil.size();
   double DT   = (galil[N-1].TimeStamp-galil[0].TimeStamp)/1E+3;
   int eventRate = N/( (int)DT ); 
   for(int i=0;i<N;i++){
      arg_ga = sf*galil[i].TimeStamp/1E+3 + dt;
      if( TMath::Abs(arg_tr-arg_ga)<thr ){
         theEvent = i;
         break;
      }
   }
   return theEvent;
}
//______________________________________________________________________________
int FindTRLYStopTimes(int probe,double angle,std::vector<trolleyAnaEvent_t> trlyData,
                      std::vector<gm2field::galilTrolley_t> trlyGalil,
                      std::vector<double> &time){
   // find the time at which the trolley stops at the azimuthal location 
   // indicated by the angle variable

   int i=0,j=0,cntr=0;
   double theTime=0; 
   double z=0,vFish=0,vSig=0;
   double dphi = 0.2; 
   double delta = 30;
   double angle_min = angle - dphi; 
   double angle_max = angle + dphi;

   std::string timeStr; 

   const int N = trlyData.size(); 
   std::cout << Form("[FindTRLYStopTimes]: Scanning through %d trolley probe %02d events...",N,probe+1) << std::endl;
   do{ 
      // get trolley speeds 
      j       = FindGalilEvent(probe,trlyData[i].time[probe],trlyData,trlyGalil);
      vFish   = trlyGalil[j].Velocities[0]; 
      vSig    = trlyGalil[j].Velocities[1]; 
      // get trolley position 
      z  = trlyData[i].phi[probe]; 
      if( (z>angle_min)&&(z<angle_max)  
         &&(TMath::Abs(vFish)==0)&&(TMath::Abs(vSig)==0) ){
	 cntr++;
         if(cntr>1){
	    delta = 55.;
         }else{
	    delta = 30.;
         }
	 // the trolley is stopped at the target location 
         theTime = trlyData[i].time[probe]/1E+9 + delta;
         timeStr = gm2fieldUtil::GetStringTimeStampFromUTC( (unsigned long)theTime ); 
	 time.push_back(theTime); 
	 std::cout << Form("Found stop at %s (%d) for event number %d",timeStr.c_str(),(int)theTime,i) << std::endl;
	 i += 300; 
      }else{
	 i++;
      }
   }while(i<N); 
 
   std::cout << "--> Done." << std::endl;

   return 0;
}
// //______________________________________________________________________________
// int FindTransitionTimes_old(int type,int gradType,double thr,double delta,std::vector<gm2field::surfaceCoils_t> data,
//                         std::vector<double> &timeOff,std::vector<double> &timeOn){
//    // find the transition times of turning off and on the surface coils 
//    // use the sum of the coils to determine if they're on or off 
//    // type: 0 (bottom), 1 (top), -1 (azi)
//    // thr: threshold above which we consider the SCC to be on 
//    // delta: how much time to delay marking the transition (for downstream analysis) 
//    int k=0; 
//    int rc=0,i=0,cntr=0,M=4,currentErr=0;
//    if(type==0||type==1) M = 100;
//    if( TMath::Abs(type)>1 ) return -1; 
//    int sign_top=0,sign_bot=0,sign_azi=0;
//    double currentThr=2.0;
//    double bot=0,top=0,azi=0;  
//    double BOT=0,TOP=0,AZI=0;  
//    double change=0,theTime=0,theTime_prev=0;
//    double diff=0,sum=0,sum2=0,sum_prev=0;
//    double SUM=0,SUM_PREV=0; 
//    const int NEV = data.size();
//    // std::cout << "Scanning " << NEV << " events..." << std::endl;
//    do{ 
//       for(int j=0;j<M;j++){
// 	 bot = data[i].BotCurrents[j]; 
// 	 top = data[i].TopCurrents[j]; 
// 	 azi = data[i].AzCurrents[0];
//          BOT = TMath::Abs(bot); 
//          TOP = TMath::Abs(top); 
//          AZI = TMath::Abs(azi); 
//          sign_bot = bot/BOT; 
//          sign_top = top/TOP; 
//          sign_azi = azi/AZI; 
//          if(gradType==0){
// 	    // x axis, sign should be the same; also require valid top current values 
// 	    if(sign_bot!=sign_top&&TOP>0.100&&BOT>0.100) diff = 1;
//          }else if(gradType==1){
// 	    // y axis, sign should be different; also require valid top current values 
// 	    if(sign_bot==sign_top&&TOP>0.100&&BOT>0.100) diff = 1;
//          }else{
// 	    diff=0;
//          }
// 	 // azi can't be on while others are on 
// 	 if(AZI>0.100 && TOP>0.100) diff = 2;
// 	 if(AZI>0.100 && BOT>0.100) diff = 2;
// 	 // sum based on type 
// 	 if(type==0)  sum += bot; 
// 	 if(type==1)  sum += top; 
// 	 if(type==-1) sum  = azi; // these coils are set symmetrically, so we don't want to sum! 
//          if(TMath::Abs(sum)>currentThr||diff!=0){
// 	    currentErr = 1;
//             theTime = data[i].TopTime[0]/1E+9; 
//             std::cout << "ERROR CODE " << diff << " AT " << gm2fieldUtil::GetStringTimeStampFromUTC(theTime) << std::endl;
// 	    break;
//          } 
//       }
//       if(currentErr==1){
// 	 // not correct currents, continue on
// 	 currentErr = 0;
// 	 // sum_prev = 0;
// 	 sum      = 0;
// 	 diff     = 0;
// 	 i++;
// 	 continue; 
//       }
//       // get the time of the event 
//       if(type==0)  theTime = data[i].BotTime[0]/1E+9; 
//       if(type==1)  theTime = data[i].TopTime[0]/1E+9; 
//       if(type==-1) theTime = data[i].TopTime[0]/1E+9; // we don't have an Azi time 
//       std::cout << "time = " << gm2fieldUtil::GetStringTimeStampFromUTC(theTime) 
//                 << " sum = " << sum << " time_prev = " << theTime_prev << " sum_prev = " << sum_prev << std::endl;
//       // now compute the change
//       change = sum - sum_prev;
//       // std::cout << gm2fieldUtil::GetStringTimeStampFromUTC(theTime) << " " << sum_prev << " " << sum << " CURRENT CHANGE = " << change << std::endl;
//       // if(k%5==0){
//       //    // hit the 5th second, take mean 
//       //    SUM /= 5.; 
//       //    change   = (SUM-SUM_PREV);
//       //    SUM_PREV = SUM;
//       //    std::cout << k << " " << gm2fieldUtil::GetStringTimeStampFromUTC(theTime) 
//       //              << " " << SUM_PREV << " " << SUM << " CURRENT CHANGE = " << change << std::endl;
//       //    SUM = 0;
//       //    k   = 0; 
//       // }else{
//       //    // keep accumulating 
//       //    SUM += sum;
//       // }
//       // k++;
//       // categorize as an on or off transition  
//       if( TMath::Abs(change)>thr ){
// 	 // found a transition
//          std::cout << "--> FOUND A TRANSITION at: " << gm2fieldUtil::GetStringTimeStampFromUTC(theTime) << " change = " << change << std::endl; 
// 	 cntr++;
// 	 // if( (type>=0&&cntr>0) || type==-1){  // FIXME: Why do we need the counter to be greater than 2 before??
// 	 //    // now is it a time off or time on? 
// 	 //    if( TMath::Abs(sum)<thr ){   
// 	 //       timeOff.push_back(theTime+delta); 
// 	 //       // std::cout << "sum = " << sum; 
// 	 //       // std::cout << " OFF transition " << timeOff.size() << ": " << gm2fieldUtil::GetStringTimeStampFromUTC(theTime+delta) << std::endl; 
// 	 //    }else{
// 	 //       timeOn.push_back(theTime+delta); 
// 	 //       // std::cout << "sum = " << sum; 
//          //       // std::cout << " ON transition " << timeOn.size() << ": " << gm2fieldUtil::GetStringTimeStampFromUTC(theTime+delta) << std::endl; 
// 	 //    }
// 	 // }
//          if(change>0) timeOn.push_back(theTime+delta);  
//          if(change<0) timeOff.push_back(theTime+delta);  
// 	 // i += 20;   // jump ahead 
//       }else{
// 	 // no transition 
//       }
//       // set up for next event 
//       sum_prev     = sum;
//       theTime_prev = theTime;
//       sum          = 0;
//       diff         = 0;
//       i++;
//    }while(i<NEV); 
// 
//    // std::cout << "--> Done." << std::endl;
// 
//    const int Non  = timeOn.size();
//    const int Noff = timeOff.size();
// 
//    // no transitions?! return an error 
//    if(Non==0 || Noff==0){
//       std::cout << "PROBLEM! No transitions!" << std::endl;
//       return -1; 
//    }
// 
//    if(Non!=Noff){
//       std::cout << "PROBLEM! Unequal on/off transitions!" << std::endl;
//       std::cout << "on: " << Non << std::endl;
//       for(int i=0;i<Non;i++){
// 	 std::cout << "on: "  << gm2fieldUtil::GetStringTimeStampFromUTC( timeOn[i] ) << std::endl;
//       }
//       std::cout << " off: " << Noff << std::endl;
//       for(int i=0;i<Noff;i++){
// 	 std::cout << "on: "  << gm2fieldUtil::GetStringTimeStampFromUTC( timeOff[i] ) << std::endl;
//       }
//       return -2;
//    }
// 
//    // determine if we start analysis with SCC on or off
//    // rc = 1, SCC on to start 
//    if(timeOn[0]<timeOff[0])  rc = 1;
// 
//    std::cout << Form("Number of SCC transition times: off = %d, on = %d",Noff,Non) << std::endl;
//    for(int i=0;i<Non;i++){
//       std::cout << "on: "  << gm2fieldUtil::GetStringTimeStampFromUTC( timeOn[i] ) << " "  
//                 << "off: " << gm2fieldUtil::GetStringTimeStampFromUTC( timeOff[i] ) << std::endl;
//    } 
// 
//    return rc;  
// }
//______________________________________________________________________________
int FindTransitionTimes(int type,int gradType,double thr,double delta,std::vector<surfaceCoilEvent_t> data,
                        std::vector<double> &timeOff,std::vector<double> &timeOn,int coilNum){
   // find the transition times of turning off and on the surface coils 
   // use the sum of the coils to determine if they're on or off 
   // type: 0 (bottom), 1 (top), -1 (azi)
   // thr: threshold above which we consider the SCC to be on 
   // delta: how much time to delay marking the transition (for downstream analysis) 
   bool debug = false; 
   int k=0; 
   int rc=0,i=0,cntr=0,M=4,currentErr=0;
   if(type==0||type==1) M = 100;
   if( TMath::Abs(type)>1 ) return -1; 
   int sign_top=0,sign_bot=0,sign_azi=0;
   double currentThr=2.0;
   double bot=0,top=0,azi=0;  
   double BOT=0,TOP=0,AZI=0;  
   double change=0,theTime=0,theTime_prev=0;
   double diff=0,sum=0,sum2=0,sum_prev=0;
   double SUM=0,SUM_PREV=0; 
   const int NEV = data.size();
   // std::cout << "Scanning " << NEV << " events..." << std::endl;
   do{ 
      for(int j=0;j<M;j++){
         if(coilNum>0){
	    bot = data[i].BotCurrents[coilNum]; 
	    top = data[i].TopCurrents[coilNum]; 
         }else{
	    bot = data[i].BotCurrents[j]; 
	    top = data[i].TopCurrents[j]; 
         }
	 azi = data[i].AzCurrents[0];
         BOT = TMath::Abs(bot); 
         TOP = TMath::Abs(top); 
         AZI = TMath::Abs(azi); 
         sign_bot = bot/BOT; 
         sign_top = top/TOP; 
         sign_azi = azi/AZI; 
         if(gradType==0){
	    // x axis, sign should be the same; also require valid top current values 
	    if(sign_bot!=sign_top&&TOP>0.100&&BOT>0.100) diff = 1;
         }else if(gradType==1){
	    // y axis, sign should be different; also require valid top current values 
	    if(sign_bot==sign_top&&TOP>0.100&&BOT>0.100) diff = 1;
         }else{
	    diff=0;
         }
	 // azi can't be on while others are on 
	 if(AZI>0.100 && TOP>0.100) diff = 2;
	 if(AZI>0.100 && BOT>0.100) diff = 2;
	 // sum based on type 
	 if(type==0)  sum += bot; 
	 if(type==1)  sum += top; 
	 if(type==-1) sum  = azi; // these coils are set symmetrically, so we don't want to sum! 
         if(TMath::Abs(sum)>currentThr||diff!=0){
	    currentErr = 1;
            theTime = data[i].TopTime[0]/1E+9; 
            if(debug) std::cout << "ERROR CODE " << diff << " AT " << gm2fieldUtil::GetStringTimeStampFromUTC(theTime) << std::endl;
	    break;
         } 
      }
      if(currentErr==1){
	 // not correct currents, continue on
	 currentErr = 0;
	 // sum_prev = 0;
	 sum      = 0;
	 diff     = 0;
	 i++;
	 continue; 
      }
      // get the time of the event 
      if(type==0)  theTime = data[i].BotTime[0]/1E+9; 
      if(type==1)  theTime = data[i].TopTime[0]/1E+9; 
      if(type==-1) theTime = data[i].TopTime[0]/1E+9; // we don't have an Azi time 
      if(debug){
	 std::cout << "time = " << gm2fieldUtil::GetStringTimeStampFromUTC(theTime) 
                << " sum = " << sum << " time_prev = " << theTime_prev << " sum_prev = " << sum_prev << std::endl;
      }
      // now compute the change
      change = sum - sum_prev;
      // categorize as an on or off transition  
      if( TMath::Abs(change)>thr ){
	 // found a transition
         if(debug) std::cout << "--> FOUND A TRANSITION at: " << gm2fieldUtil::GetStringTimeStampFromUTC(theTime) << " change = " << change << std::endl; 
	 cntr++;
         if(change>0) timeOn.push_back(theTime+delta);  
         if(change<0) timeOff.push_back(theTime+delta);  
	 // i += 20;   // jump ahead 
      }else{
	 // no transition 
      }
      // set up for next event 
      sum_prev     = sum;
      theTime_prev = theTime;
      sum          = 0;
      diff         = 0;
      i++;
   }while(i<NEV); 

   // std::cout << "--> Done." << std::endl;

   // remove duplicates
   // convert to integer first since that's a better filter for duplicates 
   std::vector<int> ton,toff; 
   int N = timeOn.size(); 
   for(int i=0;i<N;i++) ton.push_back( (int)timeOn[i] ); 
   N = timeOff.size();  
   for(int i=0;i<N;i++) toff.push_back( (int)timeOff[i] );  
   gm2fieldUtil::Algorithm::SortedRemoveDuplicates<int>(ton);  
   gm2fieldUtil::Algorithm::SortedRemoveDuplicates<int>(toff);  

   timeOn.clear(); 
   timeOff.clear(); 
   for(int i=0;i<N;i++) timeOff.push_back( (double)toff[i] );
   N = ton.size();
   for(int i=0;i<N;i++) timeOn.push_back( (double)ton[i] );

   const int Non  = timeOn.size();
   const int Noff = timeOff.size();

   // no transitions?! return an error 
   if(Non==0 || Noff==0){
      std::cout << "PROBLEM! No transitions!" << std::endl;
      return -1; 
   }

   if(Non!=Noff){
      std::cout << "PROBLEM! Unequal on/off transitions!" << std::endl;
      std::cout << "on: " << Non << std::endl;
      for(int i=0;i<Non;i++){
	 std::cout << "on: "  << gm2fieldUtil::GetStringTimeStampFromUTC( timeOn[i] ) << std::endl;
      }
      std::cout << " off: " << Noff << std::endl;
      for(int i=0;i<Noff;i++){
	 std::cout << "on: "  << gm2fieldUtil::GetStringTimeStampFromUTC( timeOff[i] ) << std::endl;
      }
      return -2;
   }

   // determine if we start analysis with SCC on or off
   // rc = 1, SCC on to start 
   if(timeOn[0]<timeOff[0])  rc = 1;

   std::cout << Form("Number of SCC transition times: off = %d, on = %d",Noff,Non) << std::endl;
   for(int i=0;i<Non;i++){
      std::cout << "on: "  << gm2fieldUtil::GetStringTimeStampFromUTC( timeOn[i] ) << " "  
                << "off: " << gm2fieldUtil::GetStringTimeStampFromUTC( timeOff[i] ) << std::endl;
   } 

   return rc;  
}
//______________________________________________________________________________
int SwapEntries(int i,std::vector<double> &x,std::vector<double> &y){
   double val_x = x[i]; 
   double val_y = y[i]; 
   y[i] = val_x; 
   x[i] = val_y; 
   return 0; 
}
//______________________________________________________________________________
int CheckDifference(std::vector<double> &x,std::vector<double> &dx){
   // check for bad entries
   int rc=-1;
   int cntr_pos=0,cntr_neg=0;
   std::vector<double> sign;
   const int N = x.size();
   for(int i=0;i<N;i++){
      if(x[i]>=0){
	 cntr_pos++;
	 sign.push_back(1); 
      }else{ 
	 cntr_neg++;
	 sign.push_back(-1); 
      } 
   }

   if(cntr_pos==0||cntr_neg==0){
      // no changes needed, exit
      return rc;
   }else{
      std::cout << "[CheckDifference]: Will adjust computed difference!" << std::endl;
   }
 
   // determine which sign is dominant 
   // if positive, domSign = 1, -1 otherwise
   int domSign = -1;  
   if(cntr_pos>cntr_neg) domSign = 1; 

   double val=0;
   for(int i=0;i<N;i++){
      if(domSign==1){
	 // OK, we're expecting a positive value.  check entries, flip sign if necessary 
	 if( x[i]<0 ){
	    val = (-1.)*x[i];
	    x.erase(x.begin()+i);
	    x.insert(x.begin()+i,val);
	    std::cout << "[CheckDifference]: Change made at entry " << i << std::endl;
	    rc = i;
	 }
      }else{
	 // OK, we're expecting a negative value.  check entries, flip sign if necessary 
	 if( x[i]>0 ){
	    val = (-1.)*x[i];
	    x.erase(x.begin()+i);
	    x.insert(x.begin()+i,val);
	    std::cout << "[CheckDifference]: Change made at entry " << i << std::endl;
	    rc = i;
	 }
      }
   }

   return rc; 

}
//______________________________________________________________________________
int GetDifference(std::vector<double> A,std::vector<double> A_err,
                  std::vector<double> B,std::vector<double> B_err,
                  std::vector<double> &diff,std::vector<double> &diff_err){
   double arg=0,arg_err=0;
   const int N = B.size();
   for(int i=0;i<N;i++){
      arg     = A[i] - B[i];
      arg_err = TMath::Sqrt(A_err[i]*A_err[i] + B_err[i]*B_err[i]);
      std::cout << Form("Trial %d: A = %.3lf +/- %.3lf Hz, B = %.3lf +/- %.3lf Hz, diff(A-B) = %.3lf +/- %.3lf Hz",
                        i,A[i],A_err[i],B[i],B_err[i],arg,arg_err) << std::endl;
      diff.push_back(arg);
      diff_err.push_back(arg_err);
   }
   return 0;
}
//______________________________________________________________________________
int GetDifference_ABA(bool useTimeWeight,
                      std::vector<double> sccTime ,std::vector<double> scc ,std::vector<double> scc_err,
                      std::vector<double> bareTime,std::vector<double> bare,std::vector<double> bare_err,
                      std::vector<double> &diff_aba,std::vector<double> &diff_aba_err){

   // WARNING: This assumes that the bare measurement comes first!
   // FIXME? THIS IS WRONG! 

   double w=0,w_prev=0;
   double dt_tot=0;
   double diff=0,diff_err=0;
   double diff_prev=0,diff_prev_err=0;
   double arg=0,arg_err=0;
   const int N = bare.size();
   for(int i=1;i<N;i++){
      dt_tot        = bareTime[i]-bareTime[i-1]; 
      w             = TMath::Abs(sccTime[i-1]-bareTime[i-1])/dt_tot;
      w_prev        = TMath::Abs(sccTime[i-1]-bareTime[i])/dt_tot;
      // first compute the difference 
      diff_prev     = scc[i-1] - bare[i-1];
      diff_prev_err = TMath::Sqrt(scc_err[i-1]*scc_err[i-1] + bare_err[i-1]*bare_err[i-1]);
      diff          = scc[i-1] - bare[i];
      diff_err      = TMath::Sqrt(scc_err[i-1]*scc_err[i-1] + bare_err[i]*bare_err[i]);
      // now get the ABA difference 
      if(useTimeWeight){
	 arg     = w*diff     + w_prev*diff_prev;
	 // arg_err = w*diff_err + w_prev*diff_prev_err;
	 arg_err = TMath::Sqrt( scc_err[i-1]*scc_err[i-1] + w*w*bare_err[i]*bare_err[i] + w_prev*w_prev*bare_err[i-1]*bare_err[i-1] );
      }else{
	 arg     = 0.5*(diff     + diff_prev);
	 // arg_err = 0.5*(diff_err + diff_prev_err);
	 arg_err = TMath::Sqrt( scc_err[i-1]*scc_err[i-1] + 0.25*bare_err[i]*bare_err[i] + 0.25*bare_err[i-1]*bare_err[i-1] );
      }
      // std::cout << Form("Trial %d: bare1 = %.3lf +/- %.3lf Hz, scc = %.3lf +/- %.3lf Hz, bare2 = %.3lf +/- %.3lf Hz, diff = %.3lf +/- %.3lf Hz",
      //                   i,bare[i-1],bare_err[i-1],scc[i-1],scc_err[i-1],bare[i],bare_err[i],arg,arg_err) << std::endl;
      // store result 
      diff_aba.push_back(arg);
      diff_aba_err.push_back(arg_err);
   }
   return 0;
}
//______________________________________________________________________________
int GetDifference_ABA_sccFirst(bool useTimeWeight,
                               std::vector<double> sccTime ,std::vector<double> scc ,std::vector<double> scc_err,
                               std::vector<double> bareTime,std::vector<double> bare,std::vector<double> bare_err,
                               std::vector<double> &diff_aba,std::vector<double> &diff_aba_err){

   // WARNING: This assumes that the scc measurement comes first!
   double w=0,w_prev=0,dt_tot=0;
   double diff=0,diff_err=0;
   double diff_prev=0,diff_prev_err=0;
   double arg=0,arg_err=0;
   const int N = bare.size();
   for(int i=1;i<N;i++){
      dt_tot        = sccTime[i]-sccTime[i-1]; 
      w             = TMath::Abs(sccTime[i-1]-bareTime[i-1])/dt_tot;
      w_prev        = TMath::Abs(sccTime[i]-bareTime[i-1])/dt_tot;
      // first compute the difference 
      diff_prev     = scc[i-1] - bare[i-1];
      diff_prev_err = TMath::Sqrt(scc_err[i-1]*scc_err[i-1] + bare_err[i-1]*bare_err[i-1]);
      diff          = scc[i] - bare[i-1];
      diff_err      = TMath::Sqrt(scc_err[i]*scc_err[i] + bare_err[i-1]*bare_err[i-1]);
      // now get the ABA difference
      if(useTimeWeight){
	 arg     = w*diff     + w_prev*diff_prev;
	 // arg_err = w*diff_err + w_prev*diff_prev_err;
	 arg_err = TMath::Sqrt( bare_err[i-1]*bare_err[i-1] + w*w*scc_err[i]*scc_err[i] + w_prev*w_prev*scc_err[i-1]*scc_err[i-1] );
      }else{ 
	 arg     = 0.5*(diff     + diff_prev);
	 // arg_err = 0.5*(diff_err + diff_prev_err);
	 arg_err = TMath::Sqrt( bare_err[i-1]*bare_err[i-1] + 0.25*scc_err[i]*scc_err[i] + 0.25*scc_err[i-1]*scc_err[i-1] );
      }
       std::cout << Form("Trial %d: scc1 = %.3lf +/- %.3lf Hz, bare = %.3lf +/- %.3lf Hz, scc2 = %.3lf +/- %.3lf Hz, diff = %.3lf +/- %.3lf Hz",
                         i,scc[i-1],scc_err[i-1],bare[i-1],bare_err[i-1],scc[i],scc_err[i],arg,arg_err) << std::endl;
      // store result 
      diff_aba.push_back(arg);
      diff_aba_err.push_back(arg_err);
   }
   return 0;
}
//______________________________________________________________________________
int GetDifference_ABA_final(bool useTimeWeight,
                            std::vector<double> A_time,std::vector<double> A,std::vector<double> A_err,
                            std::vector<double> B_time,std::vector<double> B,std::vector<double> B_err,
                            std::vector<double> &diff_aba,std::vector<double> &diff_aba_err){
   // WARNING: This assumes that the A measurement comes first!
   std::cout << "[GetDifference_ABA_final]: Computing difference..." << std::endl;
   double w=0,w_prev=0,dt_tot=0;
   double diff=0,diff_prev=0;
   double arg=0,arg_err=0,arg_err_sq=0;
   const int N = B.size();
   for(int i=1;i<N;i++){
      // set up time weights 
      dt_tot    = A_time[i] - A_time[i-1]; 
      w         = TMath::Abs(A_time[i-1]-B_time[i-1])/dt_tot;
      w_prev    = TMath::Abs(A_time[i]-B_time[i-1])/dt_tot;
      // first compute the difference 
      diff_prev = A[i-1] - B[i-1];
      diff      = A[i]   - B[i-1];
      // now get the ABA difference
      if(useTimeWeight){
	 arg        = w*diff + w_prev*diff_prev;
	 arg_err_sq = B_err[i-1]*B_err[i-1] + w*w*A_err[i]*A_err[i] + w_prev*w_prev*A_err[i-1]*A_err[i-1];
         // account for correlations in measurements 
	 // arg_err_sq = B_err[i-1]*B_err[i-1] + w*w*A_err[i]*A_err[i] + w_prev*w_prev*A_err[i-1]*A_err[i-1] 
         //            + 2.*w*w_prev*A_err[i]*A_err[i-1] - 2.*w*w_prev*A_err[i]*B_err[i-1] - 2.*w*w_prev*A_err[i-1]*B_err[i-1];
	 arg_err    = TMath::Sqrt(arg_err_sq);
      }else{ 
	 arg     = 0.5*(diff + diff_prev);
	 arg_err = TMath::Sqrt( B_err[i-1]*B_err[i-1] + 0.25*A_err[i]*A_err[i] + 0.25*A_err[i-1]*A_err[i-1] );
      }
      std::cout << Form("Trial %02d: ",i); 
      std::cout << Form("A1 = %.3lf +/- %.3lf Hz, ",A[i-1],A_err[i-1]); 
      std::cout << Form("B = %.3lf +/- %.3lf Hz, " ,B[i-1],B_err[i-1]);
      std::cout << Form("A2 = %.3lf +/- %.3lf Hz, ",A[i]  ,A_err[i]  );
      std::cout << Form("diff = %.3lf +/- %.3lf Hz",arg   ,arg_err   ) << std::endl;
      // store result 
      diff_aba.push_back(arg);
      diff_aba_err.push_back(arg_err);
   }
   std::cout << "[GetDifference_ABA_final]: Done." << std::endl;
   return 0;
}
//______________________________________________________________________________
TGraph *RemoveTrend(TGraph *g1,TF1 *func,double &mean,double &stdev){
   // subtract off a trend in the data   

   const int N1 = g1->GetN();
   double *x1 = g1->GetX();

   // int NPTS = 50;
   // double step = (max-min)/( (double)NPTS ); 
   std::vector<double> X,Y; 
   double ix=0,diff=0;
   for(int i=0;i<N1;i++){
      ix   = x1[i]; 
      diff = g1->Eval(ix) - func->Eval(ix); 
      X.push_back(ix);  
      Y.push_back(diff);
   }

   mean  = gm2fieldUtil::Math::GetMean<double>(Y); 
   stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(Y); 

   TGraph *g = gm2fieldUtil::Graph::GetTGraph(X,Y);
   return g;

}
//______________________________________________________________________________
TGraph *GetDiffPlot(TGraph *g1,TGraph *g2){
   // compute difference of two TGraphs  

   const int N1 = g1->GetN();
   const int N2 = g2->GetN();
   double *x1 = g1->GetX();
   double *x2 = g2->GetX();

   // choose which graph to use as the baseline for number of points 
   int ref  = 1; 
   int NPTS = N1;
   if(N2<N1){
      ref = 2;
      NPTS = N2;
   } 

   std::vector<double> X,Y; 
   double ix=0,diff=0;
   for(int i=0;i<NPTS;i++){
      ix   = x1[i]; // min + ( (double)i )*step;
      if(ref==2) ix = x2[i]; 
      diff = g2->Eval(ix) - g1->Eval(ix); 
      // std::cout << ix << " " << diff << std::endl; 
      X.push_back(ix);  
      Y.push_back(diff);
   }

   TGraph *g = gm2fieldUtil::Graph::GetTGraph(X,Y);
   return g;

}
//______________________________________________________________________________
TGraph *GetDiffPlot(TGraphErrors *g1,TGraphErrors *g2){
   // compute difference of two TGraphs  

   const int N1 = g1->GetN();
   const int N2 = g2->GetN();
   double *x1 = g1->GetX();
   double *x2 = g2->GetX();

   // choose which graph to use as the baseline for number of points 
   int ref  = 1; 
   int NPTS = N1;
   if(N2<N1){
      ref = 2;
      NPTS = N2;
   } 

   std::vector<double> X,Y; 
   double ix=0,diff=0;
   for(int i=0;i<NPTS;i++){
      ix   = x1[i]; // min + ( (double)i )*step;
      if(ref==2) ix = x2[i]; 
      diff = g2->Eval(ix) - g1->Eval(ix); 
      // std::cout << ix << " " << diff << std::endl; 
      X.push_back(ix);  
      Y.push_back(diff);
   }

   TGraph *g = gm2fieldUtil::Graph::GetTGraph(X,Y);
   return g;

}
//______________________________________________________________________________
int CorrectPPForDriftDuringMeasurement(int method,TF1 *fxprFit,
                                       plungingProbeAnaEvent_t ppEvent,plungingProbeAnaEvent_t &ppEventCor){
   // correct the PP data for drift as monitored by the fixed probes
   // uses a fit to the fixed probes 

   // copy over all data into the new vector 
   CopyPlungingProbe(ppEvent,ppEventCor);

   double drift=0,fxpr0=0,mean_fxpr=0,stdev_fxpr=0,freq_cor=0;
   fxpr0 = fxprFit->Eval(ppEvent.time[0]/1E+9); 
   // std::cout << "fxpr0 = " << Form("%.3lf",fxpr0) << std::endl; 

   std::cout << "CORRECTING RUN " << ppEvent.run << std::endl; 
   const int M = ppEvent.numTraces;
   for(int i=0;i<M;i++){
      mean_fxpr = fxprFit->Eval(ppEvent.time[i]/1E+9); 
      // compute drift: take difference in current fxpr freq and subtract starting fxpr freq  
      drift    = mean_fxpr - fxpr0;
      freq_cor = ppEvent.freq[i] - drift;
      // update PP frequency 
      ppEventCor.freq[i] = freq_cor;
      std::cout << "time = "          << Form("%s",gm2fieldUtil::GetStringTimeStampFromUTC(ppEvent.time[i]/1E+9 ).c_str() )  << " "
                << "temp = "          << Form("%.3lf deg C",ppEvent.temp[i])  << " "
                // << "pos (mm) = "      << Form("(%.3lf,%.3lf,%.3lf)",ppEvent.r[i],ppEvent.y[i],ppEvent.phi[i]) << " "
                << "pp freq = "       << Form("%.3lf Hz",ppEvent.freq[i]) << " " 
                << "fxpr = "          << Form("%.3lf Hz",mean_fxpr)          << " " 
                << "drift = "         << Form("%.3lf Hz",drift)              << " " 
                << "pp freq (cor) = " << Form("%.3lf Hz",freq_cor)       << std::endl; 
   }
   std::cout << "--------------------" << std::endl;

   return 0;
}
// //______________________________________________________________________________
// int CorrectPPForDriftDuringMeasurement(std::vector<fixedProbeEvent_t> fxprData,
//                                        plungingProbeAnaEvent_t ppEvent,plungingProbeAnaEvent_t &ppEventCor,bool isScan){
//    // correct the PP data for drift as monitored by the fixed probes
//    // uses a fit to the fixed probes 
// 
//    // copy over all data into the new vector 
//    CopyPlungingProbe(ppEvent,ppEventCor);
// 
//    double drift=0,fxpr0=0,mean_fxpr=0,stdev_fxpr=0,freq_cor=0;
//    GetAverageFXPR(ppEvent.time[0]/1E+9,fxprData,fxpr0,stdev_fxpr); 
// 
//    // std::cout << "CORRECTING RUN " << ppEvent.run << std::endl; 
//    const int M = ppEvent.numTraces;
//    for(int i=0;i<M;i++){
//       GetAverageFXPR(ppEvent.time[i]/1E+9,fxprData,mean_fxpr,stdev_fxpr); 
//       // compute drift: take difference in current fxpr freq and subtract starting fxpr freq  
//       drift    = mean_fxpr - fxpr0;
//       freq_cor = ppEvent.freq[i] - drift;
//       // update PP frequency 
//       ppEventCor.freq[i] = freq_cor;
//       // std::cout << "time = "          << Form("%s",gm2fieldUtil::GetStringTimeStampFromUTC(ppEvent.time[i]/1E+9 ).c_str() )  << " "
//       //           << "temp = "          << Form("%.3lf deg C",ppEvent.temp[i])  << " "
//       //           // << "pos (mm) = "      << Form("(%.3lf,%.3lf,%.3lf)",ppEvent.r[i],ppEvent.y[i],ppEvent.phi[i]) << " "
//       //           << "pp freq = "       << Form("%.3lf Hz",ppEvent.freq[i]) << " " 
//       //           << "fxpr = "          << Form("%.3lf Hz",mean_fxpr)          << " " 
//       //           << "drift = "         << Form("%.3lf Hz",drift)              << " " 
//       //           << "pp freq (cor) = " << Form("%.3lf Hz",freq_cor)       << std::endl; 
//    }
//    // std::cout << "--------------------" << std::endl;
// 
//    return 0;
// }
// //______________________________________________________________________________
// int CorrectPPForDriftDuringMeasurement(std::vector<fixedProbeEvent_t> fxprData,
//                                        std::vector<plungingProbeAnaEvent_t> ppEvent,
//                                        std::vector<plungingProbeAnaEvent_t> &ppEventCor,
//                                        bool isScan){
//    // correct the PP data for drift as monitored by the fixed probes
//    // uses the average field seen by fixed probes 
// 
//    // copy over all data into the new vector 
//    CopyPlungingProbe(ppEvent,ppEventCor);
// 
//    double drift=0,fxpr0=0,mean_fxpr=0,stdev_fxpr=0,freq_cor=0;
//    GetAverageFXPR(ppEvent[0].time[0]/1E+9,fxprData,fxpr0,stdev_fxpr); 
//    // std::cout << "fxpr0 = " << Form("%.3lf",fxpr0) << std::endl; 
// 
//    int M=0;
//    const int NPP = ppEvent.size();
//    for(int i=0;i<NPP;i++){
//       // std::cout << "RUN " << Form("%d",ppEvent[i].run)  << std::endl; 
//       M = ppEvent[i].numTraces;
//       for(int j=0;j<M;j++){ 
//          if(isScan && j==0) GetAverageFXPR(ppEvent[i].time[j]/1E+9,fxprData,fxpr0,stdev_fxpr); 
// 	 GetAverageFXPR(ppEvent[i].time[j]/1E+9,fxprData,mean_fxpr,stdev_fxpr); 
//          // compute drift: take difference in current fxpr freq and subtract starting fxpr freq  
//          drift    = mean_fxpr - fxpr0;
//          freq_cor = ppEvent[i].freq[j] - drift;
//          // update PP frequency 
//          ppEventCor[i].freq[j] = freq_cor;
//          // std::cout << "time = "          << Form("%s",gm2fieldUtil::GetStringTimeStampFromUTC(ppEvent[i].time[j]/1E+9 ).c_str() )  << " "
//          //           << "temp = "          << Form("%.3lf deg C",ppEvent[i].temp[j])  << " "
//          //           // << "pos (mm) = "      << Form("(%.3lf,%.3lf,%.3lf)",ppEvent[i].r[j],ppEvent[i].y[j],ppEvent[i].phi[j]) << " "
//          //           << "pp freq = "       << Form("%.3lf Hz",ppEvent[i].freq[j]) << " " 
//          //           << "fxpr = "          << Form("%.3lf Hz",mean_fxpr)          << " " 
//          //           << "drift = "         << Form("%.3lf Hz",drift)              << " " 
//          //           << "pp freq (cor) = " << Form("%.3lf Hz",freq_cor)       << std::endl; 
//       }
//       // std::cout << "--------------------" << std::endl;
//    }
// 
//    return 0;
// }
//______________________________________________________________________________
int CorrectPPForDriftDuringMeasurement(int method,TF1 *fxprFit,
                                       std::vector<plungingProbeAnaEvent_t> ppEvent,std::vector<plungingProbeAnaEvent_t> &ppEventCor){
   // correct the PP data for drift as monitored by the fixed probes
   // uses a fit to the fixed probes 

   // copy over all data into the new vector 
   CopyPlungingProbe(ppEvent,ppEventCor);

   double drift=0,fxpr0=0,mean_fxpr=0,stdev_fxpr=0,freq_cor=0;
   fxpr0 = fxprFit->Eval(ppEvent[0].time[0]/1E+9); 
   // std::cout << "fxpr0 = " << Form("%.3lf",fxpr0) << std::endl; 

   int M=0;
   const int NPP = ppEvent.size();
   for(int i=0;i<NPP;i++){
      // std::cout << "RUN " << Form("%d",ppEvent[i].run)  << std::endl; 
      M = ppEvent[i].numTraces;
      for(int j=0;j<M;j++){
         mean_fxpr = fxprFit->Eval(ppEvent[i].time[j]/1E+9); 
         // compute drift: take difference in current fxpr freq and subtract starting fxpr freq  
         drift    = mean_fxpr - fxpr0;
         freq_cor = ppEvent[i].freq[j] - drift;
         // update PP frequency 
         ppEventCor[i].freq[j] = freq_cor;
         // std::cout << "time = "          << Form("%s",gm2fieldUtil::GetStringTimeStampFromUTC(ppEvent[i].time[j]/1E+9 ).c_str() )  << " "
         //           << "temp = "          << Form("%.3lf deg C",ppEvent[i].temp[j])  << " "
         //           // << "pos (mm) = "      << Form("(%.3lf,%.3lf,%.3lf)",ppEvent[i].r[j],ppEvent[i].y[j],ppEvent[i].phi[j]) << " "
         //           << "pp freq = "       << Form("%.3lf Hz",ppEvent[i].freq[j]) << " " 
         //           << "fxpr = "          << Form("%.3lf Hz",mean_fxpr)          << " " 
         //           << "drift = "         << Form("%.3lf Hz",drift)              << " " 
         //           << "pp freq (cor) = " << Form("%.3lf Hz",freq_cor)       << std::endl; 
      }
      // std::cout << "--------------------" << std::endl;
   }

   return 0;
}
// //______________________________________________________________________________
// int CorrectPPForDriftDuringMeasurement(int method,std::vector<int> fxprList,std::vector<gm2field::fixedProbeFrequency_t> fxprData,
//                                        plungingProbeAnaEvent_t ppEvent,plungingProbeAnaEvent_t &ppEventCor){
//    // correct the PP data for drift as monitored by the fixed probes 
// 
//    // copy over all data into the new vector 
//    CopyPlungingProbe(ppEvent,ppEventCor);
// 
//    double drift=0,fxpr0=0,mean_fxpr=0,stdev_fxpr=0,freq_cor=0;
//    GetAverageFXPR(method,ppEvent.time[0],fxprList,fxprData,fxpr0,stdev_fxpr); 
//    // std::cout << "fxpr0 = " << Form("%.3lf",fxpr0) << std::endl; 
// 
//    const int M = ppEvent.numTraces;
//    for(int i=0;i<M;i++){
//       GetAverageFXPR(method,ppEvent.time[i],fxprList,fxprData,mean_fxpr,stdev_fxpr);
//       // compute drift: take difference in current fxpr freq and subtract starting fxpr freq  
//       drift    = mean_fxpr - fxpr0;
//       freq_cor = ppEvent.freq[i] - drift;
//       // update PP frequency 
//       ppEventCor.freq[i] = freq_cor;
//       // std::cout << "time = "          << Form("%s",gm2fieldUtil::GetStringTimeStampFromUTC(ppEvent[i].time[j]/1E+9 ).c_str() )  << " "
//       //           << "temp = "          << Form("%.3lf deg C",ppEvent[i].temp[j])  << " "
//       //           // << "pos (mm) = "      << Form("(%.3lf,%.3lf,%.3lf)",ppEvent[i].r[j],ppEvent[i].y[j],ppEvent[i].phi[j]) << " "
//       //           << "pp freq = "       << Form("%.3lf Hz",ppEvent[i].freq[j]) << " " 
//       //           << "fxpr = "          << Form("%.3lf Hz",mean_fxpr)          << " " 
//       //           << "drift = "         << Form("%.3lf Hz",drift)              << " " 
//       //           << "pp freq (cor) = " << Form("%.3lf Hz",freq_cor)       << std::endl; 
//    }
// 
//    return 0;
// }
// //______________________________________________________________________________
// int CorrectPPForDriftDuringMeasurement(int method,std::vector<int> fxprList,std::vector<gm2field::fixedProbeFrequency_t> fxprData,
//                                        std::vector<plungingProbeAnaEvent_t> ppEvent,std::vector<plungingProbeAnaEvent_t> &ppEventCor){
//    // correct the PP data for drift as monitored by the fixed probes 
// 
//    // copy over all data into the new vector 
//    CopyPlungingProbe(ppEvent,ppEventCor);
// 
//    double drift=0,fxpr0=0,mean_fxpr=0,stdev_fxpr=0,freq_cor=0;
//    GetAverageFXPR(method,ppEvent[0].time[0],fxprList,fxprData,fxpr0,stdev_fxpr); 
//    // std::cout << "fxpr0 = " << Form("%.3lf",fxpr0) << std::endl; 
// 
//    int M=0;
//    const int NPP = ppEvent.size();
//    for(int i=0;i<NPP;i++){
//       // std::cout << "RUN " << Form("%d",ppEvent[i].run)  << std::endl; 
//       M = ppEvent[i].numTraces;
//       for(int j=0;j<M;j++){
//          GetAverageFXPR(method,ppEvent[i].time[j],fxprList,fxprData,mean_fxpr,stdev_fxpr);
//          // compute drift: take difference in current fxpr freq and subtract starting fxpr freq  
//          drift    = mean_fxpr - fxpr0;
//          freq_cor = ppEvent[i].freq[j] - drift;
//          // update PP frequency 
//          ppEventCor[i].freq[j] = freq_cor;
//          // std::cout << "time = "          << Form("%s",gm2fieldUtil::GetStringTimeStampFromUTC(ppEvent[i].time[j]/1E+9 ).c_str() )  << " "
//          //           << "temp = "          << Form("%.3lf deg C",ppEvent[i].temp[j])  << " "
//          //           // << "pos (mm) = "      << Form("(%.3lf,%.3lf,%.3lf)",ppEvent[i].r[j],ppEvent[i].y[j],ppEvent[i].phi[j]) << " "
//          //           << "pp freq = "       << Form("%.3lf Hz",ppEvent[i].freq[j]) << " " 
//          //           << "fxpr = "          << Form("%.3lf Hz",mean_fxpr)          << " " 
//          //           << "drift = "         << Form("%.3lf Hz",drift)              << " " 
//          //           << "pp freq (cor) = " << Form("%.3lf Hz",freq_cor)       << std::endl; 
//       }
//       // std::cout << "--------------------" << std::endl;
//    }
// 
//    return 0;
// }
//______________________________________________________________________________
int CorrectPPForDriftDuringMeasurement(int method,std::vector<int> trlyList,std::vector<trolleyAnaEvent_t> trlyData,
                                       std::vector<plungingProbeAnaEvent_t> ppEvent,std::vector<plungingProbeAnaEvent_t> &ppEventCor){
   // correct the PP data for drift as monitored by the fixed probes 

   // copy over all data into the new vector 
   CopyPlungingProbe(ppEvent,ppEventCor);

   double drift=0,trly0=0,mean_trly=0,stdev_trly=0,freq_cor=0;
   
   // get trolley field at the START of the run      
   GetAverageTRLY(ppEvent[0].time[0],trlyList,trlyData,trly0,stdev_trly);

   int M=0;
   const int NPP = ppEvent.size();
   for(int i=0;i<NPP;i++){
      // std::cout << "RUN " << Form("%d",ppEvent[i].run)  << std::endl; 
      M = ppEvent[i].numTraces;
      for(int j=0;j<M;j++){
         GetAverageTRLY(ppEvent[i].time[j],trlyList,trlyData,mean_trly,stdev_trly);
         // event zero is the point in time we reference to
         // recall that index i is referring to a different LOCATION necessarily; index j is the first measurement at location i. 
         // if(j==0) trly0 = mean_trly;
         // compute drift: take difference in current fxpr freq and subtract starting fxpr freq  
         drift    = mean_trly - trly0;
         freq_cor = ppEvent[i].freq[j] - drift;
         // update PP frequency 
         ppEventCor[i].freq[j] = freq_cor;
         // std::cout << "time = "          << Form("%s",gm2fieldUtil::GetStringTimeStampFromUTC(ppEvent[i].time[j]/1E+9 ).c_str() )  << " "
         //           << "temp = "          << Form("%.3lf deg C",ppEvent[i].temp[j])  << " "
         //           << "pos (mm) = "      << Form("(%.3lf,%.3lf,%.3lf)",ppEvent[i].r[j],ppEvent[i].y[j],ppEvent[i].phi[j]) << " "
         //           << "pp freq = "       << Form("%.3lf Hz",ppEvent[i].freq[j]) << " " 
         //           << "drift = "         << Form("%.3lf Hz",drift)              << " " 
         //           << "pp freq (cor) = " << Form("%.3lf Hz",freq_cor)       << std::endl; 
      }
      // std::cout << "--------------------" << std::endl;
   }

   return 0;
}
// //______________________________________________________________________________
// int CorrectPPForDriftDuringMeasurementAlt(int method,std::vector<int> fxprList,std::vector<gm2field::fixedProbeFrequency_t> fxprData,
//                                           std::vector<plungingProbeAnaEvent_t> ppEvent,std::vector<plungingProbeAnaEvent_t> &ppEventCor,
//                                           std::vector<drift_t> &drift){
//    // correct the PP data for drift as monitored by the fixed probes 
//    double intercept=0,slope=0,r=0;
//    double mean_fxpr=0,stdev_fxpr=0;
//    std::vector<double> TIME,FREQ; // for fixed probe drift calc 
// 
//    drift_t drift_calc;
// 
//    double fxpr0=0;
// 
//    int M=0;
//    const int NPP = ppEvent.size();
//    for(int i=0;i<NPP;i++){
//       std::cout << "RUN " << Form("%d",ppEvent[i].run)  << std::endl;
//       M = ppEvent[i].numTraces;
//       for(int j=0;j<M;j++){
//             GetAverageFXPR(method,ppEvent[i].time[j],fxprList,fxprData,mean_fxpr,stdev_fxpr);
//             TIME.push_back(ppEvent[i].time[j]/1E+9);
//             FREQ.push_back(mean_fxpr);
//             gm2fieldUtil::Math::LeastSquaresFitting(TIME,FREQ,intercept,slope,r);
//             // std::cout << "time = "      << Form("%s",gm2fieldUtil::GetStringTimeStampFromUTC(ppEvent[i].time[j]/1E+9 ).c_str() )  << " "
//             //           << "temp = "      << Form("%.3lf deg C",ppEvent[i].temp[j])  << " "
//             //           << "pos (mm) = "  << Form("(%.3lf,%.3lf,%.3lf)",ppEvent[i].r[j],ppEvent[i].y[j],ppEvent[i].phi[j]) << " "
//             //           << "pp freq = "   << Form("%.3lf Hz",ppEvent[i].freq[j]) << " "
//             //           << "fxpr freq = " << Form("%.3lf Hz",mean_fxpr)          << std::endl;
//       }
//       drift_calc.intercept = intercept;
//       drift_calc.slope     = slope;
//       drift.push_back(drift_calc);
//       // std::cout << "--------------------" << std::endl;
//       // clean up 
//       TIME.clear();
//       FREQ.clear();
//    }
// 
//    // copy over all data into the new vector 
//    CopyPlungingProbe(ppEvent,ppEventCor);
//    double t0=0;
//    double freq_cor=0,mean_before=0,mean_after=0,stdev_before=0,stdev_after=0;
//    std::vector<double> fb,fa;
//    // apply correction for drift 
//    for(int i=0;i<NPP;i++){
//       M = ppEvent[i].numTraces;
//       std::cout << "RUN "     << ppEvent[i].run << std::endl;
//       std::cout << "DRIFT = " << Form("%.3lf Hz/sec",drift[i].slope) << std::endl;
//       for(int j=0;j<M;j++){
//          if(j==0) t0 = ppEvent[i].time[0];
//          freq_cor = ppEvent[i].freq[j] - drift[i].slope*(ppEvent[i].time[j]-t0)/1E+9;
//          std::cout << "time = "      << Form("%s",gm2fieldUtil::GetStringTimeStampFromUTC(ppEvent[i].time[j]/1E+9 ).c_str() )  << " "
//                    << "temp = "      << Form("%.3lf deg C",ppEvent[i].temp[j])  << " "
//                    << "pos (mm) = "  << Form("(%.3lf,%.3lf,%.3lf)",ppEvent[i].r[j],ppEvent[i].y[j],ppEvent[i].phi[j]) << " "
//                    << "pp freq = "   << Form("%.3lf Hz",ppEvent[i].freq[j]) << " " 
//                    << "pp freq cor = " << Form("%.3lf Hz",freq_cor)         << std::endl;
//          fb.push_back(ppEvent[i].freq[j]);
//          fa.push_back(freq_cor); 
//          // update the corrected event 
//          ppEventCor[i].freq[j] = freq_cor;
//       }  
//       mean_before  = gm2fieldUtil::Math::GetMean<double>(fb);
//       stdev_before = gm2fieldUtil::Math::GetStandardDeviation<double>(fb);
//       mean_after   = gm2fieldUtil::Math::GetMean<double>(fa);
//       stdev_after  = gm2fieldUtil::Math::GetStandardDeviation<double>(fa);
//       std::cout << "mean before = " << Form("%.3lf +/- %.3lf Hz",mean_before,stdev_before) << " "
//                 << "mean after = " << Form("%.3lf +/- %.3lf Hz",mean_after,stdev_after) << std::endl;
//       fb.clear();
//       fa.clear();
//       std::cout << "--------------------" << std::endl;
//    }  
//    return 0;
// }
//______________________________________________________________________________
int CorrectTRLYForDriftDuringMeasurement(int method,TF1 *fxprFit,
                                         std::vector<trolleyAnaEvent_t> Event,std::vector<trolleyAnaEvent_t> &EventCor){
   // correct the trolley data for drift as monitored by the fixed probes

   // copy over all data into the new vector 
   CopyTrolleyProbe(Event,EventCor);

   double drift=0,fxpr0=0,mean_fxpr=0,stdev_fxpr=0,freq_cor=0;
   std::cout << "Using fixed probe data" << std::endl;

   // for drift correction, reference the first event on probe 0 -- that's the start of the run 
   fxpr0 = fxprFit->Eval(Event[0].time[0]/1E+9); 
  
   const int NTR = Event.size();
   std::cout << "Processing " << NTR << " events" << std::endl;
   for(int i=0;i<NUM_TRLY;i++){  // trolley probe
      for(int j=0;j<NTR;j++){    // event number
            mean_fxpr = fxprFit->Eval(Event[j].time[i]/1E+9); 
            // compute drift: take difference in current fxpr freq and subtract starting fxpr freq  
            drift    = mean_fxpr - fxpr0;
            freq_cor = Event[j].freq[i] - drift;
            // update TRLY frequency 
            EventCor[j].freq[i] = freq_cor;
            // std::cout << "time = "            << Form("%s",gm2fieldUtil::GetStringTimeStampFromUTC(Event[j].time[i]/1E+9 ).c_str() )  << " "
            //           << "temp = "            << Form("%.3lf deg C",Event[j].temp[i]) << " "
            //           << "pos (mm,mm,deg) = " << Form("(%.3lf,%.3lf,%.3lf)",Event[j].r[i],Event[j].y[i],Event[j].phi[i]) << " "
            //           << "trly freq = "       << Form("%.3lf Hz",Event[j].freq[i])    << " "
            //           << "drift = "           << Form("%.3lf Hz",drift)               << " "
            //           << "trly freq (cor) = " << Form("%.3lf Hz",freq_cor)            << std::endl;
      }
      std::cout << "Finished probe " << i+1 <<" --------------------" << std::endl;
   }

   return 0;
}
// //______________________________________________________________________________
// int CorrectTRLYForDriftDuringMeasurement(std::vector<fixedProbeEvent_t> fxprData,
//                                          std::vector<trolleyAnaEvent_t> Event,std::vector<trolleyAnaEvent_t> &EventCor){
//    // correct the trolley data for drift as monitored by the fixed probes
// 
//    // copy over all data into the new vector 
//    CopyTrolleyProbe(Event,EventCor);
// 
//    double drift=0,fxpr0=0,mean_fxpr=0,stdev_fxpr=0,freq_cor=0;
//    
//    std::cout << "Using fixed probe data" << std::endl;
// 
//    // for drift correction, reference the first event on probe 0 -- that's the start of the run 
//    GetAverageFXPR(Event[0].time[0]/1E+9,fxprData,fxpr0,stdev_fxpr);
//   
//    const int NTR = Event.size();
//    std::cout << "Processing " << NTR << " events" << std::endl;
//    for(int i=0;i<NUM_TRLY;i++){  // trolley probe
//       for(int j=0;j<NTR;j++){    // event number
//             GetAverageFXPR(Event[j].time[i]/1E+9,fxprData,mean_fxpr,stdev_fxpr);
//             // compute drift: take difference in current fxpr freq and subtract starting fxpr freq  
//             drift    = mean_fxpr - fxpr0;
//             freq_cor = Event[j].freq[i] - drift;
//             EventCor[j].freq[i] = freq_cor;      // update TRLY frequency 
//             // std::cout << "time = "            << Form("%s",gm2fieldUtil::GetStringTimeStampFromUTC(Event[j].time[i]/1E+9 ).c_str() )  << " "
//             //           << "temp = "            << Form("%.3lf deg C",Event[j].temp[i]) << " "
//             //           // << "pos (mm,mm,deg) = " << Form("(%.3lf,%.3lf,%.3lf)",Event[j].r[i],Event[j].y[i],Event[j].phi[i]) << " "
//             //           << "trly freq = "       << Form("%.3lf Hz",Event[j].freq[i])    << " "
//             //           << "drift = "           << Form("%.3lf Hz",drift)               << " "
//             //           << "trly freq (cor) = " << Form("%.3lf Hz",freq_cor)            << std::endl;
//       }
//       std::cout << "Finished probe " << i+1 <<" --------------------" << std::endl;
//    }
// 
//    return 0;
// } 
// //______________________________________________________________________________
// int CorrectTRLYForDriftDuringMeasurement(int method,std::vector<int> fxprList,std::vector<gm2field::fixedProbeFrequency_t> fxprData,
//                                          std::vector<trolleyAnaEvent_t> Event,std::vector<trolleyAnaEvent_t> &EventCor){
//    // correct the trolley data for drift as monitored by the fixed probes
// 
//    // copy over all data into the new vector 
//    CopyTrolleyProbe(Event,EventCor);
// 
//    double drift=0,fxpr0=0,mean_fxpr=0,stdev_fxpr=0,freq_cor=0;
//    
//    std::cout << "Using fixed probe data" << std::endl;
// 
//    // for drift correction, reference the first event on probe 0 -- that's the start of the run 
//    GetAverageFXPR(method,Event[0].time[0],fxprList,fxprData,fxpr0,stdev_fxpr);
//   
//    const int NTR = Event.size();
//    std::cout << "Processing " << NTR << " events" << std::endl;
//    for(int i=0;i<NUM_TRLY;i++){  // trolley probe
//       for(int j=0;j<NTR;j++){    // event number
//             GetAverageFXPR(method,Event[j].time[i],fxprList,fxprData,mean_fxpr,stdev_fxpr);
//             // compute drift: take difference in current fxpr freq and subtract starting fxpr freq  
//             drift    = mean_fxpr - fxpr0;
//             freq_cor = Event[j].freq[i] - drift;
//             EventCor[j].freq[i] = freq_cor;      // update TRLY frequency 
//             // std::cout << "time = "            << Form("%s",gm2fieldUtil::GetStringTimeStampFromUTC(Event[j].time[i]/1E+9 ).c_str() )  << " "
//             //           << "temp = "            << Form("%.3lf deg C",Event[j].temp[i]) << " "
//             //           // << "pos (mm,mm,deg) = " << Form("(%.3lf,%.3lf,%.3lf)",Event[j].r[i],Event[j].y[i],Event[j].phi[i]) << " "
//             //           << "trly freq = "       << Form("%.3lf Hz",Event[j].freq[i])    << " "
//             //           << "drift = "           << Form("%.3lf Hz",drift)               << " "
//             //           << "trly freq (cor) = " << Form("%.3lf Hz",freq_cor)            << std::endl;
//       }
//       std::cout << "Finished probe " << i+1 <<" --------------------" << std::endl;
//    }
// 
//    return 0;
// }
//______________________________________________________________________________
int CorrectTRLYForDriftDuringMeasurement(int method,std::vector<int> trlyList,std::vector<trolleyAnaEvent_t> trlyData,
                                         std::vector<trolleyAnaEvent_t> Event,std::vector<trolleyAnaEvent_t> &EventCor){
   // correct the trolley data for drift as monitored by the trolley probes (!) 

   // copy over all data into the new vector 
   CopyTrolleyProbe(Event,EventCor);

   std::cout << "Using trolley data" << std::endl;

   int rc=0;
   double drift=0,trly0=0,mean_trly=0,stdev_trly=0,freq_cor=0;
            
   // for drift correction, reference the first event on probe 0 -- that's the start of the run 
   // rc = GetAverageTRLY(Event[0].time[0],trlyList,trlyData,trly0,stdev_trly);

   const int NTR = Event.size();
   std::cout << NTR << std::endl; 
   for(int i=0;i<NUM_TRLY;i++){  // trolley probe
      for(int j=0;j<NTR;j++){    // event number
	 rc = GetAverageTRLY(Event[j].time[i],trlyList,trlyData,mean_trly,stdev_trly);
	 if(rc!=0){
	    std::cout << "ERROR" << std::endl; 
	    return 1;
	 }
	 if(j==0) trly0 = mean_trly;  // event zero is the point in time we reference to 
	 // compute drift: take difference in current fxpr freq and subtract starting fxpr freq  
	 drift    = mean_trly - trly0;
	 freq_cor = Event[j].freq[i] - drift;
	 // update PP frequency 
	 EventCor[j].freq[i] = freq_cor;
	 // std::cout << "time = "            << Form("%s",gm2fieldUtil::GetStringTimeStampFromUTC(Event[j].time[i]/1E+9 ).c_str() )  << " "
	 //           << "temp = "            << Form("%.3lf deg C",Event[j].temp[i]) << " "
	 //           << "pos (mm,mm,deg) = " << Form("(%.3lf,%.3lf,%.3lf)",Event[j].r[i],Event[j].y[i],Event[j].phi[i]) << " "
	 //           << "trly freq = "       << Form("%.3lf Hz",Event[j].freq[i])    << " "
	 //           << "drift = "           << Form("%.3lf Hz",drift)               << " "
	 //           << "trly freq (cor) = " << Form("%.3lf Hz",freq_cor)            << std::endl;
      }
      std::cout << "Finished probe " << i+1 <<" --------------------" << std::endl;
   }

   return 0;
}
//______________________________________________________________________________
int CorrectTRLYForDriftAcrossRuns(unsigned long long tStart,unsigned long long tStop,
                                  TGraph *gDrift,std::vector<trolleyAnaEvent_t> &trlyData,
                                  double &drift,double &drift_err){

   // apply a correction to trolley data due to drift in changing from a previous trolley run 
   // to a new one, presumably where the SCC was turned on.  We use the gDrift TGraph to do this,
   // which has two fixed-probe data sets with bare field runs bookending the trolley run to be 
   // corrected. The times tStart and tStop correspond to the end of the previous and start of 
   // the current run times of the trolley, respectively.

   drift     = gDrift->Eval(tStop/1E+9) - gDrift->Eval(tStart/1E+9); // drift correction in Hz 
   drift_err = fabs(drift); 

   const int NEvents = trlyData.size();
   if(NEvents==0){
      std::cout << "[CorrectTRLYForDriftAcrossRuns]: No data!" << std::endl;
      return 1;
   }

   std::cout << Form("[CorrectTRLYForDriftAcrossRuns]: Applying a drift correction of %.3lf Hz",drift) << std::endl;

   for(int i=0;i<NEvents;i++){
      for(int j=0;j<NUM_TRLY;j++) trlyData[i].freq[j] -= drift; 
   }  

   return 0;
}
//______________________________________________________________________________
int CalculateAveragePP(std::vector<plungingProbeAnaEvent_t> ppData,double &B,double &B_err){
   // compute the average PP data
   // sum over all events  
   int M=0;
   const int N = ppData.size();  // number of runs 
   double arg=0;
   std::vector<double> x; 
   for(int i=0;i<N;i++){
      M = ppData[i].numTraces;
      for(int j=0;j<M;j++){
         arg = ppData[i].freq_LO[i] + ppData[i].freq[j];
	 x.push_back(arg); 
      }
   }

   // assign to output container 
   B     = gm2fieldUtil::Math::GetMean<double>(x);
   B_err = gm2fieldUtil::Math::GetStandardDeviation<double>(x);  // FIXME: does this make sense?   

   return 0;
}
// //______________________________________________________________________________
// int CalculateAveragePP(int method,std::vector<int> fxprList,std::vector<gm2field::fixedProbeFrequency_t> fxprData,
//                        std::vector<plungingProbeAnaEvent_t> ppData,double &B,double &B_err){
//    // compute the average PP data
//    // sum over all events  
//    int M=0;
//    const int N = ppData.size();  // number of runs 
//    double arg=0;
//    std::vector<double> x; 
//    for(int i=0;i<N;i++){
//       M = ppData[i].numTraces;
//       for(int j=0;j<M;j++){
//          arg = ppData[i].freq_LO[i] + ppData[i].freq[j];
// 	 x.push_back(arg); 
//       }
//    }
// 
//    // assign to output container 
//    B     = gm2fieldUtil::Math::GetMean<double>(x);
//    B_err = gm2fieldUtil::Math::GetStandardDeviation<double>(x);  // FIXME: does this make sense?   
// 
//    return 0;
// }
//______________________________________________________________________________
int CalculateAveragePP(bool isDriftCor,int method,std::vector<int> trlyList,std::vector<trolleyAnaEvent_t> trlyData,
                       std::vector<plungingProbeAnaEvent_t> ppData,double &B,double &B_err){
   // compute the average PP data 
   int M=0;
   const int N = ppData.size();  // number of runs 
   double mean=0,stdev=0,drift=0;
   double f=0,f0=0,stdev_trly=0;
   unsigned long long time = ppData[0].time[0]; 
   GetAverageTRLY(time,trlyList,trlyData,f0,stdev_trly);

   double arg=0;
   std::vector<double> x,FREQ,ERR; 
   for(int i=0;i<N;i++){
      M = ppData[i].numTraces;
      for(int j=0;j<M;j++){
         arg = ppData[i].freq_LO[i] + ppData[i].freq[j];
	 x.push_back(arg); 
      }
      // stats for a RUN 
      mean  = gm2fieldUtil::Math::GetMean<double>(x); 
      stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(x);
      // drift run-to-run 
      time = ppData[i].time[0]; // the marker for the time is the first event of the run 
      GetAverageTRLY(time,trlyList,trlyData,f,stdev_trly); 
      drift = f - f0;
      if(isDriftCor) mean -= drift; 
      // fill vectors 
      FREQ.push_back(mean); 
      ERR.push_back(stdev);  
      // set up for next run
      x.clear(); 
   }

   // assign to output container 
   B     = gm2fieldUtil::Math::GetMean<double>(FREQ);
   B_err = gm2fieldUtil::Math::GetStandardDeviation<double>(FREQ);  // FIXME: does this make sense?   

   return 0;
}
//______________________________________________________________________________
int CalculateTRLYAvg_Stationary(int probeNumber,std::vector<trolleyAnaEvent_t> Event,
                                double &mean,double &stdev){
   // get average of a given probe over all events 
   double arg=0;
   std::vector<double> x;
   const int N = Event.size(); 
   for(int i=0;i<N;i++){
      // get frequency of the event 
      arg = 61.74E+6 + Event[i].freq[probeNumber];  // the LO is added here 
      x.push_back(arg); 
   }
   // calculate the mean and standard deviation 
   mean  = gm2fieldUtil::Math::GetMean<double>(x);
   stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(x);
   return 0; 
}
//______________________________________________________________________________
int GetAverageTRLY(unsigned long long time,std::vector<int> probe,std::vector<trolleyAnaEvent_t> trlyData,double &mean,double &stdev){
   int eventNum=0;
   const int N = probe.size();
   if(N==0){
      std::cout << "[GetAverageTRLY]: No data!" << std::endl;
      return 1;
   }
   // gather the events we care about 
   std::vector<double> freq;
   for(int i=0;i<N;i++){
      eventNum = FindTrolleyEvent(time,trlyData,probe[i]);
      freq.push_back(trlyData[eventNum].freq[probe[i]]);
      // std::cout << Form("probe %d, event %d: %.3lf Hz",probe[i],eventNum,trlyData[eventNum].freq[probe[i]]) << std::endl;
   }
   // calculate the mean and standard deviation 
   mean  = gm2fieldUtil::Math::GetMean<double>(freq);
   stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(freq);
   // std::cout << "MEAN = " << mean << " +/- " << stdev << std::endl; 
   return 0; 
}
//______________________________________________________________________________       
int FindTrolleyEvent(unsigned long long time, std::vector<trolleyAnaEvent_t> trly, int TP){
   //Find first ok event (want something that is >20e9)
   //First trolley event is weird
   //Should be within first 20                                                          
   int startEvent = 0;
   for(int i=0;i<20;i++){
      if(trly[i].time[TP] > 20e9){
	 startEvent = i;
	 break;
      }
   }

   //Look for when the time is just smaller than a trolley time. Then compare between the 2 times around it  
   int highEvent = 0;
   int size = trly.size();
   for(int i=startEvent;i<size;i++){
      if(trly[i].time[TP] > time){
	 highEvent = i;
	 break;
      }
   }

   int Tevent = 0;
   if(fabs((Long64_t)(time - trly[highEvent].time[TP])) < fabs((Long64_t)(time - trly[highEvent-1].time[TP]))){
      Tevent = highEvent;
   }
   else Tevent = highEvent-1;

   return Tevent;
}
//______________________________________________________________________________       
int FindTrolleyEvent(unsigned long long time, std::vector<averageTrolleyAnaEvent_t> trly){
   //Find first ok event (want something that is >20e9)
   //First trolley event is weird
   //Should be within first 20                                                          
   int startEvent = 0;
   for(int i=0;i<20;i++){
      if(trly[i].time > 20e9){
	 startEvent = i;
	 break;
      }
   }

   //Look for when the time is just smaller than a trolley time. Then compare between the 2 times around it  
   int highEvent = 0;
   int size = trly.size();
   for(int i=startEvent;i<size;i++){
      if(trly[i].time > time){
	 highEvent = i;
	 break;
      }
   }

   int Tevent = 0;
   if(fabs((Long64_t)(time - trly[highEvent].time)) < fabs((Long64_t)(time - trly[highEvent-1].time))){
      Tevent = highEvent;
   }
   else Tevent = highEvent-1;

   return Tevent;
}
// //______________________________________________________________________________
// TGraph *GetDriftTGraph(int method,std::vector<int> driftRun,std::vector<int> fxprList,
//                        std::vector<double> &stats){
// 
//    const int ND = driftRun.size(); 
//    int rc=0,M=0,startIndex=0,stopIndex=0;
//    unsigned long long aTime;
//    std::vector<unsigned long long> runStart,runStop;
//    std::vector<gm2field::fixedProbeFrequency_t> fxprData;
// 
//    for(int i=0;i<ND;i++){
//       std::cout << "Getting drift run " << driftRun[i] << std::endl;
//       rc = gm2fieldUtil::RootHelper::GetFPFrequencies(driftRun[i],fxprData);
//       if(rc!=0){
//          std::cout << "No data!" << std::endl;
//          return NULL;
//       }
//       // finished getting a run, find its time 
//       M = fxprData.size();  // number of events 
//       stopIndex = M - 1;
//       aTime = fxprData[startIndex].GpsTimeStamp[fxprList[0]]; // take zeroth probe as the start point  
//       runStart.push_back(aTime);
//       aTime = fxprData[stopIndex].GpsTimeStamp[fxprList[0]];  // take zeroth probe as the start point  
//       runStop.push_back(aTime);
//       // set up for next run 
//       startIndex = stopIndex + 1;
//    }
//    std::cout << "--> Done." << std::endl;
// 
//    // store the data start and stop times of the interpolation  
//    // timeBound.push_back(runStop[0] ); 
//    // timeBound.push_back(runStart[1]);
// 
//    // now get a TGraph we can use in interpolating across these runs
//    // std::vector<double> stats;
//    TGraph *gDrift = GetInterpolatedTGraph(method,fxprList,runStop[0],runStart[1],1E+9,fxprData,stats);
//    return gDrift;
// }
// //______________________________________________________________________________
// TGraph *GetDriftTGraphR2R(int method,std::vector<int> driftRun,std::vector<int> fxprList,std::vector<double> &stats){
// 
//    int rc=0;
//    std::vector<gm2field::fixedProbeFrequency_t> fxprBare,fxprGrad;
// 
//    const int NFP  = fxprList.size(); 
//    int startIndex = fxprList[0]; 
//    int stopIndex  = fxprList[NFP-1];
// 
//    // get the fxpr data       
//    rc = gm2fieldUtil::RootHelper::GetFPFrequencies(driftRun[0],fxprBare);
//    rc = gm2fieldUtil::RootHelper::GetFPFrequencies(driftRun[1],fxprGrad);
// 
//    std::vector<fixedProbeEvent_t> fxprBareAvg;
//    rc = GetAverageFXPRVectorsNew(method,0,0,0,0,fxprList,fxprBare,fxprBareAvg);
// 
//    std::vector<fixedProbeEvent_t> fxprGradAvg;
//    rc = GetAverageFXPRVectorsNew(method,0,0,0,0,fxprList,fxprGrad,fxprGradAvg);
// 
//    const int NB = fxprBareAvg.size(); 
//    const int NG = fxprGradAvg.size(); 
// 
//    // get effect of turning on SCC
//    double bareField = fxprBareAvg[NB-1].freq; 
//    double gradField = fxprGradAvg[NG-1].freq; 
//    double Fs        = gradField-bareField;
//    stats.push_back(Fs); 
//    stats.push_back(0); 
// 
//    const int NEV = 30; 
//    
//    // now get slopes of both data sets 
//    double parBare[2]; 
//    rc = FitFXPR(NB-NEV,NB,fxprBareAvg,parBare);   
// 
//    double parGrad[2]; 
//    rc = FitFXPR(0,NEV,fxprGradAvg,parGrad);  
// 
//    // store the slopes in a vector 
//    std::vector<double> dFdt; 
//    dFdt.push_back(parBare[1]); 
//    dFdt.push_back(parGrad[1]); 
// 
//    // now compute estimated field drift 
//    double dt         = fxprGradAvg[0].time - fxprBareAvg[NB-1].time;  
//    double dFdt_mean  = gm2fieldUtil::Math::GetMean<double>(dFdt);
//    double dFdt_stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(dFdt);   
//    double Fd         = dFdt_mean*dt; 
//    double Fd_err     = dFdt_stdev; 
//  
//    std::cout << "df/dt: " << std::endl;
//    for(int i=0;i<2;i++) std::cout << Form("%.3lf Hz/sec",dFdt[i]) << std::endl;
//    std::cout << Form("dt = %.3lf sec",dt) << std::endl;
// 
//    stats.push_back(Fd);     // estimated drift 
//    stats.push_back(Fd_err); // error 
// 
//    // store the data from the bare field; no corrections here 
//    std::vector<double> TIME,FREQ; 
//    for(int i=0;i<NB;i++){
//       TIME.push_back( fxprBareAvg[i].time );
//       FREQ.push_back( fxprBareAvg[i].freq );
//    } 
// 
//    // now grad field   
//    double arg=0;
//    for(int i=0;i<NG;i++){
//       arg = fxprGradAvg[i].freq - Fs + Fd; 
//       TIME.push_back( fxprGradAvg[i].time ); 
//       FREQ.push_back(arg); 
//    } 
// 
//    TGraph *gDrift = gm2fieldUtil::Graph::GetTGraph(TIME,FREQ);
//    return gDrift;
// }
