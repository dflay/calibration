#include "../include/DeltaBFuncs.h"
//______________________________________________________________________________
int CalculatePPDeltaB(bool correctDrift,int method,
                      TGraph *gDrift,
                      plungingProbeAnaEvent_t bare,plungingProbeAnaEvent_t grad,
                      double &deltaB,double &deltaB_err,
                      double &drift ,double &drift_err){
   // compute delta B for PP measurements
   // this works for whether or not we applied a drift correction to the PP data

   if(correctDrift){
      std::cout << "Correcting for field drift: " << std::endl;
   }
 
   // recall: a PP event has N sub events corresponding to taking a number of shots.  
   // first, average over these. 
   std::vector<double> fb,fg; 
   
   int NB = bare.numTraces;
   for(int i=0;i<NB;i++) fb.push_back(bare.freq[i]); 
   
   int NG = grad.numTraces;
   for(int i=0;i<NG;i++) fg.push_back(grad.freq[i]); 

   double mean_bare = gm2fieldUtil::Math::GetMean<double>(fb); 
   double mean_grad = gm2fieldUtil::Math::GetMean<double>(fg); 

   double stdev_bare = gm2fieldUtil::Math::GetStandardDeviation<double>(fb); 
   double stdev_grad = gm2fieldUtil::Math::GetStandardDeviation<double>(fg); 

   double deltaB_raw     = mean_grad - mean_bare;
   double deltaB_raw_err = TMath::Sqrt( stdev_bare*stdev_bare + stdev_grad*stdev_grad ); 

   unsigned long long timeBare = bare.time[NB-1];
   unsigned long long timeGrad = grad.time[0]; 

   std::string timeStamp_bare = gm2fieldUtil::GetStringTimeStampFromUTC(timeBare/1E+9); 
   std::string timeStamp_grad = gm2fieldUtil::GetStringTimeStampFromUTC(timeGrad/1E+9); 

   std::cout << Form("%s mean bare = %.3lf +/- %.3lf Hz",timeStamp_bare.c_str(),mean_bare,stdev_bare) << std::endl;
   std::cout << Form("%s mean grad = %.3lf +/- %.3lf Hz",timeStamp_grad.c_str(),mean_grad,stdev_grad) << std::endl;
   std::cout << Form("dB (raw) = %.3lf +/- %.3lf Hz",deltaB_raw,deltaB_raw_err) << std::endl;
 
   // now to account for drift *between* the bare and grad measurements 
   double mean_fxpr_bare  = gDrift->Eval(timeBare/1E+9); 
   double mean_fxpr_bare2 = gDrift->Eval(timeGrad/1E+9); 

   if(correctDrift){
      // assign drift values 
      drift     = mean_fxpr_bare2 - mean_fxpr_bare;
      drift_err = 0; // fabs(drift);  // assume 100% uncertainty 
      std::cout << "Data from Fixed Probes: " << std::endl;
      std::cout << Form("%s: %.3lf Hz",timeStamp_bare.c_str(),mean_fxpr_bare)  << std::endl;
      std::cout << Form("%s: %.3lf Hz",timeStamp_grad.c_str(),mean_fxpr_bare2) << std::endl;
      std::cout << Form("drift = %.3lf Hz",drift) << std::endl;
   }

   // now factor in the drift if necessary 
   if(correctDrift){
      deltaB     = deltaB_raw - drift;
      deltaB_err = deltaB_raw_err; // TMath::Sqrt(deltaB_raw_err*deltaB_raw_err + drift_err*drift_err); 
   }else{
      deltaB     = deltaB_raw;
      deltaB_err = deltaB_raw_err; 
   }
   
   if(correctDrift) std::cout << Form("dB (drift-corr)  = %.3lf +/- %.3lf Hz",deltaB,deltaB_err) << std::endl;
   std::cout << "-----------------------" << std::endl; 

   return 0; 
}
//______________________________________________________________________________
int CalculatePPDeltaB(bool correctDrift,int method,std::vector<int> trlyList,
                      std::vector<trolleyAnaEvent_t> trlyData,
                      plungingProbeAnaEvent_t bare,plungingProbeAnaEvent_t grad,
                      double &deltaB,double &deltaB_err,
                      double &drift ,double &drift_err){
   // compute delta B for PP measurements
   // this works for whether or not we applied a drift correction to the PP data

   if(correctDrift){
      std::cout << "Correcting for field drift: " << std::endl;
   }
 
   // recall: a PP event has N sub events corresponding to taking a number of shots.  
   // first, average over these. 
   std::vector<double> fb,fg; 
   
   int NB = bare.numTraces;
   for(int i=0;i<NB;i++) fb.push_back(bare.freq[i]); 
   
   int NG = grad.numTraces;
   for(int i=0;i<NG;i++) fg.push_back(grad.freq[i]); 

   double mean_bare = gm2fieldUtil::Math::GetMean<double>(fb); 
   double mean_grad = gm2fieldUtil::Math::GetMean<double>(fg); 

   double stdev_bare = gm2fieldUtil::Math::GetStandardDeviation<double>(fb); 
   double stdev_grad = gm2fieldUtil::Math::GetStandardDeviation<double>(fg); 

   double deltaB_raw     = mean_grad - mean_bare;
   double deltaB_raw_err = TMath::Sqrt( stdev_bare*stdev_bare + stdev_grad*stdev_grad ); 

   std::string timeStamp_bare = gm2fieldUtil::GetStringTimeStampFromUTC(bare.time[0]/1E+9); 
   std::string timeStamp_grad = gm2fieldUtil::GetStringTimeStampFromUTC(grad.time[0]/1E+9); 

   std::cout << Form("%s mean bare = %.3lf +/- %.3lf Hz",timeStamp_bare.c_str(),mean_bare,stdev_bare) << std::endl;
   std::cout << Form("%s mean grad = %.3lf +/- %.3lf Hz",timeStamp_grad.c_str(),mean_grad,stdev_grad) << std::endl;
   std::cout << Form("dB (raw)  = %.3lf +/- %.3lf Hz",deltaB_raw,deltaB_raw_err) << std::endl;
 
   // now to account for drift *between* the bare and grad measurements 
   double mean_trly_bare,stdev_trly_bare;
   GetAverageTRLY(bare.time[0],trlyList,trlyData,mean_trly_bare,stdev_trly_bare); 

   double mean_trly_grad,stdev_trly_grad;
   GetAverageTRLY(grad.time[0],trlyList,trlyData,mean_trly_grad,stdev_trly_grad);

   if(correctDrift){
      drift     = mean_trly_grad - mean_trly_bare;
      drift_err = fabs(drift);
      std::cout << "Data from Trolley Probes: " << std::endl;
      std::cout << Form("%s: %.3lf +/- %.3lf Hz",timeStamp_bare.c_str(),mean_trly_bare,stdev_trly_bare) << std::endl;
      std::cout << Form("%s: %.3lf +/- %.3lf Hz",timeStamp_grad.c_str(),mean_trly_grad,stdev_trly_grad) << std::endl;
      std::cout << Form("drift = %.3lf +/- %.3lf Hz",drift,drift_err) << std::endl;
   }
   std::cout << "-----------------------" << std::endl; 

   // now factor in the drift if necessary 
   if(correctDrift){
      deltaB     = deltaB_raw - drift;
      deltaB_err = deltaB_raw_err; // TMath::Sqrt(deltaB_raw_err*deltaB_raw_err + drift_err*drift_err); 
   }else{
      deltaB     = deltaB_raw;
      deltaB_err = deltaB_raw_err; 
   }

   return 0; 
}
//______________________________________________________________________________
int CalculateTRLYDeltaB_Stationary(bool correctDrift,int method,int probe,
                                   TGraph *gDrift,
                                   std::vector<trolleyAnaEvent_t> bare,std::vector<trolleyAnaEvent_t> grad,
                                   double &deltaB,double &deltaB_err,
                                   double &drift ,double &drift_err){

   // for a stationary trolley run 

   if(correctDrift){
      std::cout << "Correcting for field drift: " << std::endl;
   }

   std::vector<double> fb,fg;

   const int NB = bare.size();
   std::cout << NB << " events" << std::endl;
   for(int i=0;i<NB;i++) fb.push_back( bare[i].freq[probe] );  

   const int NG = grad.size();
   for(int i=0;i<NG;i++) fg.push_back( grad[i].freq[probe] );

   double mean_bare  = gm2fieldUtil::Math::GetMean<double>(fb);  
   double mean_grad  = gm2fieldUtil::Math::GetMean<double>(fg); 

   double stdev_bare = gm2fieldUtil::Math::GetStandardDeviation<double>(fb);  
   double stdev_grad = gm2fieldUtil::Math::GetStandardDeviation<double>(fg); 

   double deltaB_raw     = mean_grad - mean_bare; 
   double deltaB_raw_err = TMath::Sqrt( stdev_bare*stdev_bare + stdev_grad*stdev_grad ); 

   unsigned long long timeBare = bare[NB-1].time[probe];  // want the last event of the bare data! 
   unsigned long long timeGrad = grad[0].time[probe];

   std::string timeStamp_bare = gm2fieldUtil::GetStringTimeStampFromUTC(timeBare/1E+9); 
   std::string timeStamp_grad = gm2fieldUtil::GetStringTimeStampFromUTC(timeGrad/1E+9); 

   std::cout << Form("%s mean bare = %.3lf +/- %.3lf Hz",timeStamp_bare.c_str(),mean_bare,stdev_bare) << std::endl;
   std::cout << Form("%s mean grad = %.3lf +/- %.3lf Hz",timeStamp_grad.c_str(),mean_grad,stdev_grad) << std::endl;
   std::cout << Form("dB (raw)  = %.3lf +/- %.3lf Hz",deltaB_raw,deltaB_raw_err) << std::endl;
 
   // now to account for drift *between* the bare and grad measurements 
   double mean_fxpr_bare  = gDrift->Eval(timeBare/1E+9); 
   double mean_fxpr_bare2 = gDrift->Eval(timeGrad/1E+9); 

   if(correctDrift){
      drift     = mean_fxpr_bare2 - mean_fxpr_bare;
      drift_err = 0; // fabs(drift); // assume 100% uncertainty for now  
      std::cout << "Data from Fixed Probes: " << std::endl;
      std::cout << Form("%s: %.3lf Hz",timeStamp_bare.c_str(),mean_fxpr_bare)  << std::endl;
      std::cout << Form("%s: %.3lf Hz",timeStamp_grad.c_str(),mean_fxpr_bare2) << std::endl;
      std::cout << Form("drift = %.3lf Hz",drift) << std::endl;
   }

   // now factor in the drift if necessary 
   if(correctDrift){
      deltaB     = deltaB_raw - drift;
      deltaB_err = deltaB_raw_err; // TMath::Sqrt(deltaB_raw_err*deltaB_raw_err + drift_err*drift_err); 
   }else{
      deltaB     = deltaB_raw;
      deltaB_err = deltaB_raw_err; 
   }
   
   if(correctDrift) std::cout << Form("dB (drift-cor) = %.3lf +/- %.3lf Hz",deltaB,deltaB_err) << std::endl;
   std::cout << "-----------------------" << std::endl; 

   return 0;
}
//______________________________________________________________________________
int CalculateTRLYDeltaB_Stationary(bool correctDrift,int method,int probe,std::vector<int> trlyList,
                                   std::vector<trolleyAnaEvent_t> trlyData,
                                   std::vector<trolleyAnaEvent_t> bare,std::vector<trolleyAnaEvent_t> grad,
                                   double &deltaB,double &deltaB_err,
                                   double &drift ,double &drift_err){
   // for a stationary trolley run 

   if(correctDrift){
      std::cout << "Correcting for field drift: " << std::endl;
   }

   std::vector<double> fb,fg;

   const int NB = bare.size();
   for(int i=0;i<NB;i++) fb.push_back( bare[i].freq[probe] );  

   const int NG = grad.size();
   for(int i=0;i<NG;i++) fg.push_back( grad[i].freq[probe] );

   double mean_bare  = gm2fieldUtil::Math::GetMean<double>(fb);  
   double mean_grad  = gm2fieldUtil::Math::GetMean<double>(fg); 

   double stdev_bare = gm2fieldUtil::Math::GetStandardDeviation<double>(fb);  
   double stdev_grad = gm2fieldUtil::Math::GetStandardDeviation<double>(fg); 

   double deltaB_raw     = mean_grad - mean_bare; 
   double deltaB_raw_err = TMath::Sqrt( stdev_bare*stdev_bare + stdev_grad*stdev_grad ); 

   unsigned long long timeBare = bare[0].time[probe];
   unsigned long long timeGrad = grad[0].time[probe];

   std::string timeStamp_bare = gm2fieldUtil::GetStringTimeStampFromUTC(timeBare/1E+9); 
   std::string timeStamp_grad = gm2fieldUtil::GetStringTimeStampFromUTC(timeGrad/1E+9); 

   std::cout << Form("%s mean bare = %.3lf +/- %.3lf Hz",timeStamp_bare.c_str(),mean_bare,stdev_bare) << std::endl;
   std::cout << Form("%s mean grad = %.3lf +/- %.3lf Hz",timeStamp_grad.c_str(),mean_grad,stdev_grad) << std::endl;
   std::cout << Form("dB (raw)  = %.3lf +/- %.3lf Hz",deltaB_raw,deltaB_raw_err) << std::endl;
 
   // now to account for drift *between* the bare and grad measurements 
   double mean_trly_bare,stdev_trly_bare;
   GetAverageTRLY(timeBare,trlyList,trlyData,mean_trly_bare,stdev_trly_bare); 

   double mean_trly_grad,stdev_trly_grad;
   GetAverageTRLY(timeGrad,trlyList,trlyData,mean_trly_grad,stdev_trly_grad);

   if(correctDrift){
      drift     = mean_trly_grad - mean_trly_bare;
      drift_err = fabs(drift);  // assume 100% uncertainty for now 
      std::cout << "Data from Trolley Probes: " << std::endl;
      std::cout << Form("%s: %.3lf Hz",timeStamp_bare.c_str(),mean_trly_bare) << std::endl;
      std::cout << Form("%s: %.3lf Hz",timeStamp_grad.c_str(),mean_trly_grad) << std::endl;
      std::cout << Form("drift = %.3lf +/- %.3lf Hz",drift,drift_err) << std::endl;
   }
   std::cout << "-----------------------" << std::endl; 

   // now factor in the drift if necessary 
   if(correctDrift){
      deltaB     = deltaB_raw - drift;
      deltaB_err = deltaB_raw_err; // TMath::Sqrt(deltaB_raw_err*deltaB_raw_err + drift_err*drift_err); 
   }else{
      deltaB     = deltaB_raw;
      deltaB_err = deltaB_raw_err; 
   }

   return 0;
}
//______________________________________________________________________________
TGraph *CalculateTRLYDeltaB_Moving(int method,int probe,std::vector<int> fxprList,std::vector<gm2field::fixedProbeFrequency_t> fxprData,
                                  std::vector<trolleyAnaEvent_t> bare,std::vector<trolleyAnaEvent_t> grad){
   // compute delta B for trolley measurements
   // since we don't know if two runs have the same position, should probably fit to a spline and subtract... 
   std::vector<double> tb,fb,tg,fg;

   // gather bare field data
   int NB = bare.size();
   for(int i=0;i<NB;i++){
      tb.push_back( bare[i].time[probe]/1E+9 ); 
      fb.push_back( bare[i].freq[probe] ); 
   }
   
   TGraph *gb = gm2fieldUtil::Graph::GetTGraph(tb,fb);
 
   // gather gradient field data
   int NG = grad.size();
   for(int i=0;i<NG;i++){
      tg.push_back( grad[i].time[probe]/1E+9 ); 
      fg.push_back( grad[i].freq[probe] ); 
   }

   TGraph *gg = gm2fieldUtil::Graph::GetTGraph(tg,fg);

   // determine bounds for new plot 
   // which data set came first?
   unsigned long long t0b = bare[0].time[probe]; 
   unsigned long long t0g = grad[0].time[probe];

   bool bare_came_first=false;
   double timeStart=0,timeStop=0,timeStep=1.; 
   if(t0b<=t0g){
      bare_came_first = true;
      timeStart = t0b/1E+9;
      timeStop = grad[NG-1].time[probe]/1E+9;
   }else{
      timeStart = t0g/1E+9;
      timeStop = bare[NB-1].time[probe]/1E+9;
   } 

   if(timeStop<timeStart){
      std::cout << "ERROR! timeStart = " << timeStart << " timeStop = " << timeStop << std::endl;
      exit(1);
   }

   // find drift between the datasets 
   double mean_fxpr_bare=0,stdev_fxpr_bare=0;
   GetAverageFXPR(method,bare[0].time[probe],fxprList,fxprData,mean_fxpr_bare,stdev_fxpr_bare); 
   
   double mean_fxpr_grad=0,stdev_fxpr_grad=0;
   GetAverageFXPR(method,grad[0].time[probe],fxprList,fxprData,mean_fxpr_grad,stdev_fxpr_grad);

   // basically assuming the bare came last here intentionally 
   double drift = mean_fxpr_bare - mean_fxpr_grad;
   if(bare_came_first) drift *= -1; 

   // find DeltaB between the graphs, correcting for field drift  
   double arg_t=0,arg_f=0;
   std::vector<double> T,F; 
   int NPTS = (timeStop-timeStart)/timeStep;
   for(int i=0;i<NPTS;i++){
      arg_t = timeStart + i*timeStep;
      arg_f = gg->Eval(arg_t) - gb->Eval(arg_t) - drift;  
      T.push_back(arg_t);
      F.push_back(arg_f);
   }
  
   TGraph *gDeltaB = gm2fieldUtil::Graph::GetTGraph(T,F);
   return gDeltaB; 
 
   return 0; 
}
