// Process the PP data  
// produce an output file based on the input run list containing 
// the time stamp and raw frequencies  

#include <cstdlib> 
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>  
#include <cmath>
#include <string>  

#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSpline.h"

#include "RootTreeStructs.h"
#include "gm2fieldMath.h"
#include "gm2fieldGraph.h"
#include "gm2fieldRootHelper.h"
#include "gm2fieldUnits.h"
#include "TemperatureSensor.h"
#include "MovingAverage.h"
#include "Blinder.h"

#include "./include/date.h"
#include "./include/Constants.h"
#include "./include/fixedProbeEvent.h"
#include "./include/trolleyAnaEvent.h"
#include "./include/sccEvent.h"
#include "./include/nmr_meas.h" 

#include "./src/CustomUtilities.C"
#include "./src/CustomImport.C"
#include "./src/CustomExport.C"
#include "./src/CustomMath.C"
#include "./src/CustomGraph.C"
#include "./src/CustomAlgorithms.C"
#include "./src/OscFuncs.C"
#include "./src/BlindFuncs.C"
#include "./src/InputManager.C"
#include "./src/FitFuncs.C"
#include "./src/TRLYFuncs.C"
// #include "./src/CalibFuncs.C"
#include "./src/FreeProton.C"

double gMarkerSize = 0.8; 

TGraphErrors *GetPPFreeGraph(TString xAxis,TString yAxis,std::vector<calibSwap_t> data); 

int GetSwapStatsForPP(int runPeriod,std::string ppID,std::vector<plungingProbeAnaEvent_t> data,
                      std::vector<calibSwap_t> &raw,std::vector<calibSwap_t> &free,double &min,double &max);

int PrintToFile(std::string outpath,std::vector<calibSwap_t> data); 

int Process_pp_prod(std::string configFile){

   std::cout << "------------------------------------" << std::endl;
   std::cout << "GET PP SWAPPING DATA" << std::endl;

   int rc=0;
   int prMethod = gm2fieldUtil::Constants::kPhaseDerivative;
   int ppMethod = plungingProbeAnalysis::kLeastSquaresPhase;

   InputManager *inputMgr = new InputManager();
   inputMgr->Load(configFile);
   inputMgr->Print();

   std::string prodVersion   = inputMgr->GetProductionTag();
   std::string nmrAnaVersion = inputMgr->GetNMRANATag();
   std::string runDate       = inputMgr->GetRunDate();
   std::string blindLabel    = inputMgr->GetBlindLabel();
   std::string cutFile       = inputMgr->GetCutFile();
   std::string ppID          = inputMgr->GetPPID();  
 
   bool isBlind              = inputMgr->IsBlind();
   bool useOscCor            = inputMgr->GetOscCorStatus();
   int probeNumber           = inputMgr->GetTrolleyProbe();
   int runPeriod             = inputMgr->GetRunPeriod(); 
   int nev                   = inputMgr->GetNumEventsToAvg();
   // systematics 
   bool isSyst             = inputMgr->GetSystStatus();
   int systDirNum          = inputMgr->GetSystDirNum();

   double tempCorValue = 0;
   bool useTempCor_pp  = inputMgr->GetTempCorStatus_pp(); 
   if(useTempCor_pp) tempCorValue = inputMgr->GetTempCor_pp();  

   char cutPath[200]; 
   sprintf(cutPath,"./input/json/run-%d/%s",runPeriod,cutFile.c_str());
   std::string cutpath = cutPath;  

   // date_t theDate;
   // GetDate(theDate);
   std::string theDate = inputMgr->GetAnaDate();

   std::string plotDir = GetPath("plots" ,isBlind,blindLabel,theDate,isSyst,systDirNum);
   std::string outDir  = GetPath("output",isBlind,blindLabel,theDate,isSyst,systDirNum);

   char outPath[500]; 
   sprintf(outPath,"%s/pp-swap-data_pr-%02d.csv",outDir.c_str(),probeNumber);
   std::string outpath_raw = outPath;  
   sprintf(outPath,"%s/pp-swap-data_free-prot_pr-%02d.csv",outDir.c_str(),probeNumber);
   std::string outpath_free = outPath; 
   sprintf(outPath,"%s/pp-swap_fpc_pr-%02d.csv",outDir.c_str(),probeNumber);
   std::string outpath_fpc = outPath; 

   int blindUnits  = inputMgr->GetBlindUnits(); 
   double blindMag = inputMgr->GetBlindScale(); 
   gm2fieldUtil::Blinder *myBlind = new gm2fieldUtil::Blinder(blindLabel,blindMag,blindUnits);
   double blindValue = myBlind->GetBlinding(1); // in Hz

   std::vector<int> run;
   std::vector<std::string> label;
   inputMgr->GetRunList(run);
   inputMgr->GetRunLabels(label);

   const int NRUNS = run.size();

   double yMin=0,yMax=0;

   TString ppPlotPath;  

   std::vector<int> fxprList; 
   inputMgr->GetFXPRList(fxprList);

   bool subtractDrift = inputMgr->GetFXPRDriftStatus();  
   std::vector<averageFixedProbeEvent_t> fxprData;
   int period = inputMgr->GetNumEventsTimeWindow(); // 10;
   for(int i=0;i<NRUNS;i++){
      rc = GetFixedProbeData_avg(run[i],prMethod,fxprList,fxprData,prodVersion,subtractDrift,period,0);
      if(rc!=0){
	 std::cout << "No data!" << std::endl;
	 return 1;
      }
   }

   // PP data 
   bool useNMRANA = true;
   std::vector<plungingProbeAnaEvent_t> ppInput,ppData;

   int NPP=0; 
   for(int i=0;i<NRUNS;i++){
      std::cout << "Getting PP data for run " << run[i] << "..." << std::endl;
      rc = GetPlungingProbeData(run[i],prMethod,ppMethod,ppInput,prodVersion,nmrAnaVersion,cutpath,useNMRANA,tempCorValue);
      if(rc!=0){
	 std::cout << "No data!" << std::endl;
	 return 1;
      }
   }
  
   if(isBlind) ApplyBlindingPP(blindValue,ppInput);

   std::vector<trolleyAnaEvent_t> trlyData;
   for(int i=0;i<NRUNS;i++){
      rc = GetTrolleyData(run[i],prMethod,trlyData,prodVersion);
      if(rc!=0){
	 std::cout << "No data!" << std::endl;
	 return 1;
      }
   }

   // get reference time 
   std::vector<double> time;
   rc = LoadTimes(probeNumber,runPeriod,prodVersion,"swap","tr",time);
   double t0 = GetSwappingT0(probeNumber-1,nev,time,fxprData,ppInput,trlyData);
   
   // oscillation correction 
   if(useOscCor){
      rc = CorrectOscillation_pp(fxprData,ppInput,ppData,t0);
   }else{
      CopyPlungingProbe(ppInput,ppData);
   }

   // load in perturbation data 
   // char inpath[200]; 
   // perturbation_t ppPert;
   // sprintf(inpath,"./input/perturbation/pp-pert_run-%d.json",runPeriod);
   // LoadPerturbationData_json(inpath,ppPert);

   int M=0;

   // get stats for the PP
   std::vector<double> fpc,fpcErr;
   std::vector<calibSwap_t> raw,free; 
   // rc = GetSwapStatsForPP(ppPert,ppData,raw,free,yMin,yMax); 
   rc = GetSwapStatsForPP(runPeriod,ppID,ppData,raw,free,yMin,yMax,fpc,fpcErr); 

   // print to file 
   rc = PrintToFile(outpath_raw ,raw ); 
   rc = PrintToFile(outpath_free,free);
   rc = PrintToFile_2dbl(outpath_fpc.c_str(),fpc,fpcErr);

   // make plot of all frequencies for all shots 
   TGraph *gPP     = GetPPTGraph1("TimeStamp","freq",ppData); 
   gm2fieldUtil::Graph::SetGraphParameters(gPP,20,kBlack);
   gPP->SetMarkerSize(gMarkerSize); 

   std::vector<double> Time,FreqFree,FreqFreeErr; 

   TGraph *gPP_free = GetPPFreeGraph("Time","Freq",free);
   gm2fieldUtil::Graph::SetGraphParameters(gPP_free,20,kBlack);
   gPP_free->SetMarkerSize(gMarkerSize); 
   
   TCanvas *c1 = new TCanvas("c1","PP Data",1200,600);
   c1->Divide(1,2);

   c1->cd(1); 
   gPP->SetTitle( Form("Run %d",run[0]) );
   gPP->Draw("alp");  
   gPP->GetYaxis()->SetRangeUser(yMin,yMax);
   gm2fieldUtil::Graph::UseTimeDisplay(gPP); 
   gPP->Draw("alp");  
   c1->Update();

   c1->cd(2); 
   gPP_free->SetTitle( Form("Run %d",run[0]) );
   gPP_free->Draw("alp");  
   gm2fieldUtil::Graph::UseTimeDisplay(gPP_free); 
   gPP_free->Draw("alp");  
   c1->Update();

   // save the plot  
   c1->cd(); 
   ppPlotPath = Form("%s/pp-swap-data_pr-%02d.png",plotDir.c_str(),probeNumber);
   c1->Print(ppPlotPath);

   return 0;
}
//______________________________________________________________________________
int GetSwapStatsForPP(int runPeriod,std::string probeID,std::vector<plungingProbeAnaEvent_t> data,
                      std::vector<calibSwap_t> &raw,std::vector<calibSwap_t> &free,double &min,double &max,
                      std::vector<double> &freeProtCor,std::vector<double> &freeProtCorErr){

   // gather stats for the PP; recall a single PP-DAQ run has N traces in it. 
   // here, a PP-DAQ run is a single swap!

   FreeProton *fp = new FreeProton(probeID,runPeriod);

   min = 50E+3; 
   max = 0; 
 
   int M=0,nmrDAQ_run=0;
   const int N = data.size();

   calibSwap_t rawEvent,freeEvent;
 
   double arg=0;
   double delta_t=0,delta_t_err=0;
   double mean_freeProt=0,mean_freeProtErr=0,stdev_freeProt=0;
   double mean_freq_free=0,stdev_freq_free=0;
   double mean_freq=0,stdev_freq=0;
   double mean_temp=0,stdev_temp=0;
   double mean_x=0,stdev_x=0;
   double mean_y=0,stdev_y=0;
   double mean_z=0,stdev_z=0;
   double freqFree=0,freqFreeErr_stat=0,freqFreeErr_syst=0;
   std::vector<double> F,FF,T,x,y,z,fpc,fpc_err;

   int rc=0; 
   char msg[200];
   std::string MSG; 

   for(int i=0;i<N;i++){
      M = data[i].numTraces;
      nmrDAQ_run = data[i].run;
      // shouldn't need this anymore, since it's in the cut file  
      // if(nmrDAQ_run==888 || nmrDAQ_run==834 || nmrDAQ_run==789){
      //    std::cout << "WARNING: Known issue in run " << nmrDAQ_run << "!  Skipping..." << std::endl;
      //    continue;
      // }
      if(M==0){
	 sprintf(msg,"[Process_pp_prod]: WARNING: No traces for PP run %d, continuing on...",nmrDAQ_run);
	 MSG = msg; 
         rc = Logger::PrintMessage(Logger::kINFO,"default",MSG,'a');  
	 continue; 
      } 
      for(int j=0;j<M;j++){
         // we want to add back in the LO since we're comparing against the TRLY eventually 
	 arg = data[i].freq[j] + data[i].freq_LO[j]; 
	 // get the free-proton frequency here
	 // GetOmegaP_free(ppPert,arg,0,data[i].temp[j],0,freqFree,freqFreeErr_stat);
         freqFree    = fp->GetOmegaP_free(arg,data[i].temp[j]);
         delta_t     = fp->GetDelta_t(data[i].temp[j])/1E-9;     // in ppb 
         delta_t_err = fp->GetDelta_t_err(data[i].temp[j])/1E-9; // in ppb  
	 F.push_back(arg);
	 FF.push_back(freqFree);
         T.push_back(data[i].temp[j]);
	 x.push_back(data[i].r[j]);  
	 y.push_back(data[i].y[j]);  
	 z.push_back(data[i].phi[j]); 
         fpc.push_back(delta_t);
         fpc_err.push_back(delta_t_err); 
         if( min>data[i].freq[j] ) min = data[i].freq[j];  
         if( max<data[i].freq[j] ) max = data[i].freq[j];  
      }
      // get average on the run (swap)
      mean_freeProt   = gm2fieldUtil::Math::GetMean<double>(fpc);
      mean_freeProtErr= gm2fieldUtil::Math::GetMean<double>(fpc_err);
      stdev_freeProt  = gm2fieldUtil::Math::GetStandardDeviation<double>(fpc); 
      mean_freq       = gm2fieldUtil::Math::GetMean<double>(F); 
      stdev_freq      = gm2fieldUtil::Math::GetStandardDeviation<double>(F);
      mean_freq_free  = gm2fieldUtil::Math::GetMean<double>(FF); 
      stdev_freq_free = gm2fieldUtil::Math::GetStandardDeviation<double>(FF);
      mean_temp       = gm2fieldUtil::Math::GetMean<double>(T); 
      stdev_temp      = gm2fieldUtil::Math::GetStandardDeviation<double>(T);
      mean_x          = gm2fieldUtil::Math::GetMean<double>(x); 
      stdev_x         = gm2fieldUtil::Math::GetStandardDeviation<double>(x);
      mean_y          = gm2fieldUtil::Math::GetMean<double>(y); 
      stdev_y         = gm2fieldUtil::Math::GetStandardDeviation<double>(y);
      mean_z          = gm2fieldUtil::Math::GetMean<double>(z); 
      stdev_z         = gm2fieldUtil::Math::GetStandardDeviation<double>(z);
      // apply free proton corrections 
      // GetOmegaP_free(ppPert,mean_freq,stdev_freq,mean_temp,stdev_temp,freqFree,freqFreeErr);
      // fill output structs 
      rawEvent.time    = data[i].time[0]/1E+9; 
      rawEvent.freq    = mean_freq;  
      rawEvent.freqErr = stdev_freq;  
      rawEvent.temp    = mean_temp;  
      rawEvent.tempErr = stdev_temp;  
      rawEvent.r       = mean_x;  
      rawEvent.rErr    = stdev_x;  
      rawEvent.y       = mean_y;  
      rawEvent.yErr    = stdev_y;  
      rawEvent.phi     = mean_z;  
      rawEvent.phiErr  = stdev_z; 
      // free proton  
      freeEvent.time    = data[i].time[0]/1E+9; 
      freeEvent.freq    = mean_freq_free;  
      freeEvent.freqErr = stdev_freq_free;  // WARNING: We'll fold in free-proton errors seprately.  we want shot errors here
      freeEvent.temp    = mean_temp;  
      freeEvent.tempErr = stdev_temp;  
      freeEvent.r       = mean_x;  
      freeEvent.rErr    = stdev_x;  
      freeEvent.y       = mean_y;  
      freeEvent.yErr    = stdev_y;  
      freeEvent.phi     = mean_z;  
      freeEvent.phiErr  = stdev_z; 
      // fill vectors 
      raw.push_back(rawEvent);  
      free.push_back(freeEvent); 
      freeProtCor.push_back(mean_freeProt);
      freeProtCorErr.push_back(mean_freeProtErr);
      // clean up for next run 
      fpc.clear();
      fpc_err.clear();
      F.clear();
      FF.clear();
      T.clear();  
      x.clear();  
      y.clear();  
      z.clear();  
   }

   min -= 2; 
   max += 2;

   delete fp; 

   return 0;
}
//______________________________________________________________________________
TGraphErrors *GetPPFreeGraph(TString xAxis,TString yAxis,std::vector<calibSwap_t> data){
   const int N = data.size(); 
   std::vector<double> x,y,ey; 
   for(int i=0;i<N;i++){
      if(xAxis=="Time") x.push_back( data[i].time ); 
      if(yAxis=="Freq"){
	 y.push_back( data[i].freq );
	 ey.push_back( data[i].freqErr );
      }else if(yAxis=="Temp"){
	 y.push_back( data[i].temp );
	 ey.push_back( data[i].tempErr );
      } 
   }
   TGraphErrors *g = gm2fieldUtil::Graph::GetTGraphErrors(x,y,ey); 
   return g;
}
//______________________________________________________________________________
int PrintToFile(std::string outpath,std::vector<calibSwap_t> data){

   const int N = data.size(); 
   char myStr[1000]; 

   std::ofstream outfile;
   outfile.open( outpath.c_str() );
   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      for(int i=0;i<N;i++){
	 sprintf(myStr,"%.0lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf",
                 data[i].time,data[i].freq,data[i].freqErr,data[i].temp,data[i].tempErr,
                 data[i].r,data[i].rErr,data[i].y,data[i].yErr,data[i].phi,data[i].phiErr); 
	 outfile << myStr << std::endl;
      } 
      std::cout << "The data has been written to file: " << outpath << std::endl;
      outfile.close();
   }  
   return 0;
}
