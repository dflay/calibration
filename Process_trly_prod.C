// Process the trolley data; identify each time the 
// trolley stopped at the calibration region and extract 
// the mean and stdev frequency for the given trolley probe 
// with the appropriate time stamp   

#include <cstdlib> 
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>  
#include <cmath> 

#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSpline.h"
#include "TLine.h"

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
#include "./src/SystFuncs.C"

double gMarkerSize = 0.8; 

TGraph *GetTRLYVelocityTGraph(int probe,TString xAxis,TString yAxis,
                              std::vector<gm2field::galilTrolley_t> galil,
                              std::vector<trolleyAnaEvent_t> trly); 

int GetTRLYSwapPlots(bool useOscCor,int probe,int nev,std::vector<double> time,
                     std::vector<averageFixedProbeEvent_t> fxpr,std::vector<trolleyAnaEvent_t> trlyData,
                     TGraph **gFreq,TGraph **gTemp); 

int GetTRLYSwapData(bool useOscCor,int probe,int nev,double time,std::vector<averageFixedProbeEvent_t> fxpr,
                    std::vector<trolleyAnaEvent_t> Data,
                    std::vector<double> &TIME,std::vector<double> &FREQ,std::vector<double> &TEMP,
		    std::vector<double> &X,std::vector<double> &Y,std::vector<double> &Z); 

int PrintToFile(std::string outpath,std::vector<calibSwap_t> Event); 

int Process_trly_prod(std::string configFile){

   std::cout << "------------------------------------" << std::endl;
   std::cout << "GET TRLY SWAPPING DATA" << std::endl;

   int rc=0;
   int method = gm2fieldUtil::Constants::kPhaseDerivative;

   InputManager *inputMgr = new InputManager();
   inputMgr->Load(configFile);
   inputMgr->Print();

   std::string prodVersion   = inputMgr->GetProductionTag();  
   std::string nmrAnaVersion = inputMgr->GetNMRANATag();  
   std::string runDate       = inputMgr->GetRunDate();
   std::string blindLabel    = inputMgr->GetBlindLabel();
   std::string cutFile       = inputMgr->GetCutFile();
 
   bool isBlind              = inputMgr->IsBlind();
   bool loadSwapTimes        = inputMgr->GetSwapTimeStatus();
   bool useTempCor           = inputMgr->GetTempCorStatus(); 
   bool useOscCor            = inputMgr->GetOscCorStatus(); 
   int axis                  = inputMgr->GetAxis();  
   int probeNumber           = inputMgr->GetTrolleyProbe(); 
   int runPeriod             = inputMgr->GetRunPeriod();
   int nev                   = inputMgr->GetNumEventsToAvg();
   // systematics 
   bool isSyst               = inputMgr->GetSystStatus();
   bool varySwap_time        = inputMgr->GetVaryTimeStatus("tr","swap");
   int systDirNum            = inputMgr->GetSystDirNum();
   double swap_delta         = inputMgr->GetDeltaTime("tr","swap");   

   // date_t theDate;
   // GetDate(theDate);
   std::string theDate = inputMgr->GetAnaDate();

   std::string plotDir = GetPath("plots" ,isBlind,blindLabel,theDate,isSyst,systDirNum); 
   std::string outDir  = GetPath("output",isBlind,blindLabel,theDate,isSyst,systDirNum);

   char cutPath[200];
   sprintf(cutPath,"./input/json/run-%d/%s",runPeriod,cutFile.c_str());
   std::string cutpath = cutPath; 

   int blindUnits  = inputMgr->GetBlindUnits(); 
   double blindMag = inputMgr->GetBlindScale(); 
   gm2fieldUtil::Blinder *myBlind = new gm2fieldUtil::Blinder(blindLabel,blindMag,blindUnits);
   double blindValue = myBlind->GetBlinding(2); // in Hz 

   std::vector<int> run;
   std::vector<std::string> label;
   inputMgr->GetRunList(run);
   inputMgr->GetRunLabels(label);

   std::vector<gm2field::galilTrolley_t> trlyGalil; 
   std::vector<trolleyAnaEvent_t> trlyData;   
 
   const int NRUNS = run.size();
   for(int i=0;i<NRUNS;i++){
      rc = GetTrolleyData(run[i],method,trlyData,prodVersion);
      if(rc!=0){
	 std::cout << "No data!" << std::endl;
	 return 1;
      }
      rc = GetTrolleyGalil<gm2field::galilTrolley_t>(run[i],trlyGalil,prodVersion);
      if(rc!=0){
	 std::cout << "No data!" << std::endl;
	 return 1;
      }
   }
   
   if(isBlind) ApplyBlindingTRLY(blindValue,trlyData);
 
   std::vector<int> fxprList;
   inputMgr->GetFXPRList(fxprList);

   bool subtractDrift = inputMgr->GetFXPRDriftStatus();  
   std::vector<averageFixedProbeEvent_t> fxprData;  
   int period = inputMgr->GetNumEventsTimeWindow();
   for(int i=0;i<NRUNS;i++){
      rc = GetFixedProbeData_avg(run[i],method,fxprList,fxprData,prodVersion,subtractDrift,period,0);
      if(rc!=0){
	 std::cout << "No data!" << std::endl;
	 return 1;
      }
   }

   // TGraph *gSig  = GetTRLYVelocityTGraph(probeNumber-1,"GpsTimeStamp","vSig" ,trlyGalil,trlyData);
   // TGraph *gFish = GetTRLYVelocityTGraph(probeNumber-1,"GpsTimeStamp","vFish",trlyGalil,trlyData);
   // 
   // gm2fieldUtil::Graph::SetGraphParameters(gSig ,20,kBlue);    
   // gm2fieldUtil::Graph::SetGraphParameters(gFish,21,kRed );    

   // TMultiGraph *mgVelo = new TMultiGraph(); 
   // mgVelo->Add(gSig ,"lp"); 
   // mgVelo->Add(gFish,"lp"); 

   std::vector<double> time;

   // find the times for each TRLY swap
   double angle = 189.160; // nominal trolley location 

   std::string dev = "tr";
   // if(axis==-1) dev = "NONE";  
   // if(axis==0) dev = "trx"; 
   // if(axis==1) dev = "try"; 
   // if(axis==2) dev = "trz"; 

   if(loadSwapTimes){
      LoadTimes(probeNumber,runPeriod,prodVersion,"swap",dev,time);
   }else{
      FindTRLYStopTimes(probeNumber-1,angle,trlyData,trlyGalil,time);
   }

   if(isSyst && varySwap_time){
      std::cout << "[Process_trly_prod]: SYSTEMATIC VARIATION! Swap times will be randomized by up to " << swap_delta << " sec" <<std::endl;
      rc = systFunc::RandomizeTimeValues(swap_delta,time);
   }

   // get reference time
   int ppMethod = plungingProbeAnalysis::kLeastSquaresPhase; 
   std::vector<plungingProbeAnaEvent_t> ppData;
   for(int i=0;i<NRUNS;i++) rc = GetPlungingProbeData(run[i],method,ppMethod,ppData,prodVersion,nmrAnaVersion,cutpath);
   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }
 
   double t0 = GetSwappingT0(probeNumber-1,nev,time,fxprData,ppData,trlyData);

   const int NL = time.size();
   std::cout << "Found " << NL << " trolley swaps" << std::endl;

   // now get average field for 30 seconds BEFORE each time stamp
   double fLO = 61.74E+6; 
   double dsigdT = inputMgr->GetTempCor_tr();
   std::vector<calibSwap_t> trlySwap; 
   rc = GetTRLYStatsAtTime(useTempCor,useOscCor,probeNumber-1,nev,fLO,time,fxprData,trlyData,trlySwap,t0,dsigdT);

   char outPath[500];
   sprintf(outPath,"%s/trly-swap-data_pr-%02d.csv",outDir.c_str(),probeNumber);
   std::string outpath = outPath;

   rc = PrintToFile(outpath,trlySwap); 

   // make plots 
   TGraph *gTR = GetTRLYTGraph(probeNumber-1,"GpsTimeStamp","freq",trlyData);
   gm2fieldUtil::Graph::SetGraphParameters(gTR,20,kBlack); 

   TGraph **gFreq = new TGraph*[NL]; 
   TGraph **gTemp = new TGraph*[NL];
   rc = GetTRLYSwapPlots(useOscCor,probeNumber-1,nev,time,fxprData,trlyData,gFreq,gTemp); 

   TMultiGraph *mgf = new TMultiGraph();
   TMultiGraph *mgt = new TMultiGraph();
   for(int i=0;i<NL;i++){
      mgf->Add(gFreq[i],"lp");
      mgt->Add(gTemp[i],"lp");
   }

   double mean=0,stdev=0;
   double mean_t=0,stdev_t=0;
   rc = GetStats("freq",probeNumber-1,trlyData,mean,stdev); 
   rc = GetStats("temp",probeNumber-1,trlyData,mean_t,stdev_t); 

   double yMin  = mean-1E+3; 
   double yMax  = mean+1E+3;

   rc = GetStats("freq",trlySwap,mean,stdev); 

   double yMin2 = (mean-fLO)-10.*stdev;
   double yMax2 = (mean-fLO)+10.*stdev;

   double yMin3 = mean_t-0.2;
   double yMax3 = mean_t+0.2;

   TLine **tSwap = new TLine*[NL]; 
   for(int i=0;i<NL;i++){
      tSwap[i] = new TLine(time[i],yMin,time[i],yMax);
      tSwap[i]->SetLineColor(kRed);
      tSwap[i]->SetLineWidth(2);  
      tSwap[i]->SetLineStyle(2); 
   } 

   std::cout << "Making plots..." << std::endl;

   TCanvas *c1 = new TCanvas("c1","TRLY Data",1200,600); 
   c1->Divide(1,3);

   c1->cd(1);
   gTR->Draw("ap"); 
   gm2fieldUtil::Graph::SetGraphLabels(gTR,Form("TRLY-%02d Swap Data",probeNumber),"","Frequency (Hz)"); 
   gm2fieldUtil::Graph::UseTimeDisplay(gTR);
   gm2fieldUtil::Graph::SetGraphLabelSizes(gTR,0.05,0.06); 
   gTR->GetYaxis()->SetRangeUser(yMin,yMax); 
   gTR->Draw("ap");
   for(int i=0;i<NL;i++) tSwap[i]->Draw("same"); 
   c1->Update();

   c1->cd(2);
   mgf->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mgf,Form("TRLY-%02d Swap Events",probeNumber),"","Frequency (Hz)");
   gm2fieldUtil::Graph::SetGraphLabelSizes(mgf,0.05,0.06);
   mgf->GetYaxis()->SetRangeUser(yMin2,yMax2);
   mgf->Draw("a");
   c1->Update();

   c1->cd(3);
   mgt->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mgt,Form("TRLY-%02d Swap Events",probeNumber),"","Temperature (#circC)");
   gm2fieldUtil::Graph::SetGraphLabelSizes(mgt,0.05,0.06);
   mgt->GetYaxis()->SetRangeUser(yMin3,yMax3);
   mgt->Draw("a");
   c1->Update();
   
   TString plotPath = Form("%s/trly-swap_pr-%02d.png",plotDir.c_str(),probeNumber);
   c1->cd(); 
   c1->Print(plotPath); 
   delete c1; 

   // TCanvas *c2 = new TCanvas("c2","TRLY Velocities",1200,600); 
   // c2->cd();

   // mgVelo->Draw("a"); 
   // gm2fieldUtil::Graph::SetGraphLabels(mgVelo,"TRLY Velocities (blue = signal, red = fishing)","","Velocity (counts/sec)"); 
   // gm2fieldUtil::Graph::UseTimeDisplay(mgVelo); 
   // gm2fieldUtil::Graph::SetGraphLabelSizes(mgVelo,0.05,0.06); 
   // mgVelo->Draw("a");
   // c2->Update();

   // plotPath = Form("%s/trly-velo_run-%05d.png",plotDir,run[0]);
   // c2->cd(); 
   // c2->Print(plotPath); 

   std::cout << "--> Done." << std::endl; 

   return 0;
}
//______________________________________________________________________________
TGraph *GetTRLYVelocityTGraph(int probe,TString xAxis,TString yAxis,
                              std::vector<gm2field::galilTrolley_t> galil,
                              std::vector<trolleyAnaEvent_t> trly){
   
   // galil is in ms; needs a GPS offset 
   double sf    = 2597./2658.;   // (43 min, 17 sec on trly/44 min, 18 sec on galil) 
   double t0_tr = trly[0].time[probe]/1E+9; 
   double t0_ga = sf*galil[0].TimeStamp/1E+3;
   double dt    = t0_tr-t0_ga;

   double arg=0;
   std::vector<double> x,y; 
   const int N = galil.size();
   for(int i=0;i<N;i++){
	 arg = sf*galil[i].TimeStamp/1E+3;
      if(xAxis=="GpsTimeStamp"){
	 arg += dt;
	 x.push_back(arg); 
      }else if(xAxis=="GalilTimeStamp"){
	 x.push_back(arg);
      } 
      if(yAxis=="vFish") y.push_back(galil[i].Velocities[0]); 
      if(yAxis=="vSig")  y.push_back(galil[i].Velocities[1]); 
   }

   const int NX = x.size();
   if(NX==0){
      std::cout << "No x data!" << std::endl;
      return NULL;
   }

   TGraph *g = gm2fieldUtil::Graph::GetTGraph(x,y);
   return g;
}
//______________________________________________________________________________
int GetTRLYSwapPlots(bool useOscCor,int probe,int nev,std::vector<double> time,
                     std::vector<averageFixedProbeEvent_t> fxpr,std::vector<trolleyAnaEvent_t> trlyData,
                     TGraph **gFreq,TGraph **gTemp){
   // for diagnostic plots
   int M=0,rc=0;
   const int NL = time.size();
   std::vector<double> EV,TIME,FREQ,TEMP,X,Y,Z;
   for(int i=0;i<NL;i++){
      rc = GetTRLYSwapData(useOscCor,probe,nev,time[i],fxpr,trlyData,TIME,FREQ,TEMP,X,Y,Z);
      M = FREQ.size();
      if(M==0) std::cout << "[GetTRLYSwapPlots]: NO EVENTS!" << std::endl;
      for(int j=0;j<M;j++) EV.push_back(j+1);
      gFreq[i] = gm2fieldUtil::Graph::GetTGraph(EV,FREQ);
      gTemp[i] = gm2fieldUtil::Graph::GetTGraph(EV,TEMP);
      gm2fieldUtil::Graph::SetGraphParameters(gFreq[i],20,kWhite+i+1);
      gm2fieldUtil::Graph::SetGraphParameters(gTemp[i],20,kWhite+i+1);
      // clean up for next time
      EV.clear();
      TIME.clear();
      FREQ.clear();
      TEMP.clear();
      X.clear();
      Y.clear();
      Z.clear();
   }
   return 0;
}
//______________________________________________________________________________
int GetTRLYSwapData(bool useOscCor,int probe,int nev,double time,
                    std::vector<averageFixedProbeEvent_t> fxpr,
                    std::vector<trolleyAnaEvent_t> Data,
                    std::vector<double> &TIME,std::vector<double> &FREQ,std::vector<double> &TEMP,
		    std::vector<double> &X,std::vector<double> &Y,std::vector<double> &Z){
   // find the mean field at the times specified in the time vector 

   int rc=0;
   // all variables 
   // rc = FilterSingle("time",probe,nev,time,Data,TIME);
   rc = FilterSingle("temp",probe,nev,time,Data,TEMP);
   rc = FilterSingle("r"   ,probe,nev,time,Data,X   );
   rc = FilterSingle("y"   ,probe,nev,time,Data,Y   );
   rc = FilterSingle("phi" ,probe,nev,time,Data,Z   );
   // now for frequency 
   // do oscillation correction and obtain ALL data associated with toggle times in time vector 
   std::vector<double> tt;
   tt.push_back(time); 
   std::vector<double> trFreq,trFreq_cor;
   rc = CorrectOscillation_trly(probe,nev,tt,fxpr,Data,TIME,trFreq,trFreq_cor);
   int M = TIME.size(); 
   for(int i=0;i<M;i++){
      if(useOscCor){
	 FREQ.push_back(trFreq_cor[i]); 
      }else{
	 FREQ.push_back(trFreq[i]); 
      }
   }

   return 0;
}
//______________________________________________________________________________
int PrintToFile(std::string outpath,std::vector<calibSwap_t> Event){

   const int N = Event.size();
   char myStr[1000];

   std::ofstream outfile;
   outfile.open( outpath.c_str() );
   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      for(int i=0;i<N;i++){
	 sprintf(myStr,"%.0lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf",
                 Event[i].time,Event[i].freq,Event[i].freqErr,Event[i].temp,Event[i].tempErr,
                 Event[i].r,Event[i].rErr,Event[i].y,Event[i].yErr,Event[i].phi,Event[i].phiErr);
         outfile << myStr << std::endl;
      }
      std::cout << "The data has been written to file: " << outpath << std::endl;
      outfile.close();
   }
   return 0;
}
