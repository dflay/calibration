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

#include "./include/date.h"
#include "./include/Constants.h"
#include "./include/fixedProbeEvent.h"
#include "./include/trolleyAnaEvent.h"
#include "./include/sccEvent.h" 

#include "./src/InputManager.C"
#include "./src/FitFuncs.C"
#include "./src/FXPRFuncs.C"
#include "./src/TRLYFuncs.C"
#include "./src/Consolidate.C"
#include "./src/CustomUtilities.C"
#include "./src/CustomMath.C"
#include "./src/CustomGraph.C"
#include "./src/CustomImport.C"
#include "./src/CustomExport.C"
#include "./src/CustomAlgorithms.C"

enum outType{
   kSWAP = 0,
   kPOS  = 1
}; 

double gMarkerSize = 0.8; 

TGraph *GetTRLYVelocityTGraph(int probe,TString xAxis,TString yAxis,
                              std::vector<gm2field::galilTrolley_t> galil,
                              std::vector<trolleyAnaEvent_t> trly); 

int GetTRLYSwapPlots(int probe,int nev,std::vector<double> time,
                     std::vector<trolleyAnaEvent_t> trlyData,TGraph **gFreq,TGraph **gTemp);

int GetTRLYSwapData(int probe,int nev,double time,std::vector<trolleyAnaEvent_t> Data,
                    std::vector<double> &TIME,std::vector<double> &FREQ,std::vector<double> &TEMP,
		    std::vector<double> &X,std::vector<double> &Y,std::vector<double> &Z); 

// int GetTRLYProbeStatsAtTime(int probe,int nev,std::vector<trolleyAnaEvent_t> trlyData,
//                             std::vector<double> time,std::vector<double> &freq,std::vector<double> &freqErr,
//                             std::vector<double> &temp,std::vector<double> &tempErr);

int PrintToFileSingle(std::string outpath,int run,double time,double t,double x,double dx,double y,double dy);

int PrintToFile(int type,std::string outpath,std::vector<trolleySwapEvent_t> Event); 
int PrintToFile(std::string outpath,std::vector<double> Time,
                std::vector<double> Freq,std::vector<double> dFreq,
                std::vector<double> Temp,std::vector<double> dTemp);

int Process_trly_prod(std::string configFile){

   int rc=0;
   int method = gm2fieldUtil::Constants::kPhaseDerivative;

   InputManager *inputMgr = new InputManager();
   inputMgr->Load(configFile);
   inputMgr->Print();

   std::string anaDate = inputMgr->GetAnalysisDate();
   bool isBlind        = inputMgr->IsBlind();
   bool loadSwapTimes  = inputMgr->GetSwapTimeStatus(); 
   int probeNumber     = inputMgr->GetTrolleyProbe(); 

   date_t theDate;
   GetDate(theDate);

   char plotDir[200];
   sprintf(plotDir,"./plots/%s",theDate.getDateString().c_str());
   rc = MakeDirectory(plotDir);

   char outDir[200];
   sprintf(outDir,"./output");
   if(isBlind)  sprintf(outDir,"%s/blinded"  ,outDir);
   if(!isBlind) sprintf(outDir,"%s/unblinded",outDir);
   sprintf(outDir,"%s/%s",outDir,theDate.getDateString().c_str());
   rc = MakeDirectory(outDir);

   blind_t blind;
   ImportBlinding(blind);
   double blindValue = blind.value_tr;

   std::vector<int> run;
   std::vector<std::string> label;
   inputMgr->GetRunList(run);
   inputMgr->GetRunLabels(label);

   std::vector<gm2field::galilTrolley_t> trlyGalil; 
   std::vector<trolleyAnaEvent_t> trlyData;   
 
   const int NRUNS = run.size();
   for(int i=0;i<NRUNS;i++){
      rc = GetTrolleyData("",run[i],method,trlyData);
      if(rc!=0){
	 std::cout << "No data!" << std::endl;
	 return 1;
      }
      rc = gm2fieldUtil::RootHelper::GetTrolleyGalil(run[i],trlyGalil);
      if(rc!=0){
	 std::cout << "No data!" << std::endl;
	 return 1;
      }
   }
   
   if(isBlind) ApplyBlindingTRLY(blindValue,trlyData);

   // TGraph *gSig  = GetTRLYVelocityTGraph(probeNumber-1,"GpsTimeStamp","vSig" ,trlyGalil,trlyData);
   // TGraph *gFish = GetTRLYVelocityTGraph(probeNumber-1,"GpsTimeStamp","vFish",trlyGalil,trlyData);
   // 
   // gm2fieldUtil::Graph::SetGraphParameters(gSig ,20,kBlue);    
   // gm2fieldUtil::Graph::SetGraphParameters(gFish,21,kRed );    

   // TMultiGraph *mgVelo = new TMultiGraph(); 
   // mgVelo->Add(gSig ,"lp"); 
   // mgVelo->Add(gFish,"lp"); 

   std::vector<double> time,freq,freqErr,temp,tempErr;

   // find the times for each TRLY swap
   double angle = 189.160; // nominal trolley location 

   if(loadSwapTimes){
      LoadTRLYTimes(probeNumber,"swap",time);
   }else{
      FindTRLYStopTimes(probeNumber-1,angle,trlyData,trlyGalil,time);
   }

   const int NL = time.size();
   std::cout << "Found " << NL << " trolley swaps" << std::endl;

   // now get average field for 30 seconds BEFORE each time stamp
   int nev = 30;
   double fLO = 61.74E+6; 

   std::vector<trolleySwapEvent_t> trlySwap; 
   rc = GetTRLYStatsAtTime(probeNumber-1,nev,fLO,time,trlyData,trlySwap);

   char outPath[500];
   sprintf(outPath,"%s/trly-swap-data_pr-%02d_%s.csv",outDir,probeNumber,anaDate.c_str());
   std::string outpath_swap = outPath;

   sprintf(outPath,"%s/trly-swap-pos_pr-%02d_%s.csv",outDir,probeNumber,anaDate.c_str());
   std::string outpath_pos = outPath;

   rc = PrintToFile(kSWAP,outpath_swap,trlySwap); 
   rc = PrintToFile(kPOS ,outpath_pos ,trlySwap); 

   // make plots 
   TGraph *gTR = GetTRLYTGraph(probeNumber-1,"GpsTimeStamp","freq",trlyData);
   gm2fieldUtil::Graph::SetGraphParameters(gTR,20,kBlack); 

   TGraph **gFreq = new TGraph*[NL]; 
   TGraph **gTemp = new TGraph*[NL];
   rc = GetTRLYSwapPlots(probeNumber-1,nev,time,trlyData,gFreq,gTemp); 

   TMultiGraph *mgf = new TMultiGraph();
   TMultiGraph *mgt = new TMultiGraph();
   for(int i=0;i<NL;i++){
      mgf->Add(gFreq[i],"lp");
      mgt->Add(gTemp[i],"lp");
   }

   double mean  = gm2fieldUtil::Math::GetMean<double>(freq); 
   double yMin = (mean-fLO)-1E+3; 
   double yMax = (mean-fLO)+1E+3;

   double stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(freq); 
   double yMin2 = (mean-fLO)-10.*stdev;
   double yMax2 = (mean-fLO)+10.*stdev;

   double mean_t = gm2fieldUtil::Math::GetMean<double>(temp); 
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
   // gm2fieldUtil::Graph::UseTimeDisplay(mgf); 
   gm2fieldUtil::Graph::SetGraphLabelSizes(mgf,0.05,0.06);
   mgf->GetYaxis()->SetRangeUser(yMin2,yMax2);
   mgf->Draw("a");
   // L->Draw("same"); 
   c1->Update();

   std::cout << __LINE__ << std::endl;

   c1->cd(3);
   mgt->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mgt,Form("TRLY-%02d Swap Events",probeNumber),"","Temperature (#circC)");
   // gm2fieldUtil::Graph::UseTimeDisplay(mgt); 
   gm2fieldUtil::Graph::SetGraphLabelSizes(mgt,0.05,0.06);
   mgt->GetYaxis()->SetRangeUser(yMin3,yMax3);
   mgt->Draw("a");
   c1->Update();
   
   TString plotPath = Form("%s/trly-swap_pr-%02d.png",plotDir,probeNumber);
   c1->cd(); 
   std::cout << __LINE__ << std::endl;
   c1->Print(plotPath); 
   std::cout << __LINE__ << std::endl;
   delete c1; 

   // TCanvas *c2 = new TCanvas("c2","TRLY Velocities",1200,600); 
   // c2->cd();

   // c1->cd(2); 
   // mgVelo->Draw("a"); 
   // gm2fieldUtil::Graph::SetGraphLabels(mgVelo,"TRLY Velocities (blue = signal, red = fishing)","","Velocity (counts/sec)"); 
   // gm2fieldUtil::Graph::UseTimeDisplay(mgVelo); 
   // gm2fieldUtil::Graph::SetGraphLabelSizes(mgVelo,0.05,0.06); 
   // mgVelo->Draw("a");
   // c1->Update();


   std::cout << __LINE__ << std::endl;

   std::cout << "--> Done." << std::endl; 

   // plotPath = Form("%s/trly-velo_run-%05d.png",plotDir,run[0]);
   // c2->cd(); 
   // c2->Print(plotPath); 

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
int GetTRLYSwapPlots(int probe,int nev,std::vector<double> time,
                     std::vector<trolleyAnaEvent_t> trlyData,TGraph **gFreq,TGraph **gTemp){
   // for diagnostic plots
   int M=0,rc=0;
   const int NL = time.size();
   std::vector<double> EV,TIME,FREQ,TEMP,X,Y,Z;
   for(int i=0;i<NL;i++){
      rc = GetTRLYSwapData(probe,nev,time[i],trlyData,TIME,FREQ,TEMP,X,Y,Z);
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
int GetTRLYSwapData(int probe,int nev,double time,std::vector<trolleyAnaEvent_t> Data,
                    std::vector<double> &TIME,std::vector<double> &FREQ,std::vector<double> &TEMP,
		    std::vector<double> &X,std::vector<double> &Y,std::vector<double> &Z){
   // find the mean field at the times specified in the time vector 
   int rc=0;
   rc = FilterSingle("time",probe,nev,time,Data,TIME);
   rc = FilterSingle("freq",probe,nev,time,Data,FREQ);
   rc = FilterSingle("temp",probe,nev,time,Data,TEMP);
   rc = FilterSingle("r"   ,probe,nev,time,Data,X   );
   rc = FilterSingle("y"   ,probe,nev,time,Data,Y   );
   rc = FilterSingle("phi" ,probe,nev,time,Data,Z   );
   return 0;
}
//______________________________________________________________________________
int PrintToFile(int type,std::string outpath,std::vector<trolleySwapEvent_t> Event){

   const int N = Event.size();
   char myStr[1000];

   std::ofstream outfile;
   outfile.open( outpath.c_str() );
   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      for(int i=0;i<N;i++){
	 if(type==kSWAP){
	    sprintf(myStr,"%.0lf,%.3lf,%.3lf,%.3lf,%.3lf",Event[i].time,Event[i].freq,Event[i].freq_err,Event[i].temp,Event[i].temp_err);
	 }else if(type==kPOS){
	    sprintf(myStr,"%.0lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf",
                    Event[i].time,Event[i].r,Event[i].r_err,Event[i].y,Event[i].y_err,Event[i].phi,Event[i].phi_err);
	 }
         outfile << myStr << std::endl;
      }
      std::cout << "The data has been written to file: " << outpath << std::endl;
      outfile.close();
   }
   return 0;
}
//______________________________________________________________________________
int PrintToFile(std::string outpath,std::vector<double> Time,
                std::vector<double> Freq,std::vector<double> dFreq,
                std::vector<double> Temp,std::vector<double> dTemp){

   const int N = Time.size();
   char myStr[1000];

   std::ofstream outfile;
   outfile.open( outpath.c_str() );
   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      for(int i=0;i<N;i++){
         sprintf(myStr,"%.0lf,%.3lf,%.3lf,%.3lf,%.3lf",Time[i],Freq[i],dFreq[i],Temp[i],dTemp[i]);
         outfile << myStr << std::endl;
      }
      std::cout << "The data has been written to file: " << outpath << std::endl;
      outfile.close();
   }
   return 0;
}