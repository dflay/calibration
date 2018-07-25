// Process the PP data; apply a drift correction during the run, 
// produce an output file based on the input run list containing 
// the raw frequencies and the drift-corrected frequencies for all probes  

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

double gMarkerSize = 0.8; 

int GetStatsForPP(std::vector<plungingProbeAnaEvent_t> data,double &mean,double &stdev,double &temp,double &min,double &max); 

int PrintToFile(std::string date,int run,double time,double temp,double x,double dx,double y,double dy);

int ProcessPP(std::string configFile){

   int rc=0;
   int method = gm2fieldUtil::Constants::kPhaseDerivative;

   InputManager *inputMgr = new InputManager();
   inputMgr->Load(configFile);
   inputMgr->Print();

   std::string anaDate = inputMgr->GetAnalysisDate();
   bool isBlind        = inputMgr->IsBlind();
   int fxpr_tag        = inputMgr->GetFixedProbeListTag(); 

   date_t theDate;
   GetDate(theDate);

   char plotDir[200];
   sprintf(plotDir,"./plots/%02d-%02d-%02d",theDate.month,theDate.day,theDate.year-2000);
   rc = MakeDirectory(plotDir);

   char outDir[200];
   sprintf(outDir,"./output"); 
   if(isBlind)  sprintf(outDir,"%s/blinded"  ,outDir);
   if(!isBlind) sprintf(outDir,"%s/unblinded",outDir);
   sprintf(outDir,"%s/%02d-%02d-%02d",outDir,theDate.month,theDate.day,theDate.year-2000); 
   rc = MakeDirectory(outDir);

   char outPath[500]; 
   sprintf(outPath,"%s/pp-data_%s.csv",outDir,anaDate.c_str());
   std::string outpath = outPath;  

   blind_t blind;
   ImportBlinding(blind);
   double blindValue = blind.value_pp;

   std::vector<int> run;
   std::vector<std::string> label;
   inputMgr->GetRunList(run);
   inputMgr->GetRunLabels(label);

   const int NRUNS = run.size();
   
   char fxpr_path[500]; 
   sprintf(fxpr_path,"./input/probe-lists/fxpr-list_set-%d.csv",fxpr_tag); 
   std::string fxprPath = fxpr_path;
   std::vector<int> fxprList;
   gm2fieldUtil::Import::ImportData1<int>(fxprPath,"csv",fxprList);

   std::vector<gm2field::fixedProbeFrequency_t> fxprData;

   // details for start and stop indices for fxpr data  
   int NN=0;
   int NFP = fxprList.size();
   int startIndex = fxprList[0];
   int stopIndex  = fxprList[NFP-1];
   unsigned long long t0     = 0;
   unsigned long long tStart=0,tStop=0,tStep=1E+9;

   std::vector<fixedProbeEvent_t> fxprDataAvg;

   TGraph **gFPAVG = new TGraph*[NRUNS];
   TString fxprPlots = Form("%s/fxpr_runs-%04d-%04d.png",plotDir,run[0],run[NRUNS-1]);

   double yMin=0,yMax=0;

   double time,freq,freqErr,freq_cor,freqErr_cor,temp,temp_cor; 
   std::vector<double> FREQ,FREQERR,FREQ_cor,FREQERR_cor; 
   std::vector<trolleyAnaEvent_t> trlyData,trlyDataCor;     

   TCanvas *c1 = new TCanvas("c1","TRLY plots",1200,1000);
   c1->Divide( (NRUNS+1)/2,2 );

   TCanvas *c2 = new TCanvas("c2","FXPR plots",1200,1000);
   c2->Divide((NRUNS+1)/2,2);

   TString ppPlots;  

   // PP data 
   std::vector<plungingProbeAnaEvent_t> ppData,ppDataCor;

   // main analysis loop 
   for(int i=0;i<NRUNS;i++){
      std::cout << "------------------ RUN " << run[i] << " ------------------" << std::endl;
      std::cout << "Getting fixed probe data..." << std::endl;
      rc = gm2fieldUtil::RootHelper::GetFPFrequencies(run[i],fxprData);
      NN = fxprData.size();
      tStart = fxprData[0].GpsTimeStamp[startIndex];
      tStop  = fxprData[NN-1].GpsTimeStamp[stopIndex];
      GetAverageFXPRVectorsNew(method,t0,tStart,tStop,tStep,fxprList,fxprData,fxprDataAvg);
      gFPAVG[i] = GetTGraphNew(fxprDataAvg);
      gm2fieldUtil::Graph::SetGraphParameters(gFPAVG[i],20,kBlack);
      std::cout << "Getting PP data..." << std::endl;
      rc = GetPlungingProbeData(run[i],method,ppData);
      if(rc!=0){
	 std::cout << "No data!" << std::endl;
	 return 1;
      }
      if(isBlind) ApplyBlindingPP(blindValue,ppData);
      std::cout << "--> Done." << std::endl;
      // apply P2P drift correction
      std::cout << "Applying drift correction..." << std::endl;
      rc = CorrectPPForDriftDuringMeasurement(fxprDataAvg,ppData,ppDataCor);
      std::cout << "--> Done." << std::endl;
      // get stats for the PP 
      time = ppData[0].time[0]/1E+9; // time defined as that of the first trace
      rc = GetStatsForPP(ppData   ,freq    ,freqErr    ,temp    ,yMin,yMax); 
      rc = GetStatsForPP(ppDataCor,freq_cor,freqErr_cor,temp_cor,yMin,yMax); 
      // make plots
      TGraph *gPP     = GetPPTGraph1("TimeStamp","freq",ppData); 
      TGraph *gPP_cor = GetPPTGraph1("TimeStamp","freq",ppDataCor);
      gm2fieldUtil::Graph::SetGraphParameters(gPP    ,21,kBlack);
      gm2fieldUtil::Graph::SetGraphParameters(gPP_cor,20,kRed);
      gPP->SetMarkerSize(gMarkerSize); 
      gPP_cor->SetMarkerSize(gMarkerSize);
      // put in canvas 
      c1->cd(i+1); 
      gPP->SetTitle( Form("Run %d",run[i]) );
      gPP->Draw("alp");  
      gPP->GetYaxis()->SetRangeUser(yMin,yMax);
      gm2fieldUtil::Graph::UseTimeDisplay(gPP); 
      gPP->Draw("alp");  
      gPP_cor->Draw("same"); 
      c1->Update();
      // fxpr plots 
      c2->cd(i+1);
      gFPAVG[i]->Draw("alp");
      gm2fieldUtil::Graph::SetGraphLabels(gFPAVG[i],Form("FXPR Data Run %d",run[i]),"","Frequency (Hz)");
      gm2fieldUtil::Graph::UseTimeDisplay(gFPAVG[i]);
      gFPAVG[i]->Draw("alp");
      c2->Update(); 
      // print to file 
      rc = PrintToFile(outpath,run[i],time,temp,freq,freqErr,freq_cor,freqErr_cor); 
      // set up for next run
      ppData.clear();
      ppDataCor.clear(); 
      fxprData.clear(); 
      fxprDataAvg.clear(); 
   }

   // save the plot  
   c1->cd(); 
   ppPlots = Form("%s/pp-data.png",plotDir);
   c1->Print(ppPlots);

   c2->cd();
   c2->Print(fxprPlots);

   return 0;
}
//______________________________________________________________________________
int GetStatsForPP(std::vector<plungingProbeAnaEvent_t> data,double &mean,double &stdev,double &temp,double &min,double &max){
   min = 50E+3; 
   max = 0; 
   // gather stats for the PP; recall a single run has N traces in it
   double arg=0;
   int M=0;
   const int N = data.size();
   std::vector<double> x,xt; 
   for(int i=0;i<N;i++){
      M = data[i].numTraces;
      for(int j=0;j<M;j++){
         // we want to add back in the LO since we're comparing against the TRLY eventually 
	 arg = data[i].freq[j] + data[i].freq_LO[j]; 
	 x.push_back(arg);
         xt.push_back(data[i].temp[j]); 
         if( min>data[i].freq[j] ) min = data[i].freq[j];  
         if( max<data[i].freq[j] ) max = data[i].freq[j];  
      }
   }
   mean  = gm2fieldUtil::Math::GetMean<double>(x); 
   stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(x);
   temp  = gm2fieldUtil::Math::GetMean<double>(xt);  

   min -= 2; 
   max += 2;

   return 0;
}
//______________________________________________________________________________
int PrintToFile(std::string outpath,int run,double time,double temp,double x,double dx,double y,double dy){
   // x = uncorrected, y = drift corrected 
   char myStr[1000]; 
   // prepare the line
   sprintf(myStr,"%d,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf",run,time,temp,x,dx,y,dy);

   std::ofstream outfile;
   outfile.open(outpath.c_str(),ios::app);
   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      outfile << myStr << std::endl;
      std::cout << "The data has been written to file: " << outpath << std::endl;
   }  
   return 0;
}
