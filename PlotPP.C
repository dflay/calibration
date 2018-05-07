// Plot the Plunging Probe data   

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

#include "./include/Constants.h"
#include "./include/drift.h"
#include "./include/fixedProbeEvent.h"
#include "./include/plungingProbeAnaEvent.h"
#include "./include/trolleyAnaEvent.h"

#include "./src/FXPRFuncs.C"
#include "./src/Consolidate.C"
#include "./src/CustomMath.C"
#include "./src/CustomGraph.C"
#include "./src/CustomImport.C"
#include "./src/CustomExport.C"
#include "./src/CustomAlgorithms.C"
#include "./src/DeltaBFuncs.C"

int PlotPP(){

   int rc=0;
   int method = gm2fieldUtil::Constants::kPhaseDerivative;

   int runNumber=0;

   std::cout << "Enter run number: ";
   std::cin  >> runNumber; 

   // PP data 
   std::vector<plungingProbeAnaEvent_t> ppEvent,ppEventCor; 
   std::cout << "Getting run " << runNumber << std::endl; 
   rc = GetPlungingProbeData(runNumber,method,ppEvent);
   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   const int NPP = ppEvent.size();
   for(int i=0;i<NPP;i++) rc = ModifyPlungingProbeData(kLeastSquaresPhase,ppEvent[i]);

   // fixed probe data
   std::string fxprPath = "./input/probe-lists/fxpr-list.csv";
   std::vector<int> fxprList;
   gm2fieldUtil::Import::ImportData1<int>(fxprPath,"csv",fxprList);

   std::vector<gm2field::fixedProbeFrequency_t> fxprData;
   rc = gm2fieldUtil::RootHelper::GetFPFrequencies(runNumber,fxprData);
   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   const int NN = fxprData.size();
   const int NFP = fxprList.size(); 
   unsigned long long t0     = 0;
   unsigned long long tStart = fxprData[0].GpsTimeStamp[fxprList[0]];  
   unsigned long long tStop  = fxprData[NN-1].GpsTimeStamp[fxprList[NFP-1]]; 
   unsigned long long tStep  = 1E+9;  

   std::vector<fixedProbeEvent_t> fxprDataAvg; 
   GetAverageFXPRVectorsNew(method,t0,tStart,tStop,tStep,fxprList,fxprData,fxprDataAvg);

   rc = CorrectPPForDriftDuringMeasurement(fxprDataAvg,ppEvent,ppEventCor);

   std::vector<double> fb,fa;
   const int NB = ppEvent[0].numTraces;  
   for(int i=0;i<NB;i++) fb.push_back( ppEvent[0].freq[i] ); 
   double mean_before  = gm2fieldUtil::Math::GetMean<double>(fb); 
   double stdev_before = gm2fieldUtil::Math::GetStandardDeviation<double>(fb); 

   const int NA = ppEventCor[0].numTraces;  
   for(int i=0;i<NA;i++) fa.push_back( ppEventCor[0].freq[i] ); 
   double mean_after  = gm2fieldUtil::Math::GetMean<double>(fa); 
   double stdev_after = gm2fieldUtil::Math::GetStandardDeviation<double>(fa); 
  
   std::cout << Form("mean (before) = %.3lf +/- %.3lf Hz",mean_before,stdev_before) << std::endl; 
   std::cout << Form("mean (after)  = %.3lf +/- %.3lf Hz",mean_after ,stdev_after)  << std::endl; 

   std::cout << "Getting PP plots..." << std::endl;

   // Plunging probe plots 
   TGraph *g1 = GetPPTGraph1("TimeStamp","freq",ppEvent);
   gm2fieldUtil::Graph::SetGraphParameters(g1,20,kBlack);

   TGraph *g2 = GetPPTGraph1("TimeStamp","freq",ppEventCor);
   gm2fieldUtil::Graph::SetGraphParameters(g2,20,kRed);
 
   std::cout << "--> Done." << std::endl;

   TMultiGraph *mg = new TMultiGraph();
   mg->Add(g1,"lp"); 
   mg->Add(g2,"lp");

   std::cout << "Getting FXPR plot..." << std::endl; 

   // Fixed probe plot
   std::cout << "FXPR start = " << gm2fieldUtil::GetStringTimeStampFromUTC(tStart/1E+9) << std::endl; 
   std::cout << "FXPR end   = " << gm2fieldUtil::GetStringTimeStampFromUTC(tStop/1E+9)  << std::endl; 
 
   TGraph *gFXPR = GetTGraphNew(fxprDataAvg);
   gm2fieldUtil::Graph::SetGraphParameters(gFXPR,20,kBlack);  
   
   std::cout << "--> Done." << std::endl;

   double mean_during=0,stdev=0;
   std::vector<double> X,Y; 

   TLegend *L = new TLegend(0.6,0.6,0.8,0.8); 
   L->AddEntry(g1,"Before Drift Correction","p");  
   L->AddEntry(g2,"After Drift Correction","p");  
   
   TString Title = Form("PP Data [Run %d]",runNumber);

   TMultiGraph *mg2 = new TMultiGraph();
   mg2->Add(gFXPR,"lp"); 

   double xMin = tStart/1E+9;  
   double xMax = tStop/1E+9;  

   TCanvas *c1 = new TCanvas("c1","PP & FXPR Data",1200,600);
   c1->Divide(1,2);
  
   c1->cd(1); 
   mg->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mg,Title,"","Frequency (Hz)"); 
   gm2fieldUtil::Graph::UseTimeDisplay(mg);  
   gm2fieldUtil::Graph::SetGraphLabelSizes(mg,0.05,0.06); 
   mg->GetXaxis()->SetLimits(xMin,xMax);
   mg->Draw("a");
   L->Draw("same"); 
   c1->Update(); 

   c1->cd(2); 
   mg2->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mg2,"Fixed Probe Average","","Frequency (Hz)"); 
   gm2fieldUtil::Graph::UseTimeDisplay(mg2);  
   gm2fieldUtil::Graph::SetGraphLabelSizes(mg2,0.05,0.06); 
   mg2->GetXaxis()->SetLimits(xMin,xMax);
   mg2->Draw("a");
   // gTest->Draw("same"); 
   c1->Update(); 

   return 0;
}

