// Compute DeltaB for TRLY data sets  

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
#include "./include/plungingProbeAnaEvent.h"
#include "./include/trolleyAnaEvent.h"

#include "./src/FitFuncs.C"
#include "./src/FXPRFuncs.C"
#include "./src/Consolidate.C"
#include "./src/CustomMath.C"
#include "./src/CustomGraph.C"
#include "./src/CustomImport.C"
#include "./src/CustomExport.C"
#include "./src/CustomAlgorithms.C"
#include "./src/DeltaBFuncs.C"

int PlotTRLY(){

   int rc=0;
   int method = gm2fieldUtil::Constants::kPhaseDerivative;

   int run=0;
   std::cout << "Enter run: ";
   std::cin  >> run; 

   // fixed probe data
   std::string fxprPath = "./input/probe-lists/fxpr-list.csv"; 
   std::vector<int> fxprList; 
   gm2fieldUtil::Import::ImportData1<int>(fxprPath,"csv",fxprList);

   // get data for drift correction *during* a run 
   std::vector<gm2field::fixedProbeFrequency_t> fxprData;
   rc = gm2fieldUtil::RootHelper::GetFPFrequencies(run,fxprData);
   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   // Trolley data 
   std::vector<trolleyAnaEvent_t> trlyData,trlyDataCor;   

   rc = GetTrolleyData("",run,method,trlyData);
   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   std::cout << "Applying field drift corrections (during measurement)..." << std::endl;
   // correct for field drift [fxpr]  
   rc = CorrectTRLYForDriftDuringMeasurement(method,fxprList,fxprData,trlyData,trlyDataCor);
   if(rc!=0) return 1;
   std::cout << "--> Done" << std::endl; 

   std::vector<double> fb,fa;
   const int NB = trlyData.size();
   for(int i=0;i<NB;i++) fb.push_back( trlyData[i].freq[0] );

   const int NA = trlyData.size();
   for(int i=0;i<NA;i++) fa.push_back( trlyDataCor[i].freq[0] );

   double mean_before  = gm2fieldUtil::Math::GetMean<double>(fb); 
   double stdev_before = gm2fieldUtil::Math::GetStandardDeviation<double>(fb); 

   double mean_after  = gm2fieldUtil::Math::GetMean<double>(fa); 
   double stdev_after = gm2fieldUtil::Math::GetStandardDeviation<double>(fa); 

   std::cout << Form("mean before = %.3lf +/- %.3lf Hz",mean_before,stdev_before) << std::endl; 
   std::cout << Form("mean after  = %.3lf +/- %.3lf Hz",mean_after ,stdev_after)  << std::endl; 

   TGraph *gBefore = GetTRLYTGraph(0,"GpsTimeStamp","freq",trlyData);
   gm2fieldUtil::Graph::SetGraphParameters(gBefore,20,kBlack);
 
   TGraph *gAfter  = GetTRLYTGraph(0,"GpsTimeStamp","freq",trlyDataCor);
   gm2fieldUtil::Graph::SetGraphParameters(gAfter,20,kRed);

   const int NFXPR = fxprList.size();
   int startIndex  = fxprList[0];
   int stopIndex   = fxprList[NFXPR-1];

   // Fixed probe average plot 
   const int NFP = fxprData.size();
   TGraph *gFPAVG = GetTGraph(method,fxprData[0].GpsTimeStamp[startIndex],fxprData[NFP-1].GpsTimeStamp[stopIndex],1E+9,fxprList,fxprData);
   gm2fieldUtil::Graph::SetGraphParameters(gFPAVG,20,kBlack); 

   const int NPAR = 6;
   TF1 *fxprFit = GetFitToFXPR("myFit",MyFitFunc_pol4,NPAR,method,fxprList,fxprData);

   TMultiGraph *mg = new TMultiGraph();
   mg->Add(gBefore,"lp"); 
   mg->Add(gAfter ,"lp");

   TLegend *L = new TLegend(0.6,0.6,0.8,0.8);
   L->AddEntry(gBefore,"Before Drift Correction","p"); 
   L->AddEntry(gAfter ,"After Drift Correction" ,"p");

   TMultiGraph *mg2 = new TMultiGraph();
   TLegend *L2      = new TLegend(0.6,0.6,0.8,0.8); 

   gm2fieldUtil::Graph::FillMultiGraph(fxprList,method,"GpsTimeStamp","Frequency",fxprData,mg2,L2);

   TCanvas *c1 = new TCanvas("c1","TRLY Data",1200,600);
   c1->Divide(1,2);  
   
   c1->cd(1);
   mg->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mg,"TRLY Data","","Frequency (Hz)");
   gm2fieldUtil::Graph::SetGraphLabelSizes(mg,0.05,0.06); 
   gm2fieldUtil::Graph::UseTimeDisplay(mg); 
   mg->Draw("a");
   L->Draw("same"); 
   c1->Update();

   c1->cd(2);
   gFPAVG->Draw("alp"); 
   gm2fieldUtil::Graph::SetGraphLabels(gFPAVG,"FXPR Avg","","Frequency (Hz)"); 
   gm2fieldUtil::Graph::SetGraphLabelSizes(gFPAVG,0.05,0.06); 
   gm2fieldUtil::Graph::UseTimeDisplay(gFPAVG); 
   gFPAVG->Draw("alp"); 
   fxprFit->Draw("same"); 
   c1->Update();

   // c1->cd(3);
   // mg2->Draw("a"); 
   // gm2fieldUtil::Graph::SetGraphLabels(mg2,"Fixed Probes","","Frequency (Hz)"); 
   // gm2fieldUtil::Graph::SetGraphLabelSizes(mg2,0.05,0.06); 
   // gm2fieldUtil::Graph::UseTimeDisplay(mg2); 
   // mg2->Draw("a");
   // L2->Draw("same");  
   // c1->Update();

   return 0;
}

