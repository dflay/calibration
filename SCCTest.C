// Test out the FXPR plotting functions   

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

#include "./src/FitFuncs.C"
#include "./src/FXPRFuncs.C"
#include "./src/Consolidate.C"
#include "./src/CustomMath.C"
#include "./src/CustomGraph.C"
#include "./src/CustomImport.C"
#include "./src/CustomExport.C"
#include "./src/CustomAlgorithms.C"
#include "./src/DeltaBFuncs.C"

int SCCTest(){

   int rc=0;
   int method = gm2fieldUtil::Constants::kPhaseDerivative;

   std::vector<int> run;
   gm2fieldUtil::Import::GetRunList(run);
   const int NRUNS = run.size();

   // fixed probe data
   std::string fxprPath = "./input/probe-lists/fxpr-list.csv"; 
   std::vector<int> fxprList; 
   gm2fieldUtil::Import::ImportData1<int>(fxprPath,"csv",fxprList);

   // get fixed probe data 
   std::vector<gm2field::fixedProbeFrequency_t> fxprData;
   for(int i=0;i<NRUNS;i++) rc = gm2fieldUtil::RootHelper::GetFPFrequencies(run[i],fxprData);
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

   double thr = 30.; 
   std::vector<double> tLo,tHi; 
   rc = FindTransitionTimes(thr,fxprDataAvg,tLo,tHi);

   double yMin = 55E+3; 
   double yMax = 60E+3;

   const int NL = tLo.size();
   TLine **LO = new TLine*[NL];
   for(int i=0;i<NL;i++){
      LO[i] = new TLine(tLo[i],yMin,tLo[i],yMax);
      LO[i]->SetLineColor(kGreen+1); 
      LO[i]->SetLineWidth(2); 
      LO[i]->SetLineStyle(2); 
   }   

   const int NH = tHi.size();
   TLine **HI= new TLine*[NH];
   for(int i=0;i<NH;i++){
      HI[i] = new TLine(tHi[i],yMin,tHi[i],yMax);
      HI[i]->SetLineColor(kRed); 
      HI[i]->SetLineWidth(2); 
      HI[i]->SetLineStyle(2); 
   }

   int nev=30;
   std::vector<fixedProbeEvent_t> lo1,hi1; 
   std::vector<fixedProbeEvent_t> lo2,hi2; 
   std::vector<fixedProbeEvent_t> lo3,hi3; 
   rc = FilterSingle(nev,tLo[0],fxprDataAvg,lo1);
   rc = FilterSingle(nev,tHi[0],fxprDataAvg,hi1);
   rc = FilterSingle(nev,tLo[1],fxprDataAvg,lo2);
   rc = FilterSingle(nev,tHi[1],fxprDataAvg,hi2);
   rc = FilterSingle(nev,tLo[2],fxprDataAvg,lo3);
   rc = FilterSingle(nev,tHi[2],fxprDataAvg,hi3);
  
   // Fixed probe average plot 
   TGraph *gFPAVG = GetTGraphNew(fxprDataAvg);
   gm2fieldUtil::Graph::SetGraphParameters(gFPAVG,21,kBlack); 

   TGraph *gLo1 = GetTGraphNew(lo1);
   gm2fieldUtil::Graph::SetGraphParameters(gLo1,20,kGreen+2); 

   TGraph *gHi1 = GetTGraphNew(hi1);
   gm2fieldUtil::Graph::SetGraphParameters(gHi1,20,kRed+2); 

   TGraph *gLo2 = GetTGraphNew(lo2);
   gm2fieldUtil::Graph::SetGraphParameters(gLo2,20,kGreen+2); 

   TGraph *gHi2 = GetTGraphNew(hi2);
   gm2fieldUtil::Graph::SetGraphParameters(gHi2,20,kRed+2);
 
   TGraph *gLo3 = GetTGraphNew(lo3);
   gm2fieldUtil::Graph::SetGraphParameters(gLo3,20,kGreen+2); 

   TGraph *gHi3 = GetTGraphNew(hi3);
   gm2fieldUtil::Graph::SetGraphParameters(gHi3,20,kRed+2); 

   TMultiGraph *mg = new TMultiGraph();
   mg->Add(gFPAVG,"lp"); 
   mg->Add(gLo1  ,"lp"); 
   mg->Add(gHi1  ,"lp"); 
   mg->Add(gLo2  ,"lp"); 
   mg->Add(gHi2  ,"lp"); 
   mg->Add(gLo3  ,"lp"); 
   mg->Add(gHi3  ,"lp"); 

   TCanvas *c1 = new TCanvas("c1","FXPR Data",1200,600);

   c1->cd();
   mg->Draw("a"); 
   gm2fieldUtil::Graph::SetGraphLabels(mg,"FXPR Avg","","Frequency (Hz)"); 
   // gm2fieldUtil::Graph::SetGraphLabelSizes(gFPAVG,0.05,0.06); 
   gm2fieldUtil::Graph::UseTimeDisplay(mg); 
   mg->GetYaxis()->SetRangeUser(yMin,yMax); 
   mg->Draw("a");
   for(int i=0;i<NL;i++) LO[i]->Draw("same"); 
   for(int i=0;i<NH;i++) HI[i]->Draw("same"); 
   c1->Update();

   return 0;
}

