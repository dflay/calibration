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

int PlotFXPR(){

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

   // Fixed probe average plot 
   TGraph *gFPAVG = GetTGraphNew(fxprDataAvg);
   gm2fieldUtil::Graph::SetGraphParameters(gFPAVG,20,kBlack); 

   const int nOrder = 2;
   std::vector< std::vector<double> > eps;
   TF1 *fxprFit = GetPolyFitToFXPR("myFit",nOrder,gFPAVG,eps);

   const int NPAR = fxprFit->GetNpar(); 
   double par[NPAR],parErr[NPAR];
   for(int i=0;i<NPAR;i++){
      par[i]    = fxprFit->GetParameter(i); 
      parErr[i] = fxprFit->GetParError(i);
      std::cout << Form("p[%d] = %.3lf +/- %.3lf",i,par[i],parErr[i]) << std::endl; 
   }

   // get a graph of all probes 
   // TMultiGraph *mg = new TMultiGraph(); 
   // TLegend *L      = new TLegend(0.6,0.6,0.8,0.8); 
   // rc = gm2fieldUtil::Graph::FillMultiGraph(fxprList,method,"GpsTimeStamp","Frequency",fxprData,mg,L);

   TCanvas *c1 = new TCanvas("c1","FXPR Data",1200,600);
   // c1->Divide(1,2);  
   
   // c1->cd(1);
   // mg->Draw("al");
   // gm2fieldUtil::Graph::SetGraphLabels(mg,"Individual Fixed Probes","","Frequency (Hz)");
   // gm2fieldUtil::Graph::SetGraphLabelSizes(mg,0.05,0.06); 
   // gm2fieldUtil::Graph::UseTimeDisplay(mg); 
   // mg->Draw("al");
   // L->Draw("same"); 
   // c1->Update();

   c1->cd();
   gFPAVG->Draw("alp"); 
   gm2fieldUtil::Graph::SetGraphLabels(gFPAVG,"FXPR Avg","","Frequency (Hz)"); 
   // gm2fieldUtil::Graph::SetGraphLabelSizes(gFPAVG,0.05,0.06); 
   gm2fieldUtil::Graph::UseTimeDisplay(gFPAVG); 
   gFPAVG->Draw("alp");
   // gFPAVG->Fit("pol1");  
   // fxprFit->Draw("same"); 
   c1->Update();

   return 0;
}

