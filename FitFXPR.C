// Compute DeltaB for TRLY data sets  

#include <cstdlib> 
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>  
#include <cmath> 

#include "TLine.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSpline.h"
#include "TROOT.h"

#include "RootTreeStructs.h"
#include "gm2fieldMath.h"
#include "gm2fieldGraph.h"
#include "gm2fieldRootHelper.h"
#include "gm2fieldUnits.h"
#include "gm2fieldImport.h"
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

int FitFXPR(){

   int rc=0;
   int method = gm2fieldUtil::Constants::kPhaseDerivative;

   std::cout << "----------------------------------" << std::endl;
   std::cout << "Fit Test " << std::endl;

   std::vector<int> run;
   rc = gm2fieldUtil::Import::GetRunList(run);
   const int NRUNS = run.size(); 
   std::cout << "Processing runs ";
   for(int i=0;i<NRUNS;i++) std::cout << run[i] << " ";
   std::cout << std::endl;

   // Fixed Probe data
   std::string fxprPath = "./input/probe-lists/fxpr-list.csv";
   std::vector<int> fxprList;
   gm2fieldUtil::Import::ImportData1<int>(fxprPath,"csv",fxprList);

   std::vector<unsigned long long int> runStart,runStop;
   std::vector<gm2field::fixedProbeFrequency_t> fxprData;
   for(int i=0;i<NRUNS;i++){
      rc = gm2fieldUtil::RootHelper::GetFPFrequencies(run[i],fxprData);
   }

   if (rc!=0) {
      std::cout << "No data.  Exiting..." << std::endl;
      return 1;
   }
 
   const int N     = fxprData.size();
   const int NFXPR = fxprList.size();
   int startIndex  = fxprList[0]; 
   int stopIndex   = fxprList[NFXPR-1];

   unsigned long long timeStart = fxprData[0].GpsTimeStamp[startIndex]; 
   unsigned long long timeStop  = fxprData[N-1].GpsTimeStamp[stopIndex];
   unsigned long long timeStep  = 1E+9; 

   const int NPTS = (timeStop-timeStart)/timeStep;
   double arg=0;
   std::vector<double> TIME;
   for(int i=0;i<NPTS;i++){
      arg = (timeStart + ( (double)i )*timeStep )/1E+9;
      TIME.push_back(arg);
   }
 
   // Fixed probe average plot 
   TGraph *gFPAVG = GetTGraphNew(method,timeStart,timeStop,timeStep,fxprList,fxprData);
   gm2fieldUtil::Graph::SetGraphParameters(gFPAVG,20,kBlack);

   TSpline3 *gSp = new TSpline3("mySpline",gFPAVG);
   gSp->SetLineColor(kBlue);
   gSp->SetLineWidth(2);
   
   const int nOrder = 1;
   std::vector< std::vector<double> > eps;  
   TF1 *myFit = GetPolyFitToFXPR("myFit",nOrder,gFPAVG,eps);
   TGraphErrors *myErrorBand = GetFitErrorBand(TIME,myFit,eps,1);

   TString fitLabel; 
   if(nOrder==3){
      fitLabel = Form("p_{0} + p_{1}x + p_{2}x^{2} + p_{3}x^{3}");
   }else if(nOrder==4){
      fitLabel = Form("p_{0} + p_{1}x + p_{2}x^{2} + p_{3}x^{3} + p_{4}x^{4}");
   }

   TLegend *L = new TLegend(0.6,0.6,0.8,0.8);
   L->AddEntry(gFPAVG,"Fixed Probe Data","p");
   L->AddEntry(myFit ,fitLabel,"l");
   L->AddEntry(gSp   ,"3^{rd} Order Spline","l");
 
   TCanvas *c1 = new TCanvas("c1","Fixed Probe Data",1200,600); 

   c1->cd();
   gFPAVG->Draw("ap"); 
   gm2fieldUtil::Graph::SetGraphLabels(gFPAVG,"Fixed Probe Avg","","Frequency (Hz)");
   gm2fieldUtil::Graph::UseTimeDisplay(gFPAVG); 
   gm2fieldUtil::Graph::SetGraphLabelSizes(gFPAVG,0.04,0.04);
   gFPAVG->Draw("ap");
   myErrorBand->Draw("same 4");  
   myFit->Draw("same");
   gSp->Draw("same"); 
   L->Draw("same"); 
   c1->Update();  

   std::cout << "FIT PARAMETERS" << std::endl;
   const int NPAR = nOrder+2;   
   double par[NPAR],parErr[NPAR];   
   for(int i=0;i<NPAR;i++){
      par[i]    = myFit->GetParameter(i); 
      parErr[i] = myFit->GetParError(i); 
      std::cout << Form("%.3E +/- %.3E",par[i],parErr[i]) << std::endl;
   }

   double mean=0,stdev=0;
   TGraph *gFPAVG_alt = RemoveTrend(gFPAVG,myFit,mean,stdev); 
   gm2fieldUtil::Graph::SetGraphParameters(gFPAVG_alt,20,kBlack);

   std::cout << Form("TREND SUBTRACTED: MEAN = %.3lf, STDEV = %.3lf",mean,stdev) << std::endl; 

   TCanvas *c2 = new TCanvas("c2","Fixed Probe Data [Trend Subtracted]",1200,600); 

   c2->cd();
   gFPAVG_alt->Draw("ap"); 
   gm2fieldUtil::Graph::SetGraphLabels(gFPAVG_alt,"Fixed Probe Avg","","Frequency (Hz)");
   gm2fieldUtil::Graph::UseTimeDisplay(gFPAVG_alt); 
   gm2fieldUtil::Graph::SetGraphLabelSizes(gFPAVG_alt,0.04,0.04);
   gFPAVG_alt->Draw("ap");
   c2->Update();  

   return 0;
}
