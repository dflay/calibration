// Compute azimuthal gradients in trolley data (after shimming) 

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
#include "TStyle.h"

#include "RootTreeStructs.h"
#include "gm2fieldMath.h"
#include "gm2fieldGraph.h"
#include "gm2fieldRootHelper.h"
#include "gm2fieldUnits.h"
#include "TemperatureSensor.h"

#include "./include/date.h"
#include "./include/Constants.h"
#include "./include/drift.h"
#include "./include/plungingProbeAnaEvent.h"
#include "./include/trolleyAnaEvent.h"

#include "./src/FXPRFuncs.C"
#include "./src/Consolidate.C"
#include "./src/CustomMath.C"
#include "./src/CustomGraph.C"
#include "./src/CustomImport.C"
#include "./src/CustomExport.C"
#include "./src/CustomAlgorithms.C"
#include "./src/CustomUtilities.C"
#include "./src/DeltaBFuncs.C"

int GetShimmedAziGrad(std::string date,std::string fitFuncStr,int probeNumber,bool isBlind){
   
   TString fitFunc = Form("%s",fitFuncStr.c_str());
   // TString fitFunc = Form("pol1");

   int rc=0;
   int method = gm2fieldUtil::Constants::kPhaseDerivative;

   date_t theDate;
   GetDate(theDate);

   char plotDir[200];
   sprintf(plotDir,"./plots/%02d-%02d-%02d",theDate.month,theDate.day,theDate.year-2000);
   MakeDirectory(plotDir); 

   blind_t blind;
   ImportBlinding(blind);
   double blindValue = blind.value_tr; 

   std::cout << "----------------------------------" << std::endl;
   std::cout << "SHIMMED AZI GRADIENT CALCULATION" << std::endl;

   char inpath[200];
   sprintf(inpath,"./input/runlists/%s/trly-shimmed_azi_%s.csv",date.c_str(),date.c_str());

   std::vector<int> run;
   std::vector<double> sf;
   std::vector<std::string> label;
   ImportDeltaBFileList_csv(inpath,run,label,sf);

   const int NRUN = run.size();
   if(NRUN==0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   // Trolley data
   std::vector<trolleyAnaEvent_t> trlyData; // for drift  
   std::vector<trolleyAnaEvent_t> Data,DataCor,DataCorAlt;    

   for(int i=0;i<NRUN;i++) rc = GetTrolleyData(date,run[i],method,Data);
   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   if(isBlind) rc = ApplyBlindingTRLY(blindValue,Data); 
   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   // fixed probe data
   std::string fxprPath = "./input/probe-lists/fxpr-list.csv";
   std::vector<int> fxprList;
   gm2fieldUtil::Import::ImportData1<int>(fxprPath,"csv",fxprList);

   std::vector<gm2field::fixedProbeFrequency_t> fxprData;
   for(int i=0;i<NRUN;i++) rc = gm2fieldUtil::RootHelper::GetFPFrequencies(run[i],fxprData);
   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   // trolley data for drift correction
   std::vector<int> inList,trlyList; 
   std::string trlyPath = "./input/probe-lists/trly-list.csv"; 
   gm2fieldUtil::Import::ImportData1<int>(trlyPath,"csv",inList);  
   for(int i=0;i<NRUN;i++) rc = GetTrolleyData(date,run[i],method,trlyData); 
   if(rc!=0) return 1;

   // remove the trolley probe of interest 
   int NL = inList.size();
   for(int i=0;i<NL;i++) if(probeNumber!=inList[i]) trlyList.push_back(inList[i]); 

   std::cout << "Applying field drift corrections (during measurement)..." << std::endl;
   // correct for field drift 
   // use fxpr 
   rc = CorrectTRLYForDriftDuringMeasurement(method,fxprList,fxprData,Data,DataCor);
   if(rc!=0) return 1;
   // use trly 
   rc = CorrectTRLYForDriftDuringMeasurement(method,trlyList,trlyData,Data,DataCorAlt);
   if(rc!=0) return 1;
   std::cout << "--> Done" << std::endl; 

   // find gradients 
   // azimuthal 
   TGraph *gAB         = GetTRLYTGraph(probeNumber,"phi","freq",Data); 
   TGraph *gAB_fxpr    = GetTRLYTGraph(probeNumber,"phi","freq",DataCor); 
   TGraph *gAB_trly    = GetTRLYTGraph(probeNumber,"phi","freq",DataCorAlt); 

   gm2fieldUtil::Graph::SetGraphParameters(gAB     ,20,kBlack); 
   gm2fieldUtil::Graph::SetGraphParameters(gAB_fxpr,20,kBlack); 
   gm2fieldUtil::Graph::SetGraphParameters(gAB_trly,20,kBlack); 
   
   TString xAxisTitle = Form("Azimuthal Position (deg)");
   TString yAxisTitle = Form("Field Strength (Hz)");

   double xmin = 189.15; 
   double xmax = 189.35;  

   TCanvas *c1 = new TCanvas("c1","Azimuthal",1200,800);
   c1->Divide(1,3);

   gStyle->SetOptFit(0); 

   c1->cd(1); 
   gAB->Draw("ap");
   gm2fieldUtil::Graph::SetGraphLabels(gAB,"No Drift Cor",xAxisTitle,yAxisTitle);
   gm2fieldUtil::Graph::SetGraphLabelSizes(gAB,0.06,0.06); 
   gAB->Draw("ap");
   gAB->Fit(fitFunc,"QR","",xmin,xmax); 
   c1->Update();

   c1->cd(2); 
   gAB_fxpr->Draw("ap");
   gm2fieldUtil::Graph::SetGraphLabels(gAB_fxpr,"Drift Cor (FXPR)",xAxisTitle,yAxisTitle);
   gm2fieldUtil::Graph::SetGraphLabelSizes(gAB_fxpr,0.06,0.06); 
   gAB_fxpr->Draw("ap");
   gAB_fxpr->Fit(fitFunc,"QR","",xmin,xmax); 
   c1->Update();

   c1->cd(3); 
   gAB_trly->Draw("ap");
   gm2fieldUtil::Graph::SetGraphLabels(gAB_trly,"Drift Cor (TRLY)",xAxisTitle,yAxisTitle);
   gm2fieldUtil::Graph::SetGraphLabelSizes(gAB_trly,0.06,0.06); 
   gAB_trly->Draw("ap");
   gAB_trly->Fit(fitFunc,"QR","",xmin,xmax); 
   c1->Update();

   TString plotPath = Form("%s/shimmed-azi-grad_%s.png",plotDir,date.c_str());
   c1->cd();
   c1->Print(plotPath); 

   double PA[3],EA[3];
   TF1 *fitAB      = gAB->GetFunction(fitFunc);
   TF1 *fitAB_fxpr = gAB_fxpr->GetFunction(fitFunc);
   TF1 *fitAB_trly = gAB_trly->GetFunction(fitFunc);

   PA[0] = fitAB->GetParameter(1);
   PA[1] = fitAB_fxpr->GetParameter(1);
   PA[2] = fitAB_trly->GetParameter(1);

   EA[0] = fitAB->GetParError(1);
   EA[1] = fitAB_fxpr->GetParError(1);
   EA[2] = fitAB_trly->GetParError(1);

   // convert to Hz/mm
   // originally in deg, where 1 deg = 124 mm 
   for(int i=0;i<3;i++){
      PA[i] /= 124.;  
      EA[i] /= 124.;  
   }
  
   double drift[3]     = {0,0,0}; 
   double drift_err[3] = {0,0,0}; 

   std::cout << "============================ RESULTS ============================" << std::endl; 
   std::cout << "Azimuthal:        " << Form("%.3lf +/- %.3lf Hz/mm",PA[0],EA[0]) << std::endl;
   std::cout << "Azimuthal [fxpr]: " << Form("%.3lf +/- %.3lf Hz/mm",PA[1],EA[1]) << std::endl;
   std::cout << "Azimuthal [trly]: " << Form("%.3lf +/- %.3lf Hz/mm",PA[2],EA[2]) << std::endl;

   char outpath[200],outdir[200];
   if(isBlind)  sprintf(outdir,"./output/blinded/%02d-%02d-%02d"  ,theDate.month,theDate.day,theDate.year-2000); 
   if(!isBlind) sprintf(outdir,"./output/unblinded/%02d-%02d-%02d",theDate.month,theDate.day,theDate.year-2000); 
   rc = MakeDirectory(outdir);
   sprintf(outpath,"%s/azi-grad_final-location_pr-%02d_%s.csv"  ,outdir,probeNumber,date.c_str()); 
   PrintToFile(outpath,"azi-grad",PA,EA,drift,drift_err);

   // Fixed probe plot 
   const int NTR = Data.size();
   TGraph *gFXPR = GetTGraph(method,Data[0].time[0],Data[NTR-1].time[3],1E+9,fxprList,fxprData);
   gm2fieldUtil::Graph::SetGraphParameters(gFXPR,20,kBlack); 
 
   TCanvas *c2 = new TCanvas("c2","Average FXPR Data",1200,600);
   
   c2->cd();
   gFXPR->Draw("alp");
   gm2fieldUtil::Graph::SetGraphLabels(gFXPR,"Fixed Probe Average","","Frequency (Hz)");
   gm2fieldUtil::Graph::UseTimeDisplay(gFXPR);
   // gm2fieldUtil::Graph::SetGraphLabelSizes(gFXPR,0.05,0.06);
   gFXPR->Draw("alp");
   c2->Update();
 
   plotPath = Form("%s/shimmed-azi-grad_avg-fxpr_pr-%02d_%s.png",plotDir,probeNumber,date.c_str()); 
   c2->Print(plotPath);   

   return 0;
}
