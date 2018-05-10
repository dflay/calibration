// Compute final transverse gradients in trolley data (after shimming) 

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

#include "./src/InputManager.C"
#include "./src/FXPRFuncs.C"
#include "./src/Consolidate.C"
#include "./src/CustomMath.C"
#include "./src/CustomGraph.C"
#include "./src/CustomImport.C"
#include "./src/CustomExport.C"
#include "./src/CustomAlgorithms.C"
#include "./src/CustomUtilities.C"
#include "./src/DeltaBFuncs.C"

int GetShimmedTransGrad(std::string configFile){

   std::cout << "----------------------------------" << std::endl;
   std::cout << "SHIMMED TRANSVERSE GRADIENT CALCULATION" << std::endl;

   int rc=0;
   int method = gm2fieldUtil::Constants::kPhaseDerivative;

   InputManager *inputMgr = new InputManager();
   inputMgr->Load(configFile);
   inputMgr->Print();

   std::string date    = inputMgr->GetAnalysisDate();
   std::string fitFunc = inputMgr->GetFitFunction();
   bool isBlind        = inputMgr->IsBlind();
   int probeNumber     = inputMgr->GetTrolleyProbe();

   date_t theDate; 
   GetDate(theDate);

   char plotDir[200];
   sprintf(plotDir,"./plots/%02d-%02d-%02d",theDate.month,theDate.day,theDate.year-2000);
   MakeDirectory(plotDir); 

   blind_t blind;
   ImportBlinding(blind);
   double blindValue = blind.value_tr; 

   std::vector<int> run;
   std::vector<std::string> label;
   inputMgr->GetRunList(run);
   inputMgr->GetRunLabels(label);

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

   const int NN = fxprData.size();
   const int NFP = fxprList.size();
   unsigned long long t0     = 0;
   unsigned long long tStart = fxprData[0].GpsTimeStamp[fxprList[0]];
   unsigned long long tStop  = fxprData[NN-1].GpsTimeStamp[fxprList[NFP-1]];
   unsigned long long tStep  = 1E+9;

   std::vector<fixedProbeEvent_t> fxprDataAvg;
   GetAverageFXPRVectorsNew(method,t0,tStart,tStop,tStep,fxprList,fxprData,fxprDataAvg);

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
   rc = CorrectTRLYForDriftDuringMeasurement(fxprDataAvg,Data,DataCor);
   if(rc!=0) return 1;
   // use trly 
   rc = CorrectTRLYForDriftDuringMeasurement(method,trlyList,trlyData,Data,DataCorAlt);
   if(rc!=0) return 1;
   std::cout << "--> Done" << std::endl; 

   // find gradients 

   // radial
   TGraph *gRB         = GetSlicePlot('r',Data); 
   TGraph *gRB_fxpr    = GetSlicePlot('r',DataCor); 
   TGraph *gRB_trly    = GetSlicePlot('r',DataCorAlt); 

   // vertical  
   TGraph *gVB         = GetSlicePlot('v',Data); 
   TGraph *gVB_fxpr    = GetSlicePlot('v',DataCor); 
   TGraph *gVB_trly    = GetSlicePlot('v',DataCorAlt);

   gm2fieldUtil::Graph::SetGraphParameters(gRB     ,20,kBlack); 
   gm2fieldUtil::Graph::SetGraphParameters(gRB_fxpr,20,kBlack); 
   gm2fieldUtil::Graph::SetGraphParameters(gRB_trly,20,kBlack); 

   gm2fieldUtil::Graph::SetGraphParameters(gVB     ,20,kBlack); 
   gm2fieldUtil::Graph::SetGraphParameters(gVB_fxpr,20,kBlack); 
   gm2fieldUtil::Graph::SetGraphParameters(gVB_trly,20,kBlack); 
   
   TString xAxisTitle = Form("r (cm)");
   TString yAxisTitle = Form("Field Strength (Hz)");
 
   TCanvas *c1 = new TCanvas("c1","Radial",1200,800);
   c1->Divide(1,3);
   
   gStyle->SetOptFit(0); 

   c1->cd(1); 
   gRB->Draw("ap");
   gm2fieldUtil::Graph::SetGraphLabels(gRB,"No Drift Cor",xAxisTitle,yAxisTitle);
   gm2fieldUtil::Graph::SetGraphLabelSizes(gRB,0.05,0.06); 
   gRB->GetXaxis()->SetLimits(-4.5,4.5); 
   gRB->Draw("ap");
   gRB->Fit(fitFunc.c_str(),"Q");
   c1->Update();

   c1->cd(2); 
   gRB_fxpr->Draw("ap");
   gm2fieldUtil::Graph::SetGraphLabels(gRB_fxpr,"Drift Cor (FXPR)",xAxisTitle,yAxisTitle);
   gm2fieldUtil::Graph::SetGraphLabelSizes(gRB_fxpr,0.05,0.06); 
   gRB_fxpr->GetXaxis()->SetLimits(-4.5,4.5); 
   gRB_fxpr->Draw("ap");
   gRB_fxpr->Fit(fitFunc.c_str(),"Q");
   c1->Update();

   c1->cd(3); 
   gRB_trly->Draw("ap");
   gm2fieldUtil::Graph::SetGraphLabels(gRB_trly,"Drift Cor (TRLY)",xAxisTitle,yAxisTitle);
   gm2fieldUtil::Graph::SetGraphLabelSizes(gRB_trly,0.05,0.06); 
   gRB_trly->GetXaxis()->SetLimits(-4.5,4.5); 
   gRB_trly->Draw("ap");
   gRB_trly->Fit(fitFunc.c_str(),"Q");
   c1->Update();

   TString plotPath = Form("%s/trly-shimmed_rad-grad_run-%d_%s.png",plotDir,run[0],date.c_str());
   c1->cd(); 
   c1->Print(plotPath); 

   xAxisTitle = Form("y (cm)");

   TCanvas *c2 = new TCanvas("c2","Vertical",1200,800);
   c2->Divide(1,3);

   gStyle->SetOptFit(0); 

   c2->cd(1); 
   gVB->Draw("ap");
   gm2fieldUtil::Graph::SetGraphLabels(gVB,"No Drift Cor",xAxisTitle,yAxisTitle);
   gm2fieldUtil::Graph::SetGraphLabelSizes(gVB,0.05,0.06); 
   gVB->GetXaxis()->SetLimits(-4.5,4.5); 
   gVB->Draw("ap");
   gVB->Fit(fitFunc.c_str(),"Q"); 
   c2->Update();

   c2->cd(2); 
   gVB_fxpr->Draw("ap");
   gm2fieldUtil::Graph::SetGraphLabels(gVB_fxpr,"Drift Cor (FXPR)",xAxisTitle,yAxisTitle);
   gm2fieldUtil::Graph::SetGraphLabelSizes(gVB_fxpr,0.05,0.06); 
   gVB_fxpr->GetXaxis()->SetLimits(-4.5,4.5); 
   gVB_fxpr->Draw("ap");
   gVB_fxpr->Fit(fitFunc.c_str(),"Q"); 
   c2->Update();

   c2->cd(3); 
   gVB_trly->Draw("ap");
   gm2fieldUtil::Graph::SetGraphLabels(gVB_trly,"Drift Cor (TRLY)",xAxisTitle,yAxisTitle);
   gm2fieldUtil::Graph::SetGraphLabelSizes(gVB_trly,0.05,0.06); 
   gVB_trly->GetXaxis()->SetLimits(-4.5,4.5); 
   gVB_trly->Draw("ap");
   gVB_trly->Fit(fitFunc.c_str(),"Q"); 
   c2->Update();

   plotPath = Form("%s/trly-shimmed_vert-grad_run-%d_%s.png",plotDir,run[0],date.c_str());
   c2->cd(); 
   c2->Print(plotPath); 

   double PR[3],ER[3];
   TF1 *fitRB      = gRB->GetFunction( fitFunc.c_str() );
   TF1 *fitRB_fxpr = gRB_fxpr->GetFunction( fitFunc.c_str() );
   TF1 *fitRB_trly = gRB_trly->GetFunction( fitFunc.c_str() );

   double PV[3],EV[3];
   TF1 *fitVB      = gVB->GetFunction( fitFunc.c_str() );
   TF1 *fitVB_fxpr = gVB_fxpr->GetFunction( fitFunc.c_str() );
   TF1 *fitVB_trly = gVB_trly->GetFunction( fitFunc.c_str() );

   PR[0] = fitRB->GetParameter(1);
   PR[1] = fitRB_fxpr->GetParameter(1);
   PR[2] = fitRB_trly->GetParameter(1);

   ER[0] = fitRB->GetParError(1);
   ER[1] = fitRB_fxpr->GetParError(1);
   ER[2] = fitRB_trly->GetParError(1);

   PV[0] = fitVB->GetParameter(1);
   PV[1] = fitVB_fxpr->GetParameter(1);
   PV[2] = fitVB_trly->GetParameter(1);

   EV[0] = fitVB->GetParError(1);
   EV[1] = fitVB_fxpr->GetParError(1);
   EV[2] = fitVB_trly->GetParError(1);

   // convert to Hz/mm 
   for(int i=0;i<3;i++){
      PR[i] /= 10.;
      ER[i] /= 10.;
      PV[i] /= 10.;
      EV[i] /= 10.;
   }

   std::cout << "============================ RESULTS ============================" << std::endl; 
   std::cout << "Radial:          " << Form("%.3lf +/- %.3lf Hz/mm",PR[0],ER[0]) << std::endl;
   std::cout << "Radial [fxpr]:   " << Form("%.3lf +/- %.3lf Hz/mm",PR[1],ER[1]) << std::endl;
   std::cout << "Radial [trly]:   " << Form("%.3lf +/- %.3lf Hz/mm",PR[2],ER[2]) << std::endl;
   std::cout << "-----------------------------------------------------------------" << std::endl;
   std::cout << "Vertical:        " << Form("%.3lf +/- %.3lf Hz/mm",PV[0],EV[0]) << std::endl;
   std::cout << "Vertical [fxpr]: " << Form("%.3lf +/- %.3lf Hz/mm",PV[1],EV[1]) << std::endl;
   std::cout << "Vertical [trly]: " << Form("%.3lf +/- %.3lf Hz/mm",PV[2],EV[2]) << std::endl;

   // WARNING: since we do only in-run drift corrections, we're putting a zero error on the drift 
   //          between runs, since it doesn't exist here 
   double drift[3]     = {0,0,0}; 
   double drift_err[3] = {0,0,0}; 

   char outpath[200],outdir[200];
   if(isBlind)  sprintf(outdir,"./output/blinded/%02d-%02d-%02d"  ,theDate.month,theDate.day,theDate.year-2000); 
   if(!isBlind) sprintf(outdir,"./output/unblinded/%02d-%02d-%02d",theDate.month,theDate.day,theDate.year-2000); 
   rc = MakeDirectory(outdir);
   sprintf(outpath,"%s/rad-grad_final-location_pr-%02d_%s.csv"    ,outdir,probeNumber,date.c_str()); 
   PrintToFile(outpath,"rad-grad",PR,ER,drift,drift_err); 
   sprintf(outpath,"%s/vert-grad_final-location_pr-%02d_%s.csv"   ,outdir,probeNumber,date.c_str()); 
   PrintToFile(outpath,"vert-grad",PV,EV,drift,drift_err); 

   delete inputMgr; 

   return 0;
}

