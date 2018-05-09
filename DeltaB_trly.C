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

#include "./include/date.h"
#include "./include/Constants.h"
#include "./include/drift.h"
#include "./include/plungingProbeAnaEvent.h"
#include "./include/trolleyAnaEvent.h"

#include "./src/InputManager.C"
#include "./src/FitFuncs.C"
#include "./src/FXPRFuncs.C"
#include "./src/Consolidate.C"
#include "./src/CustomMath.C"
#include "./src/CustomGraph.C"
#include "./src/CustomImport.C"
#include "./src/CustomExport.C"
#include "./src/CustomAlgorithms.C"
#include "./src/CustomUtilities.C"
#include "./src/DeltaBFuncs.C"

int DeltaB_trly(std::string configFile){

   std::cout << "----------------------------------" << std::endl;
   std::cout << "TRLY DELTA B CALCULATION" << std::endl;

   int rc=0;
   int method = gm2fieldUtil::Constants::kPhaseDerivative;

   InputManager *inputMgr = new InputManager();
   inputMgr->Load(configFile);
   inputMgr->Print();

   std::string date    = inputMgr->GetAnalysisDate();
   bool isBlind        = inputMgr->IsBlind();
   bool finalPos       = inputMgr->IsFinalLocation();
   bool isFullAnalysis = inputMgr->IsFullAnalysis();
   bool useP2PFit      = inputMgr->UseP2PFit();
   int probeNumber     = inputMgr->GetTrolleyProbe();
   int axis            = inputMgr->GetAxis();

   date_t theDate; 
   GetDate(theDate); 

   blind_t blind;
   ImportBlinding(blind);
   double blindValue = blind.value_tr;

   std::string gradName;
   if(axis==0) gradName = "rad"; 
   if(axis==1) gradName = "vert"; 
   if(axis==2) gradName = "azi"; 

   char tag[100];
   if(finalPos==1){
      sprintf(tag,"_final-location_");
   }else if(finalPos==0){
      sprintf(tag,"_");
   }else{
      std::cout << "Invalid type! finalPos = " << finalPos <<std::endl;
      return 1;
   }

   // char inpath[200];
   // sprintf(inpath,"./input/runlists/%s/dB-trly%s%s-grad_%s.csv",date.c_str(),tag,gradName.c_str(),date.c_str());
   // std::vector<int> allRuns;
   // std::vector<double> sf;
   // std::vector<std::string> label;
   // ImportDeltaBFileList_csv(inpath,allRuns,label,sf);

   std::vector<int> allRuns;
   std::vector<std::string> label;
   inputMgr->GetRunList(allRuns);
   inputMgr->GetRunLabels(label);

   const int NRUN = allRuns.size();
   if(NRUN==0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   std::vector<int> index,run,driftRun;
   SortRuns(label,allRuns,run,driftRun,index);

   int bareIndex  = index[0];
   int gradIndex  = index[1];
   int bare2Index = index[2];

   // fixed probe data
   std::string fxprPath = "./input/probe-lists/fxpr-list.csv"; 
   std::vector<int> fxprList; 
   gm2fieldUtil::Import::ImportData1<int>(fxprPath,"csv",fxprList);

   // now get a TGraph we can use in interpolating across these runs
   std::vector<double> stats; 
   TGraph *gDrift = GetDriftTGraph(method,driftRun,fxprList,stats);

   // get data for drift correction *during* a run 
   // all data (not using P2P fitting) 
   const int N3 = run.size();
   std::vector<gm2field::fixedProbeFrequency_t> fxprData;
   for(int i=0;i<N3;i++) rc = gm2fieldUtil::RootHelper::GetFPFrequencies(run[i],fxprData);
   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   // using P2P fitting 
   // bare run 
   std::vector<gm2field::fixedProbeFrequency_t> fxprData_bare;
   rc = gm2fieldUtil::RootHelper::GetFPFrequencies(run[bareIndex],fxprData_bare);
   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }
   // grad run
   std::vector<gm2field::fixedProbeFrequency_t> fxprData_grad;
   rc = gm2fieldUtil::RootHelper::GetFPFrequencies(run[gradIndex],fxprData_grad);
   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   int NN  = fxprData_bare.size();
   int NFP = fxprList.size();
   int startIndex = fxprList[0];
   int stopIndex  = fxprList[NFP-1];
   unsigned long long t0     = 0;
   unsigned long long tStart = fxprData_bare[0].GpsTimeStamp[startIndex];
   unsigned long long tStop  = fxprData_bare[NN-1].GpsTimeStamp[stopIndex];
   unsigned long long tStep  = 1E+9;

   // for graphing 
   double xMin_bare = tStart/1E+9;
   double xMax_bare = tStop/1E+9;

   std::vector<fixedProbeEvent_t> fxprDataAvg_bare;
   GetAverageFXPRVectorsNew(method,t0,tStart,tStop,tStep,fxprList,fxprData_bare,fxprDataAvg_bare);

   NN = fxprData_grad.size();
   tStart = fxprData_grad[0].GpsTimeStamp[startIndex];
   tStop  = fxprData_grad[NN-1].GpsTimeStamp[stopIndex];
   std::vector<fixedProbeEvent_t> fxprDataAvg_grad;
   GetAverageFXPRVectorsNew(method,t0,tStart,tStop,tStep,fxprList,fxprData_grad,fxprDataAvg_grad);

   // for graphing 
   double xMin_grad = tStart/1E+9;
   double xMax_grad = tStop/1E+9;

   std::cout << "--> Done." << std::endl;

   // Fixed probe average plot 
   TGraph *gFPAVG_bare = GetTGraphNew(fxprDataAvg_bare);
   gm2fieldUtil::Graph::SetGraphParameters(gFPAVG_bare,20,kBlack);

   TGraph *gFPAVG_grad = GetTGraphNew(fxprDataAvg_grad);
   gm2fieldUtil::Graph::SetGraphParameters(gFPAVG_grad,20,kBlack);

   // get a fit function to these FXPR data
   const int nOrder = 4;
   std::vector< std::vector<double> > eps_bare,eps_grad;
   TF1 *fxprFit_bare = GetPolyFitToFXPR("myFitBare",nOrder,gFPAVG_bare,eps_bare);
   TF1 *fxprFit_grad = GetPolyFitToFXPR("myFitGrad",nOrder,gFPAVG_grad,eps_grad);

   double mean_p2p=0;
   double stdev_p2p_bare=0,stdev_p2p_grad=0;
   double STDEV_P2P_BARE=0,STDEV_P2P_GRAD=0;
   TGraph *gFPAVG_bare_trm = RemoveTrend(gFPAVG_bare,fxprFit_bare,mean_p2p,STDEV_P2P_BARE);
   TGraph *gFPAVG_grad_trm = RemoveTrend(gFPAVG_grad,fxprFit_grad,mean_p2p,STDEV_P2P_GRAD);

   if(useP2PFit){
      stdev_p2p_bare = STDEV_P2P_BARE;
      stdev_p2p_grad = STDEV_P2P_GRAD;
   }

   // Trolley data 
   std::vector<trolleyAnaEvent_t> trlyData;   // for drift correction  
   std::vector<trolleyAnaEvent_t> EventBare,EventBareCor,EventBareCorAlt; 
   std::vector<trolleyAnaEvent_t> EventGrad,EventGradCor,EventGradCorAlt;

   std::cout << "Getting run " << run[bareIndex] << std::endl; 
   rc = GetTrolleyData(date,run[bareIndex],method,EventBare);
   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   std::cout << "Getting run " << run[gradIndex] << std::endl; 
   rc = GetTrolleyData(date,run[gradIndex],method,EventGrad);
   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   if(isBlind) ApplyBlindingTRLY(blindValue,EventBare); 
   if(isBlind) ApplyBlindingTRLY(blindValue,EventGrad); 

   // trolley data for drift correction
   std::vector<int> inList,trlyList; 
   std::string trlyPath = "./input/probe-lists/trly-list.csv"; 
   gm2fieldUtil::Import::ImportData1<int>(trlyPath,"csv",inList);  
   for(int i=0;i<N3;i++) rc = GetTrolleyData(date,run[i],method,trlyData); 
   if(rc!=0) return 1;

   // remove the trolley probe of interest 
   int NL = inList.size();
   for(int i=0;i<NL;i++) if(probeNumber!=inList[i]) trlyList.push_back(inList[i]); 

   std::cout << "Applying field drift corrections (during measurement)..." << std::endl;
   // correct for field drift 
   // use fxpr 
   if(useP2PFit){
      rc = CorrectTRLYForDriftDuringMeasurement(method,fxprFit_bare,EventBare,EventBareCor);
      rc = CorrectTRLYForDriftDuringMeasurement(method,fxprFit_grad,EventGrad,EventGradCor);
   }else{
      // rc = CorrectTRLYForDriftDuringMeasurement(method,fxprList,fxprData,EventBare,EventBareCor);
      // rc = CorrectTRLYForDriftDuringMeasurement(method,fxprList,fxprData,EventGrad,EventGradCor);
      rc = CorrectTRLYForDriftDuringMeasurement(fxprDataAvg_bare,EventBare,EventBareCor);
      rc = CorrectTRLYForDriftDuringMeasurement(fxprDataAvg_grad,EventGrad,EventGradCor);
   }
   if(rc!=0) return 1;

   // use trly 
   rc = CorrectTRLYForDriftDuringMeasurement(method,trlyList,trlyData,EventBare,EventBareCorAlt);
   rc = CorrectTRLYForDriftDuringMeasurement(method,trlyList,trlyData,EventGrad,EventGradCorAlt);
   if(rc!=0) return 1;
   std::cout << "--> Done" << std::endl; 

   // now do delta B calcs
   double dB[3]        = {0,0,0}; 
   double dB_err[3]    = {0,0,0}; 
   double drift[3]     = {0,0,0};  
   double drift_err[3] = {0,0,0};  

   // double dB_grad=0,dB_grad_err=0;
   CalculateTRLYDeltaB_Stationary(false,method,probeNumber,gDrift,EventBare,EventGrad,dB[0],dB_err[0],drift[0],drift_err[0]);
   // drift corrected (fxpr)  
   CalculateTRLYDeltaB_Stationary(true,method,probeNumber,gDrift,EventBareCor,EventGradCor,dB[1],dB_err[1],drift[1],drift_err[1]);
   // drift corrected (trly)  
   CalculateTRLYDeltaB_Stationary(true,method,probeNumber,trlyList,trlyData,EventBareCorAlt,EventGradCorAlt,dB[2],dB_err[2],drift[2],drift_err[2]);

   // add the stdev from the anchor points as the drift error
   // also add in stdev of P2P stdev after subtracting off the fit for both bare and grad runs
   double arg_sq=0;
   for(int i=1;i<3;i++){
      arg_sq       = stats[1]*stats[1] + stats[3]*stats[3] + stdev_p2p_bare*stdev_p2p_bare + stdev_p2p_grad*stdev_p2p_grad;
      drift_err[i] = TMath::Sqrt(arg_sq);
   }

   std::cout << Form("===================== RESULTS FOR PROBE %02d =====================",probeNumber+1) << std::endl;
   std::cout << "Raw results: " << std::endl; 
   std::cout << Form("dB (%s): %.3lf +/- %.3lf Hz (%.3lf +/- %.3lf ppb)",label[gradIndex].c_str(),
                     dB[0],dB_err[0],dB[0]/0.06179,dB_err[0]/0.06179) << std::endl;
   std::cout << "-----------------------------------------------------------------------" << std::endl;  
   std::cout << "Drift corrected [FXPR]: " << std::endl;
   std::cout << Form("dB (%s): %.3lf +/- %.3lf +/- %.3lf Hz (%.3lf +/- %.3lf +/- %.3lf ppb)",label[gradIndex].c_str(),
                     dB[1],dB_err[1],drift_err[1],dB[1]/0.06179,dB_err[1]/0.06179,drift_err[1]/0.06179) << std::endl;
   std::cout << "-----------------------------------------------------------------------" << std::endl;  
   std::cout << "Drift corrected [TRLY]: " << std::endl;
   std::cout << Form("dB (%s): %.3lf +/- %.3lf +/- %.3lf Hz (%.3lf +/- %.3lf +/- %.3lf ppb)",label[gradIndex].c_str(),
                     dB[2],dB_err[2],drift_err[2],dB[2]/0.06179,dB_err[2]/0.06179,drift_err[2]/0.06179) << std::endl;

   char outpath[200],outdir[200];
   if(isBlind)  sprintf(outdir,"./output/blinded/%02d-%02d-%02d"  ,theDate.month,theDate.day,theDate.year-2000); 
   if(!isBlind) sprintf(outdir,"./output/unblinded/%02d-%02d-%02d",theDate.month,theDate.day,theDate.year-2000); 
   if(!isFullAnalysis){
      // we're doing a delta-b calculation during a measurement
      // store in a different location 
      sprintf(outdir,"./output/delta-b/%02d-%02d-%02d",theDate.month,theDate.day,theDate.year-2000);
   }

   rc = MakeDirectory(outdir);
   sprintf(outpath,"%s/dB-trly%s%s-grad_pr-%02d_%s.csv"  ,outdir,tag,gradName.c_str(),probeNumber,date.c_str()); 
   rc = PrintToFile(outpath,gradName.c_str(),dB,dB_err,drift,drift_err); 

   // make some plots 
   // Trolley plots 
   TGraph *gTR_bare = GetTRLYTGraph(probeNumber,"GpsTimeStamp","freq",EventBare);
   gm2fieldUtil::Graph::SetGraphParameters(gTR_bare,20,kBlack);

   TGraph *gTR_bare_fxpr = GetTRLYTGraph(probeNumber,"GpsTimeStamp","freq",EventBareCor);
   gm2fieldUtil::Graph::SetGraphParameters(gTR_bare_fxpr,20,kRed);

   TGraph *gTR_grad = GetTRLYTGraph(probeNumber,"GpsTimeStamp","freq",EventGrad);
   gm2fieldUtil::Graph::SetGraphParameters(gTR_grad,20,kBlack);

   TGraph *gTR_grad_fxpr = GetTRLYTGraph(probeNumber,"GpsTimeStamp","freq",EventGradCor);
   gm2fieldUtil::Graph::SetGraphParameters(gTR_grad_fxpr,20,kRed);

   TMultiGraph *mgTR_bare = new TMultiGraph();
   mgTR_bare->Add(gTR_bare     ,"lp");
   mgTR_bare->Add(gTR_bare_fxpr,"lp");  

   TMultiGraph *mgTR_grad = new TMultiGraph();
   mgTR_grad->Add(gTR_grad     ,"lp");
   mgTR_grad->Add(gTR_grad_fxpr,"lp");  

   TLegend *LT = new TLegend(0.6,0.6,0.8,0.8);
   LT->AddEntry(gTR_bare     ,"Raw","p");
   LT->AddEntry(gTR_bare_fxpr,"FXPR Drift Corrected","p");

   char plotDir[200];
   sprintf(plotDir,"./plots/%02d-%02d-%02d",theDate.month,theDate.day,theDate.year-2000);
   rc = MakeDirectory(plotDir);

   TCanvas *c1 = new TCanvas("c1","TRLY & FXPR Data",1200,600);
   c1->Divide(1,2);

   c1->cd(1);
   mgTR_bare->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mgTR_bare,Form("TRLY Run %d",run[bareIndex]),"","Frequency (Hz)");
   gm2fieldUtil::Graph::UseTimeDisplay(mgTR_bare);
   gm2fieldUtil::Graph::SetGraphLabelSizes(mgTR_bare,0.05,0.06);
   mgTR_bare->GetXaxis()->SetLimits(xMin_bare,xMax_bare); 
   mgTR_bare->Draw("a");
   LT->Draw("same");
   c1->Update();

   c1->cd(2);
   gFPAVG_bare->Draw("alp");
   gm2fieldUtil::Graph::SetGraphLabels(gFPAVG_bare,"Fixed Probe Average","","Frequency (Hz)");
   gm2fieldUtil::Graph::UseTimeDisplay(gFPAVG_bare);
   gm2fieldUtil::Graph::SetGraphLabelSizes(gFPAVG_bare,0.05,0.06);
   gFPAVG_bare->GetXaxis()->SetLimits(xMin_bare,xMax_bare); 
   gFPAVG_bare->Draw("alp");
   if(useP2PFit) fxprFit_bare->Draw("same");
   c1->Update();

   TString plotPath = Form("%s/trly_pr-%02d_bare-run-%d.png",plotDir,probeNumber,run[bareIndex]);
   c1->cd();
   c1->Print(plotPath);

   TCanvas *c2 = new TCanvas("c2","TRLY & FXPR Data",1200,600);
   c2->Divide(1,2);

   c2->cd(1);
   mgTR_grad->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mgTR_grad,Form("TRLY Run %d",run[gradIndex]),"","Frequency (Hz)");
   gm2fieldUtil::Graph::UseTimeDisplay(mgTR_grad);
   gm2fieldUtil::Graph::SetGraphLabelSizes(mgTR_grad,0.05,0.06);
   mgTR_grad->GetXaxis()->SetLimits(xMin_grad,xMax_grad); 
   mgTR_grad->Draw("a");
   LT->Draw("same");
   c2->Update();

   c2->cd(2);
   gFPAVG_grad->Draw("alp");
   gm2fieldUtil::Graph::SetGraphLabels(gFPAVG_grad,"Fixed Probe Average","","Frequency (Hz)");
   gm2fieldUtil::Graph::UseTimeDisplay(gFPAVG_grad);
   gm2fieldUtil::Graph::SetGraphLabelSizes(gFPAVG_grad,0.05,0.06);
   gFPAVG_grad->GetXaxis()->SetLimits(xMin_grad,xMax_grad); 
   gFPAVG_grad->Draw("alp");
   if(useP2PFit) fxprFit_grad->Draw("same");
   c2->Update();

   plotPath = Form("%s/trly_pr-%02d_grad-run-%d.png",plotDir,probeNumber,run[gradIndex]);
   c2->cd();
   c2->Print(plotPath);

   return 0;
}

