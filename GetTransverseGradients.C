// Compute imposed transverse gradients using trolley data   

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
#include "./src/FXPRFuncs.C"
#include "./src/Consolidate.C"
#include "./src/CustomMath.C"
#include "./src/CustomGraph.C"
#include "./src/CustomImport.C"
#include "./src/CustomExport.C"
#include "./src/CustomAlgorithms.C"
#include "./src/CustomUtilities.C"
#include "./src/DeltaBFuncs.C"

int GetTransverseGradients(std::string configFile){

   std::cout << "----------------------------------" << std::endl;
   std::cout << "SCC TRANSVERSE GRADIENT CALCULATION" << std::endl;

   int rc=0;
   int method = gm2fieldUtil::Constants::kPhaseDerivative;

   InputManager *inputMgr = new InputManager();
   inputMgr->Load(configFile);
   inputMgr->Print();

   std::string date    = inputMgr->GetAnalysisDate();
   std::string fitFunc = inputMgr->GetFitFunction();
   bool isFullAnalysis = inputMgr->IsFullAnalysis();
   bool isBlind        = inputMgr->IsBlind();
   int probeNumber     = inputMgr->GetTrolleyProbe();

   date_t theDate;
   GetDate(theDate);

   char plotDir[200];
   sprintf(plotDir,"./plots/%02d-%02d-%02d",theDate.month,theDate.day,theDate.year-2000);
   if(!isFullAnalysis){
      // we're doing a delta-b calculation during a measurement
      // store in a different location 
      sprintf(plotDir,"./plots/delta-b/%02d-%02d-%02d",theDate.month,theDate.day,theDate.year-2000);
   }

   MakeDirectory(plotDir); 

   blind_t blind;
   ImportBlinding(blind);
   double blindValue = blind.value_tr;

   // char inpath[200];
   // sprintf(inpath,"./input/runlists/%s/trans-grad_%s.csv",date.c_str(),date.c_str());
   // std::vector<int> allRuns,run,vRun,rRun,aRun,bRun,b2Run;
   // std::vector<double> sf;
   // std::vector<std::string> label;
   // ImportDeltaBFileList_csv(inpath,allRuns,label,sf);

   std::vector<int> allRuns,run,vRun,rRun,aRun,bRun,b2Run;
   std::vector<std::string> label;
   inputMgr->GetRunList(allRuns);
   inputMgr->GetRunLabels(label);

   const int NRUN = allRuns.size();
   if(NRUN==0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   for(int i=0;i<NRUN;i++){
      if(label[i].compare("skew-quad")==0) vRun.push_back(allRuns[i]); 
      if(label[i].compare("norm-quad")==0) rRun.push_back(allRuns[i]); 
      if(label[i].compare("bare")==0)      bRun.push_back(allRuns[i]); 
      if(label[i].compare("bare-2")==0)    b2Run.push_back(allRuns[i]); 
      if(label[i].compare("azi")==0)       aRun.push_back(allRuns[i]); 
      if(label[i].compare("bare-2")!=0){
	 // save all runs that aren't the last bare run
	 run.push_back(allRuns[i]);
      }
   }

   std::cout << "Radial gradient runs: " << std::endl;
   const int NR = rRun.size();
   for(int i=0;i<NR;i++) std::cout << rRun[i] << std::endl;  

   std::cout << "Vertical gradient runs: " << std::endl;
   const int NV = vRun.size();
   for(int i=0;i<NV;i++) std::cout << vRun[i] << std::endl; 

   std::cout << "Azimuthal gradient runs: " << std::endl;
   const int NA = aRun.size();
   for(int i=0;i<NA;i++) std::cout << aRun[i] << std::endl; 

   std::cout << "Bare field runs: " << std::endl;
   const int NB = bRun.size();
   for(int i=0;i<NB;i++) std::cout << bRun[i] << std::endl; 

   std::cout << "Bare-2 field runs: " << std::endl;
   const int NB2 = b2Run.size();
   for(int i=0;i<NB2;i++) std::cout << b2Run[i] << std::endl; 

   std::vector<int> driftRun;
   for(int i=0;i<NB;i++)  driftRun.push_back(bRun[i]); 
   for(int i=0;i<NB2;i++) driftRun.push_back(b2Run[i]); 

   // Trolley data
   std::vector<trolleyAnaEvent_t> trlyData; // for drift  
   std::vector<trolleyAnaEvent_t> bareData,bareDataCor,bareDataCorAlt;    
   std::vector<trolleyAnaEvent_t> radData ,radDataCor ,radDataCorAlt;    
   std::vector<trolleyAnaEvent_t> vertData,vertDataCor,vertDataCorAlt;    

   for(int i=0;i<NR;i++) rc = GetTrolleyData(date,rRun[i],method,radData);
   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   for(int i=0;i<NV;i++) rc = GetTrolleyData(date,vRun[i],method,vertData);
   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   for(int i=0;i<NB;i++) rc = GetTrolleyData(date,bRun[i],method,bareData);
   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   if(isBlind){
      rc = ApplyBlindingTRLY(blindValue,radData);
      rc = ApplyBlindingTRLY(blindValue,vertData);
      rc = ApplyBlindingTRLY(blindValue,bareData);
   }

   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   // fixed probe data
   std::string fxprPath = "./input/probe-lists/fxpr-list.csv";
   std::vector<int> fxprList;
   gm2fieldUtil::Import::ImportData1<int>(fxprPath,"csv",fxprList);

   // get graph of fxpr data for drift correction *across* runs 
   std::vector<double> stats;
   TGraph *gDrift = GetDriftTGraph(method,driftRun,fxprList,stats); 

   // for drift correction *during* a run
   const int NZ = run.size();
   std::vector<gm2field::fixedProbeFrequency_t> fxprData;
   for(int i=0;i<NZ;i++) rc = gm2fieldUtil::RootHelper::GetFPFrequencies(run[i],fxprData);
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
   for(int i=0;i<NZ;i++) rc = GetTrolleyData(date,run[i],method,trlyData); 
   if(rc!=0) return 1;

   // remove the trolley probe of interest 
   int NL = inList.size();
   for(int i=0;i<NL;i++) if(probeNumber!=inList[i]) trlyList.push_back(inList[i]); 

   std::cout << "Applying field drift corrections (during measurement)..." << std::endl;
   // correct for field drift 
   // use fxpr 
   rc = CorrectTRLYForDriftDuringMeasurement(fxprDataAvg,bareData,bareDataCor);
   rc = CorrectTRLYForDriftDuringMeasurement(fxprDataAvg,radData ,radDataCor);
   rc = CorrectTRLYForDriftDuringMeasurement(fxprDataAvg,vertData,vertDataCor);
   if(rc!=0) return 1;
   // use trly 
   rc = CorrectTRLYForDriftDuringMeasurement(method,trlyList,trlyData,bareData,bareDataCorAlt);
   rc = CorrectTRLYForDriftDuringMeasurement(method,trlyList,trlyData,radData ,radDataCorAlt);
   rc = CorrectTRLYForDriftDuringMeasurement(method,trlyList,trlyData,vertData,vertDataCorAlt);
   if(rc!=0) return 1;
   std::cout << "--> Done" << std::endl; 

   // correct the radial data 
   double drift_rad[3]     = {0,0,0}; 
   double drift_rad_err[3] = {0,0,0}; 
   const int NBC = bareDataCor.size(); 
   tStart = bareDataCor[NBC-1].time[NUM_TRLY-1];  // last event, last probe time  
   tStop  = radDataCor[0].time[0];               // first event, first probe time 
   rc = CorrectTRLYForDriftAcrossRuns(tStart,tStop,gDrift,radDataCor,drift_rad[1],drift_rad_err[1]); 
   if(rc!=0) return 1;

   // correct the vertical data -- need to redefine the stop time to 
   // sync up with the start of the vertical data set 
   double drift_vert[3]     = {0,0,0}; 
   double drift_vert_err[3] = {0,0,0}; 
   tStop = vertDataCor[0].time[0]; 
   rc = CorrectTRLYForDriftAcrossRuns(tStart,tStop,gDrift,vertDataCor,drift_vert[1],drift_vert_err[1]); 
   if(rc!=0) return 1;

   // find gradients 

   // radial
   TGraph *gRB         = GetSlicePlot('r',bareData); 
   TGraph *gRG         = GetSlicePlot('r',radData);
   TGraph *gRDiff      = GetDiffPlot(gRB,gRG); 

   TGraph *gRB_fxpr    = GetSlicePlot('r',bareDataCor); 
   TGraph *gRG_fxpr    = GetSlicePlot('r',radDataCor);
   TGraph *gRDiff_fxpr = GetDiffPlot(gRB_fxpr,gRG_fxpr); 

   TGraph *gRB_trly    = GetSlicePlot('r',bareDataCorAlt); 
   TGraph *gRG_trly    = GetSlicePlot('r',radDataCorAlt);
   TGraph *gRDiff_trly = GetDiffPlot(gRB_trly,gRG_trly); 

   // vertical  
   TGraph *gVB         = GetSlicePlot('v',bareData); 
   TGraph *gVG         = GetSlicePlot('v',vertData);
   TGraph *gVDiff      = GetDiffPlot(gVB,gVG); 

   TGraph *gVB_fxpr    = GetSlicePlot('v',bareDataCor); 
   TGraph *gVG_fxpr    = GetSlicePlot('v',vertDataCor);
   TGraph *gVDiff_fxpr = GetDiffPlot(gVB_fxpr,gVG_fxpr); 

   TGraph *gVB_trly    = GetSlicePlot('v',bareDataCorAlt); 
   TGraph *gVG_trly    = GetSlicePlot('v',vertDataCorAlt);
   TGraph *gVDiff_trly = GetDiffPlot(gVB_trly,gVG_trly); 

   gm2fieldUtil::Graph::SetGraphParameters(gVB   ,20,kBlack); 
   gm2fieldUtil::Graph::SetGraphParameters(gVG   ,20,kRed); 
   gm2fieldUtil::Graph::SetGraphParameters(gVDiff,20,kBlack); 

   gm2fieldUtil::Graph::SetGraphParameters(gVB_fxpr   ,20,kBlack); 
   gm2fieldUtil::Graph::SetGraphParameters(gVG_fxpr   ,20,kRed); 
   gm2fieldUtil::Graph::SetGraphParameters(gVDiff_fxpr,20,kBlack); 

   gm2fieldUtil::Graph::SetGraphParameters(gVB_trly   ,20,kBlack); 
   gm2fieldUtil::Graph::SetGraphParameters(gVG_trly   ,20,kRed); 
   gm2fieldUtil::Graph::SetGraphParameters(gVDiff_trly,20,kBlack); 

   gm2fieldUtil::Graph::SetGraphParameters(gRB   ,20,kBlack); 
   gm2fieldUtil::Graph::SetGraphParameters(gRG   ,20,kRed); 
   gm2fieldUtil::Graph::SetGraphParameters(gRDiff,20,kBlack); 

   gm2fieldUtil::Graph::SetGraphParameters(gRB_fxpr   ,20,kBlack); 
   gm2fieldUtil::Graph::SetGraphParameters(gRG_fxpr   ,20,kRed); 
   gm2fieldUtil::Graph::SetGraphParameters(gRDiff_fxpr,20,kBlack); 

   gm2fieldUtil::Graph::SetGraphParameters(gRB_trly   ,20,kBlack); 
   gm2fieldUtil::Graph::SetGraphParameters(gRG_trly   ,20,kRed); 
   gm2fieldUtil::Graph::SetGraphParameters(gRDiff_trly,20,kBlack); 

   TMultiGraph *mgr = new TMultiGraph();
   mgr->Add(gRB,"p"); 
   mgr->Add(gRG,"p"); 

   TMultiGraph *mgr_fxpr = new TMultiGraph();
   mgr_fxpr->Add(gRB_fxpr,"p"); 
   mgr_fxpr->Add(gRG_fxpr,"p"); 
 
   TMultiGraph *mgr_trly = new TMultiGraph();
   mgr_trly->Add(gRB_trly,"p"); 
   mgr_trly->Add(gRG_trly,"p"); 
 
   TMultiGraph *mgv = new TMultiGraph();
   mgv->Add(gVB,"p"); 
   mgv->Add(gVG,"p"); 

   TMultiGraph *mgv_fxpr = new TMultiGraph();
   mgv_fxpr->Add(gVB_fxpr,"p"); 
   mgv_fxpr->Add(gVG_fxpr,"p"); 
 
   TMultiGraph *mgv_trly = new TMultiGraph();
   mgv_trly->Add(gVB_trly,"p"); 
   mgv_trly->Add(gVG_trly,"p"); 
 
   TCanvas *c1 = new TCanvas("c1","Radial Gradients",1200,800);
   c1->Divide(3,2);

   c1->cd(1); 
   mgr->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mgr,"Radial Gradient","r (cm)","Frequency (Hz)");
   mgr->GetXaxis()->SetLimits(-4.5,4.5); 
   mgr->Draw("a");
   c1->Update();

   c1->cd(2); 
   mgr_fxpr->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mgr_fxpr,"Radial Gradient (Drift Cor FXPR)","r (cm)","Frequency (Hz)");
   mgr_fxpr->GetXaxis()->SetLimits(-4.5,4.5); 
   mgr_fxpr->Draw("a");
   c1->Update();

   c1->cd(3); 
   mgr_trly->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mgr_trly,"Radial Gradient (Drift Cor TRLY)","r (cm)","Frequency (Hz)");
   mgr_trly->GetXaxis()->SetLimits(-4.5,4.5); 
   mgr_trly->Draw("a");
   c1->Update();

   c1->cd(4); 
   gRDiff->Draw("ap");
   gm2fieldUtil::Graph::SetGraphLabels(gRDiff,"Radial Gradient","r (cm)","Frequency (Hz)");
   gRDiff->GetXaxis()->SetLimits(-4.5,4.5); 
   gRDiff->Draw("ap");
   gRDiff->Fit(fitFunc.c_str(),"Q");
   c1->Update();

   c1->cd(5); 
   gRDiff_fxpr->Draw("ap");
   gm2fieldUtil::Graph::SetGraphLabels(gRDiff_fxpr,"Radial Gradient (Drift Cor FXPR)","r (cm)","Frequency (Hz)");
   gRDiff_fxpr->GetXaxis()->SetLimits(-4.5,4.5); 
   gRDiff_fxpr->Draw("ap");
   gRDiff_fxpr->Fit(fitFunc.c_str(),"Q");
   c1->Update();

   c1->cd(6); 
   gRDiff_trly->Draw("ap");
   gm2fieldUtil::Graph::SetGraphLabels(gRDiff_trly,"Radial Gradient (Drift Cor TRLY)","r (cm)","Frequency (Hz)");
   gRDiff_trly->GetXaxis()->SetLimits(-4.5,4.5); 
   gRDiff_trly->Draw("ap");
   gRDiff_trly->Fit(fitFunc.c_str(),"Q");
   c1->Update();

   TString rPlotPath = Form("%s/rad-grad_bare-run-%d_grad-run-%d_%s.png",plotDir,bRun[0],rRun[0],date.c_str());
   c1->cd(); 
   c1->Print(rPlotPath);  

   TCanvas *c2 = new TCanvas("c2","Vertical Gradients",1200,800);
   c2->Divide(3,2);

   c2->cd(1); 
   mgv->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mgv,"Vertical Gradient","y (cm)","Frequency (Hz)");
   mgv->GetXaxis()->SetLimits(-4.5,4.5); 
   mgv->Draw("a");
   c2->Update();

   c2->cd(2); 
   mgv_fxpr->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mgv_fxpr,"Vertical Gradient (Drift Cor FXPR)","y (cm)","Frequency (Hz)");
   mgv_fxpr->GetXaxis()->SetLimits(-4.5,4.5); 
   mgv_fxpr->Draw("a");
   c2->Update();

   c2->cd(3); 
   mgv_trly->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mgv_trly,"Vertical Gradient (Drift Cor TRLY)","y (cm)","Frequency (Hz)");
   mgv_trly->GetXaxis()->SetLimits(-4.5,4.5); 
   mgv_trly->Draw("a");
   c2->Update();

   c2->cd(4); 
   gVDiff->Draw("ap");
   gm2fieldUtil::Graph::SetGraphLabels(gVDiff,"Vertical Gradient","r (cm)","Frequency (Hz)");
   gVDiff->GetXaxis()->SetLimits(-4.5,4.5); 
   gVDiff->Draw("ap");
   gVDiff->Fit(fitFunc.c_str(),"Q");
   c2->Update();

   c2->cd(5); 
   gVDiff_fxpr->Draw("ap");
   gm2fieldUtil::Graph::SetGraphLabels(gVDiff_fxpr,"Vertical Gradient (Drift Cor FXPR)","y (cm)","Frequency (Hz)");
   gVDiff_fxpr->GetXaxis()->SetLimits(-4.5,4.5); 
   gVDiff_fxpr->Draw("ap");
   gVDiff_fxpr->Fit(fitFunc.c_str(),"Q");
   c2->Update();

   c2->cd(6); 
   gVDiff_trly->Draw("ap");
   gm2fieldUtil::Graph::SetGraphLabels(gVDiff_trly,"Vertical Gradient (Drift Cor TRLY)","y (cm)","Frequency (Hz)");
   gVDiff_trly->GetXaxis()->SetLimits(-4.5,4.5); 
   gVDiff_trly->Draw("ap");
   gVDiff_trly->Fit(fitFunc.c_str(),"Q");
   c2->Update();

   TString vPlotPath = Form("%s/vert-grad_bare-run-%d_grad-run-%d_%s.png",plotDir,bRun[0],vRun[0],date.c_str());
   c2->cd(); 
   c2->Print(vPlotPath);  

   double PR[3],ER[3];
   TF1 *fitR      = gRDiff->GetFunction(fitFunc.c_str());
   TF1 *fitR_fxpr = gRDiff_fxpr->GetFunction(fitFunc.c_str());
   TF1 *fitR_trly = gRDiff_trly->GetFunction(fitFunc.c_str());

   double PV[3],EV[3];
   TF1 *fitV      = gVDiff->GetFunction(fitFunc.c_str());
   TF1 *fitV_fxpr = gVDiff_fxpr->GetFunction(fitFunc.c_str());
   TF1 *fitV_trly = gVDiff_trly->GetFunction(fitFunc.c_str());

   PR[0] = fitR->GetParameter(1);
   PR[1] = fitR_fxpr->GetParameter(1);
   PR[2] = fitR_trly->GetParameter(1);

   ER[0] = fitR->GetParError(1);
   ER[1] = fitR_fxpr->GetParError(1);
   ER[2] = fitR_trly->GetParError(1);

   PV[0] = fitV->GetParameter(1);
   PV[1] = fitV_fxpr->GetParameter(1);
   PV[2] = fitV_trly->GetParameter(1);

   EV[0] = fitV->GetParError(1);
   EV[1] = fitV_fxpr->GetParError(1);
   EV[2] = fitV_trly->GetParError(1);

   // convert to Hz/mm 
   for(int i=0;i<3;i++){
      PR[i] /= 10.;
      ER[i] /= 10.;
      PV[i] /= 10.;
      EV[i] /= 10.;
   }

   // estimate the error due to drift as the difference 
   // between uncorrected data and drift-corrected data
   // adding in the stdev of the anchor points as an additional uncertainty  
   drift_rad_err[1]  = TMath::Sqrt( TMath::Power(PR[1]-PR[0],2.) + TMath::Power(stats[1],2.) + TMath::Power(stats[3],2.)); 
   drift_rad_err[2]  = TMath::Sqrt( TMath::Power(PR[2]-PR[0],2.) + TMath::Power(stats[1],2.) + TMath::Power(stats[3],2.)); 
   drift_vert_err[1] = TMath::Sqrt( TMath::Power(PV[1]-PV[0],2.) + TMath::Power(stats[1],2.) + TMath::Power(stats[3],2.)); 
   drift_vert_err[2] = TMath::Sqrt( TMath::Power(PV[2]-PV[0],2.) + TMath::Power(stats[1],2.) + TMath::Power(stats[3],2.)); 

   std::cout << "============================ RESULTS ============================" << std::endl; 
   std::cout << "Radial:          " << Form("%.3lf +/- %.3lf Hz/mm",PR[0],ER[0]) << std::endl;
   std::cout << "Radial [fxpr]:   " << Form("%.3lf +/- %.3lf +/- %.3lf Hz/mm",PR[1],ER[1],drift_rad_err[1]) << std::endl;
   std::cout << "Radial [trly]:   " << Form("%.3lf +/- %.3lf +/- %.3lf Hz/mm",PR[2],ER[2],drift_rad_err[2]) << std::endl;
   std::cout << "-----------------------------------------------------------------" << std::endl;
   std::cout << "Vertical:        " << Form("%.3lf +/- %.3lf Hz/mm",PV[0],EV[0]) << std::endl;
   std::cout << "Vertical [fxpr]: " << Form("%.3lf +/- %.3lf +/- %.3lf Hz/mm",PV[1],EV[1],drift_vert_err[1]) << std::endl;
   std::cout << "Vertical [trly]: " << Form("%.3lf +/- %.3lf +/- %.3lf Hz/mm",PV[2],EV[2],drift_vert_err[2]) << std::endl;

   char outpath[200],outdir[200];
   if(isBlind)  sprintf(outdir,"./output/blinded/%02d-%02d-%02d"  ,theDate.month,theDate.day,theDate.year-2000); 
   if(!isBlind) sprintf(outdir,"./output/unblinded/%02d-%02d-%02d",theDate.month,theDate.day,theDate.year-2000);
   if(!isFullAnalysis){
      // we're doing a delta-b calculation during a measurement
      // store in a different location 
      sprintf(outdir,"./output/delta-b/%02d-%02d-%02d",theDate.month,theDate.day,theDate.year-2000);
   } 
   rc = MakeDirectory(outdir);
   // radial 
   sprintf(outpath,"%s/rad-grad_pr-%02d_%s.csv"  ,outdir,probeNumber,date.c_str()); 
   PrintToFile(outpath,"rad-grad",PR,ER,drift_rad,drift_rad_err); 
   // vertical 
   sprintf(outpath,"%s/vert-grad_pr-%02d_%s.csv"  ,outdir,probeNumber,date.c_str()); 
   PrintToFile(outpath,"vert-grad",PV,EV,drift_vert,drift_vert_err); 

   delete inputMgr; 

   return 0;
}

