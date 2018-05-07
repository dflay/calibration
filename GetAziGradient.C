// Compute imposed azimuthal gradient using trolley data   

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

int GetAziGradient(std::string date,std::string fitFuncStr,int probeNumber,bool isBlind,bool isFullAnalysis){

   TString fitFunc = Form("%s",fitFuncStr.c_str());
   // TString fitFunc = Form("pol1");

   int rc=0;
   int method = gm2fieldUtil::Constants::kPhaseDerivative;

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

   std::cout << "----------------------------------" << std::endl;
   std::cout << "SCC AZIMUTHAL GRADIENT CALCULATION" << std::endl;

   char inpath[200];
   sprintf(inpath,"./input/runlists/%s/azi-grad_%s.csv",date.c_str(),date.c_str());

   std::vector<int> allRuns,run,aRun,bRun,b2Run;
   std::vector<double> sf;
   std::vector<std::string> label;
   ImportDeltaBFileList_csv(inpath,allRuns,label,sf);

   const int NRUN = allRuns.size();
   if(NRUN==0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   for(int i=0;i<NRUN;i++){
      if(label[i].compare("bare")==0)   bRun.push_back(allRuns[i]); 
      if(label[i].compare("bare-2")==0) b2Run.push_back(allRuns[i]); 
      if(label[i].compare("azi")==0)    aRun.push_back(allRuns[i]);
      if(label[i].compare("bare-2")!=0){
	 run.push_back(allRuns[i]); 
      } 
   }
  
   const int NZ = run.size(); 
 
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
   std::vector<trolleyAnaEvent_t> trlyData;    
   std::vector<trolleyAnaEvent_t> aziData,bareData;    
   std::vector<trolleyAnaEvent_t> bareDataCor,bareDataCorAlt;    
   std::vector<trolleyAnaEvent_t> aziDataCor,aziDataCorAlt;    

   for(int i=0;i<NB;i++) rc = GetTrolleyData(date,bRun[i],method,bareData);
   if(rc!=0){
      std::cout << "No TRLY data!" << std::endl;
      return 1;
   }

   for(int i=0;i<NA;i++) rc = GetTrolleyData(date,aRun[i],method,aziData);
   if(rc!=0){
      std::cout << "No TRLY data!" << std::endl;
      return 1;
   }

   if(isBlind){
      rc = ApplyBlindingTRLY(blindValue,bareData);
      rc = ApplyBlindingTRLY(blindValue,aziData);
   }

   if(rc!=0){
      std::cout << "No TRLY data!" << std::endl;
      return 1;
   }

   // fixed probe data
   std::string fxprPath = "./input/probe-lists/fxpr-list.csv";
   std::vector<int> fxprList;
   gm2fieldUtil::Import::ImportData1<int>(fxprPath,"csv",fxprList);

   // get graph of fxpr data for drift correction *across* runs 
   std::vector<double> stats;
   TGraph *gDrift = GetDriftTGraph(method,driftRun,fxprList,stats);

   std::vector<gm2field::fixedProbeFrequency_t> fxprData;
   for(int i=0;i<NZ;i++) rc = gm2fieldUtil::RootHelper::GetFPFrequencies(run[i],fxprData);
   if(rc!=0){
      std::cout << "No FXPR data!" << std::endl;
      return 1;
   }
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
   rc = CorrectTRLYForDriftDuringMeasurement(method,fxprList,fxprData,bareData,bareDataCor);
   rc = CorrectTRLYForDriftDuringMeasurement(method,fxprList,fxprData,aziData ,aziDataCor);
   if(rc!=0) return 1;
   // use trly 
   rc = CorrectTRLYForDriftDuringMeasurement(method,trlyList,trlyData,bareData,bareDataCorAlt);
   rc = CorrectTRLYForDriftDuringMeasurement(method,trlyList,trlyData,aziData ,aziDataCorAlt);
   if(rc!=0) return 1;
   std::cout << "--> Done" << std::endl;

   // correct the azimuthal data
   double drift[3]     = {0,0,0};
   double drift_err[3] = {0,0,0}; 
   const int NN = bareDataCor.size();
   unsigned long long tStart = bareDataCor[NN-1].time[NUM_TRLY-1];  // last event, last probe time  
   unsigned long long tStop  = aziDataCor[0].time[0];               // first event, first probe time 
   rc = CorrectTRLYForDriftAcrossRuns(tStart,tStop,gDrift,aziDataCor,drift[1],drift_err[1]);
   if(rc!=0) return 1;

   // find gradient 
   TGraph *gPRA  = GetTRLYTGraph(0,"phi","freq",aziData);
   TGraph *gPRB  = GetTRLYTGraph(0,"phi","freq",bareData);
   TGraph *gDiff = GetDiffPlot(gPRB,gPRA); 

   TGraph *gPRA_fxpr  = GetTRLYTGraph(0,"phi","freq",aziDataCor);
   TGraph *gPRB_fxpr  = GetTRLYTGraph(0,"phi","freq",bareDataCor);
   TGraph *gDiff_fxpr = GetDiffPlot(gPRB_fxpr,gPRA_fxpr); 

   TGraph *gPRA_trly  = GetTRLYTGraph(0,"phi","freq",aziDataCorAlt);
   TGraph *gPRB_trly  = GetTRLYTGraph(0,"phi","freq",bareDataCorAlt);
   TGraph *gDiff_trly = GetDiffPlot(gPRB_trly,gPRA_trly); 

   gm2fieldUtil::Graph::SetGraphParameters(gPRA ,20,kRed);
   gm2fieldUtil::Graph::SetGraphParameters(gPRB ,20,kBlack);
   gm2fieldUtil::Graph::SetGraphParameters(gDiff,20,kBlack);

   gm2fieldUtil::Graph::SetGraphParameters(gPRA_fxpr ,20,kRed);
   gm2fieldUtil::Graph::SetGraphParameters(gPRB_fxpr ,20,kBlack);
   gm2fieldUtil::Graph::SetGraphParameters(gDiff_fxpr,20,kBlack);

   gm2fieldUtil::Graph::SetGraphParameters(gPRA_trly ,20,kRed);
   gm2fieldUtil::Graph::SetGraphParameters(gPRB_trly ,20,kBlack);
   gm2fieldUtil::Graph::SetGraphParameters(gDiff_trly,20,kBlack);

   TMultiGraph *mg = new TMultiGraph();
   mg->Add(gPRB,"p"); 
   mg->Add(gPRA,"p"); 

   TMultiGraph *mg2 = new TMultiGraph();
   mg2->Add(gPRB_fxpr,"p"); 
   mg2->Add(gPRA_fxpr,"p"); 

   TMultiGraph *mg3 = new TMultiGraph();
   mg3->Add(gPRB_trly,"p"); 
   mg3->Add(gPRA_trly,"p"); 

   double xmin = 188.9;
   double xmax = 189.3;

   TCanvas *c1 = new TCanvas("c1","Azimuthal Gradient",1200,600);
   c1->Divide(3,2); 

   c1->cd(1);
   mg->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mg,"Azimuthal Field","Azimuthal Position (deg)","Frequency (Hz)");
   mg->Draw("a");
   c1->Update(); 

   c1->cd(2);
   mg2->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mg2,"Azimuthal Field (FXPR Drift Cor)","Azimuthal Position (deg)","Frequency (Hz)");
   mg2->Draw("a");
   c1->Update(); 

   c1->cd(3);
   mg3->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mg3,"Azimuthal Field (TRLY Drift Cor)","Azimuthal Position (deg)","Frequency (Hz)");
   mg3->Draw("a");
   c1->Update(); 

   c1->cd(4);
   gStyle->SetOptFit(111);
   gDiff->Draw("ap");
   gm2fieldUtil::Graph::SetGraphLabels(gDiff,"Azimuthal Gradient","Azimuthal Position (deg)","Frequency (Hz)");
   gDiff->Draw("ap");
   gDiff->Fit(fitFunc,"QR","",xmin,xmax);
   c1->Update(); 

   c1->cd(5);
   gStyle->SetOptFit(111);
   gDiff_fxpr->Draw("ap");
   gm2fieldUtil::Graph::SetGraphLabels(gDiff_fxpr,"Azimuthal Gradient (FXPR Drift Cor)","Azimuthal Position (deg)","Frequency (Hz)");
   gDiff_fxpr->Draw("ap");
   gDiff_fxpr->Fit(fitFunc,"QR","",xmin,xmax);
   c1->Update(); 

   c1->cd(6);
   gStyle->SetOptFit(111);
   gDiff_trly->Draw("ap");
   gm2fieldUtil::Graph::SetGraphLabels(gDiff_trly,"Azimuthal Gradient (TRLY Drift Cor)","Azimuthal Position (deg)","Frequency (Hz)");
   gDiff_trly->Draw("ap");
   gDiff_trly->Fit(fitFunc,"QR","",xmin,xmax);
   c1->Update(); 
 
   TString plotPath = Form("%s/azi-grad_%s.png",plotDir,date.c_str());

   c1->cd();
   c1->Print(plotPath); 
 
   TF1 *fitZ      = gDiff->GetFunction(fitFunc);
   TF1 *fitZ_fxpr = gDiff_fxpr->GetFunction(fitFunc);
   TF1 *fitZ_trly = gDiff_trly->GetFunction(fitFunc);
   double P[3],err[3];
   P[0] = fitZ->GetParameter(1); 
   P[1] = fitZ_fxpr->GetParameter(1); 
   P[2] = fitZ_trly->GetParameter(1);

   err[0] = fitZ->GetParError(1); 
   err[1] = fitZ_fxpr->GetParError(1); 
   err[2] = fitZ_trly->GetParError(1);

   // convert to Hz/mm
   for(int i=0;i<3;i++){
      P[i]   /= 124.;
      err[i] /= 124.;
   }

   // estimate drift error by taking difference in drift-corrected
   // and uncorrected data
   // adding in the stdev of the anchor points as an additional uncertainty  
   drift_err[1] = TMath::Sqrt( TMath::Power(P[1]-P[0],2.) + TMath::Power(stats[1],2.) + TMath::Power(stats[3],2.)); 
   drift_err[2] = TMath::Sqrt( TMath::Power(P[2]-P[0],2.) + TMath::Power(stats[1],2.) + TMath::Power(stats[3],2.)); 

   std::cout << "====================================== RESULTS ======================================" << std::endl;
   std::cout << "azi gradient        = " << Form("%.3lf +/- %.3lf Hz/mm"          ,P[0],err[0])              << std::endl;
   std::cout << "azi gradient [fxpr] = " << Form("%.3lf +/- %.3lf +/- %.3lf Hz/mm",P[1],err[1],drift_err[1]) << std::endl;
   std::cout << "azi gradient [trly] = " << Form("%.3lf +/- %.3lf +/- %.3lf Hz/mm",P[2],err[2],drift_err[2]) << std::endl;

   char outpath[200],outdir[200];
   if(isBlind)  sprintf(outdir,"./output/blinded/%02d-%02d-%02d"  ,theDate.month,theDate.day,theDate.year-2000); 
   if(!isBlind) sprintf(outdir,"./output/unblinded/%02d-%02d-%02d",theDate.month,theDate.day,theDate.year-2000);
   if(!isFullAnalysis){
      // we're doing a delta-b calculation during a measurement
      // store in a different location 
      sprintf(outdir,"./output/delta-b/%02d-%02d-%02d",theDate.month,theDate.day,theDate.year-2000);
   }
 
   rc = MakeDirectory(outdir);
   sprintf(outpath,"%s/azi-grad_pr-%02d_%s.csv"  ,outdir,probeNumber,date.c_str());
   PrintToFile(outpath,"azi-grad",P,err,drift,drift_err);

   return 0;
}

