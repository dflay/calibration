// Get shimmed field for trolley data (stationary run)   

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

int GetShimmedField_trly(std::string configFile){

   std::cout << "----------------------------------" << std::endl;
   std::cout << "TRLY SHIMMED FIELD CALCULATION" << std::endl;

   int rc=0;
   int method = gm2fieldUtil::Constants::kPhaseDerivative;

   InputManager *inputMgr = new InputManager();
   inputMgr->Load(configFile);
   inputMgr->Print();

   std::string date    = inputMgr->GetAnalysisDate();
   bool isBlind        = inputMgr->IsBlind();
   bool isFullAnalysis = inputMgr->IsFullAnalysis();
   bool useP2PFit      = inputMgr->UseP2PFit();
   int probeNumber     = inputMgr->GetTrolleyProbe();

   date_t theDate;
   GetDate(theDate);

   char plotDir[200];
   sprintf(plotDir,"./plots/%02d-%02d-%02d",theDate.month,theDate.day,theDate.year-2000);
   MakeDirectory(plotDir); 

   blind_t blind;
   ImportBlinding(blind);
   double blindValue = blind.value_tr; 

   // char inpath[200];
   // sprintf(inpath,"./input/runlists/%s/trly-shimmed_%s.csv",date.c_str(),date.c_str());
   // std::vector<int> run;
   // std::vector<double> sf;
   // std::vector<std::string> label;
   // ImportDeltaBFileList_csv(inpath,run,label,sf);

   std::vector<int> run;
   std::vector<std::string> label;
   inputMgr->GetRunList(run);
   inputMgr->GetRunLabels(label);

   const int NRUN = run.size();
   if(NRUN==0){
      std::cout << "No data!" << std::endl;
      return 1;
   }
 
   // int bareRun=0,bareIndex=0,gradRun=0,gradIndex=0;
   // for(int i=0;i<NRUN;i++){
   //    if(label[i].compare("bare")==0){
   //       bareRun   = run[i];
   //       bareIndex = i;
   //    }else{
   //       gradRun   = run[i];
   //       gradIndex = i;
   //    }
   // } 

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

   int NN  = fxprData.size();
   int NFP = fxprList.size();
   int startIndex = fxprList[0];
   int stopIndex  = fxprList[NFP-1];
   unsigned long long t0     = 0;
   unsigned long long tStart = fxprData[0].GpsTimeStamp[startIndex];
   unsigned long long tStop  = fxprData[NN-1].GpsTimeStamp[stopIndex];
   unsigned long long tStep  = 1E+9;

   std::vector<fixedProbeEvent_t> fxprDataAvg;
   GetAverageFXPRVectorsNew(method,t0,tStart,tStop,tStep,fxprList,fxprData,fxprDataAvg);

   // Fixed probe average plot 
   TGraph *gFPAVG = GetTGraphNew(fxprDataAvg);
   gm2fieldUtil::Graph::SetGraphParameters(gFPAVG,20,kBlack);

   // get a fit function to these FXPR data
   const int nOrder = 4;
   std::vector< std::vector<double> > eps;
   TF1 *fxprFit = GetPolyFitToFXPR("myFit",nOrder,gFPAVG,eps);

   double mean_p2p=0;
   double stdev_p2p=0;
   double STDEV_P2P=0;
   TGraph *gFPAVG_trm = RemoveTrend(gFPAVG,fxprFit,mean_p2p,STDEV_P2P);

   if(useP2PFit){
      stdev_p2p = STDEV_P2P;
   }

   // Trolley data 
   std::vector<trolleyAnaEvent_t> trlyData;   // for drift correction  
   std::vector<trolleyAnaEvent_t> Event,EventCor,EventCorAlt; 

   for(int i=0;i<NRUN;i++){
      std::cout << "Getting run " << run[i] << std::endl; 
      rc = GetTrolleyData(date,run[i],method,Event);
      if(rc!=0){
	 std::cout << "No data!" << std::endl;
	 return 1;
      }
   }

   if(isBlind) rc = ApplyBlindingTRLY(blindValue,Event); 
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
   if(useP2PFit){
      rc = CorrectTRLYForDriftDuringMeasurement(method,fxprFit,Event,EventCor);
   }else{
      rc = CorrectTRLYForDriftDuringMeasurement(fxprDataAvg,Event,EventCor);
   }
   if(rc!=0) return 1;
   // use trly 
   rc = CorrectTRLYForDriftDuringMeasurement(method,trlyList,trlyData,Event,EventCorAlt);
   if(rc!=0) return 1;
   std::cout << "--> Done" << std::endl; 

   // now do averaging  
   double B[3]      = {0,0,0}; 
   double B_err[3]  = {0,0,0};
   double P2P[3]    = {0,0,0};  
   double P2PErr[3] = {0,0,0};  

   // raw 
   CalculateTRLYAvg_Stationary(probeNumber,Event      ,B[0],B_err[0]);
   // drift corrected (fxpr)  
   CalculateTRLYAvg_Stationary(probeNumber,EventCor   ,B[1],B_err[1]);
   // drift corrected (trly)   
   CalculateTRLYAvg_Stationary(probeNumber,EventCorAlt,B[2],B_err[2]);

   if(useP2PFit){
      for(int i=1;i<3;i++) P2PErr[i] = STDEV_P2P; 
   }

   std::cout << Form("===================== RESULTS FOR PROBE %02d =====================",probeNumber) << std::endl;
   std::cout << "Raw results: " << std::endl; 
   std::cout << Form("%.3lf +/- %.3lf Hz (%.3lf ppb)",
                     B[0],B_err[0],B_err[0]/0.06179) << std::endl;
   std::cout << "-----------------------------------------------------------------------" << std::endl;  
   std::cout << "Drift corrected [FXPR]: " << std::endl;
   std::cout << Form("%.3lf +/- %.3lf Hz (%.3lf ppb) +/- %.3lf (%.3lf ppb)",
                     B[1],B_err[1],B_err[1]/0.06179,P2PErr[1],P2PErr[1]/0.06179) << std::endl;
   std::cout << "-----------------------------------------------------------------------" << std::endl;  
   std::cout << "Drift corrected [TRLY]: " << std::endl;
   std::cout << Form("%.3lf +/- %.3lf Hz (%.3lf ppb)",
                     B[2],B_err[2],B_err[2]/0.06179) << std::endl;

   char outpath[200],outdir[200];
   if(isBlind)  sprintf(outdir,"./output/blinded/%02d-%02d-%02d"  ,theDate.month,theDate.day,theDate.year-2000); 
   if(!isBlind) sprintf(outdir,"./output/unblinded/%02d-%02d-%02d",theDate.month,theDate.day,theDate.year-2000); 
   rc = MakeDirectory(outdir);
   sprintf(outpath,"%s/trly_shimmed-field_pr-%02d_%s.csv",outdir,probeNumber,date.c_str()); 
   rc = PrintToFile(outpath,"shim",B,B_err,P2P,P2PErr); 

   // make some plots 

   // Trolley plots 
   TGraph *g1 = GetTRLYTGraph(probeNumber,"GpsTimeStamp","freq",Event);
   gm2fieldUtil::Graph::SetGraphParameters(g1,20,kBlack);

   TGraph *g2 = GetTRLYTGraph(probeNumber,"GpsTimeStamp","freq",EventCor);
   gm2fieldUtil::Graph::SetGraphParameters(g2,20,kRed);

   TMultiGraph *mgTRLY = new TMultiGraph();
   mgTRLY->Add(g1,"lp");
   mgTRLY->Add(g2,"lp");

   // Fixed probe plot 
   // const int N = Event.size();
   // TGraph *gFPAVG = GetTGraph(method,Event[0].time[0],Event[N-1].time[3],1E+9,fxprList,fxprData);
   // gm2fieldUtil::Graph::SetGraphParameters(gFPAVG,20,kBlack);

   TLegend *LT = new TLegend(0.6,0.6,0.8,0.8);
   LT->AddEntry(g1,"Raw","p");
   LT->AddEntry(g2,"Drift Corrected (FXPR)","p");

   double xMin = tStart/1E+9;
   double xMax = tStop/1E+9;

   TCanvas *c2 = new TCanvas("c2","TRLY & FXPR Data",1200,600);
   c2->Divide(1,2);

   c2->cd(1);
   mgTRLY->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mgTRLY,"TRLY Data","","Frequency (Hz)");
   gm2fieldUtil::Graph::UseTimeDisplay(mgTRLY);
   gm2fieldUtil::Graph::SetGraphLabelSizes(mgTRLY,0.05,0.06);
   mgTRLY->GetXaxis()->SetLimits(xMin,xMax); 
   mgTRLY->Draw("a");
   LT->Draw("same");
   c2->Update();

   c2->cd(2);
   gFPAVG->Draw("alp");
   gm2fieldUtil::Graph::SetGraphLabels(gFPAVG,"Fixed Probe Average","","Frequency (Hz)");
   gm2fieldUtil::Graph::UseTimeDisplay(gFPAVG);
   gm2fieldUtil::Graph::SetGraphLabelSizes(gFPAVG,0.05,0.06);
   gFPAVG->GetXaxis()->SetLimits(xMin,xMax); 
   gFPAVG->Draw("alp");
   if(useP2PFit) fxprFit->Draw("same"); 
   c2->Update();

   TString plotPath = Form("%s/trly-shimmed_drift-cor_pr-%02d_run-%d_%s.png",plotDir,probeNumber,run[0],date.c_str());

   c2->cd();
   c2->Print(plotPath);

   delete inputMgr; 

   return 0;
}

