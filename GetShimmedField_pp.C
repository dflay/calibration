// Compute the shimmed field for the PP   

#include <cstdlib> 
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>  
#include <cmath> 

#include "TROOT.h"
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

#include "./src/FitFuncs.C"
#include "./src/FXPRFuncs.C"
#include "./src/Consolidate.C"
#include "./src/CustomMath.C"
#include "./src/CustomGraph.C"
#include "./src/CustomImport.C"
#include "./src/CustomExport.C"
#include "./src/CustomAlgorithms.C"
#include "./src/CustomUtilities.C"

int GetShimmedField_pp(std::string date,bool isBlind,bool useFXPRFit){

   int rc=0;
   int method = gm2fieldUtil::Constants::kPhaseDerivative;

   date_t theDate; 
   GetDate(theDate);

   char plotDir[200];
   sprintf(plotDir,"./plots/%02d-%02d-%02d",theDate.month,theDate.day,theDate.year-2000);
   MakeDirectory(plotDir); 

   blind_t blind;
   ImportBlinding(blind);
   double blindValue = blind.value_pp;

   std::cout << "----------------------------------" << std::endl;
   std::cout << "PP SHIMMED FIELD" << std::endl;
   
   char inpath[200];   
   sprintf(inpath,"./input/runlists/%s/pp-shimmed_%s.csv",date.c_str(),date.c_str());

   std::vector<int> run;
   std::vector<double> sf;
   std::vector<std::string> label;
   ImportDeltaBFileList_csv(inpath,run,label,sf);

   const int NRUN = run.size();

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

   double MEAN_P2P=0,STDEV_P2P=0;
   double mean_p2p=0,stdev_p2p=0;
   TGraph *gFPAVG_alt = RemoveTrend(gFPAVG,fxprFit,MEAN_P2P,STDEV_P2P);

   if(useFXPRFit){
      stdev_p2p = STDEV_P2P;
   }

   // TString tag;
   // TMultiGraph *mg = new TMultiGraph();
   // TLegend *L = new TLegend(0.6,0.6,0.8,0.8); 
   // const int NFP = fxprList.size();
   // TGraph **gFXPR = new TGraph*[NFP];
   // for(int i=0;i<NFP;i++){
   //    gFXPR[i] = gm2fieldUtil::Graph::GetTGraph(fxprList[i],method,"GpsTimeStamp","Frequency",fxprData);
   //    gm2fieldUtil::Graph::SetGraphParameters(gFXPR[i],20,kWhite+i+1);
   //    mg->Add(gFXPR[i],"p");
   //    tag = Form("fxpr %03d",fxprList[i]); 
   //    L->AddEntry(gFXPR[i],tag.Data(),"p"); 
   // }

   // TCanvas *c1 = new TCanvas("c1","Fixed and Trolley Probes",1200,600);
   // c1->Divide(1,2); 

   // c1->cd(1);
   // mg->Draw("a");
   // gm2fieldUtil::Graph::SetGraphLabels(mg,"Fixed Probe Data","","Frequency (Hz)");
   // gm2fieldUtil::Graph::UseTimeDisplay(mg);  
   // mg->Draw("a");
   // L->Draw("same");  
   // c1->Update();

   // trolley data 
   std::string trlyPath = "./input/probe-lists/trly-list.csv"; 
   std::vector<int> trlyList;
   gm2fieldUtil::Import::ImportData1<int>(trlyPath,"csv",trlyList);

   std::vector<trolleyAnaEvent_t> trlyData;
   for(int i=0;i<NRUN;i++) rc = GetTrolleyData(date,run[i],method,trlyData);

   int NTR = trlyData.size();
   std::cout << "Found " << NTR << " trolley events" << std::endl; 
  
   // if(NTR==0){
   //    std::cout << "ERROR: No trolley data!" << std::endl;
   //    return 1;
   // }

   // TMultiGraph *mg2 = new TMultiGraph();
   // TLegend *L2 = new TLegend(0.6,0.6,0.8,0.8); 
   // const int NTR = trlyList.size();
   // TGraph **gTRLY = new TGraph*[NTR];
   // for(int i=0;i<NTR;i++){
   //    gTRLY[i] = GetTRLYTGraph(trlyList[i],"GpsTimeStamp","freq",trlyData);
   //    gm2fieldUtil::Graph::SetGraphParameters(gTRLY[i],20,kWhite+i+1);
   //    mg2->Add(gTRLY[i],"p");
   //    tag = Form("trly %03d",trlyList[i]); 
   //    L2->AddEntry(gTRLY[i],tag.Data(),"p"); 
   // }

   // c1->cd(2);
   // mg2->Draw("a");
   // gm2fieldUtil::Graph::SetGraphLabels(mg2,"Trolley Probe Data","","Frequency (Hz)");
   // gm2fieldUtil::Graph::UseTimeDisplay(mg2);  
   // mg2->Draw("a");
   // L2->Draw("same");  
   // c1->Update();

   // PP data 
   std::vector<plungingProbeAnaEvent_t> ppEvent,ppEvent_tr,ppEventCor,ppEventCor_tr; 
   for(int i=0;i<NRUN;i++){
      std::cout << "Getting run " << run[i] << std::endl; 
      rc = GetPlungingProbeData(run[i],method,ppEvent);
      if(rc!=0){
	 std::cout << "No data!" << std::endl;
	 return 1;
      }
   }

   const int NPP = ppEvent.size();
   for(int i=0;i<NPP;i++) rc = ModifyPlungingProbeData(kLeastSquaresPhase,ppEvent[i]);

   if(isBlind){
      ApplyBlindingPP(blindValue,ppEvent);
   }
 
   CopyPlungingProbe(ppEvent,ppEvent_tr);

   if(useFXPRFit){
      rc = CorrectPPForDriftDuringMeasurement(method,fxprFit,ppEvent,ppEventCor);
   }else{
      rc = CorrectPPForDriftDuringMeasurement(fxprDataAvg,ppEvent,ppEventCor);
      // rc = CorrectPPForDriftDuringMeasurement(method,fxprList,fxprData,ppEvent,ppEventCor);
   }
   // rc = CorrectPPForDriftDuringMeasurement(method,trlyList,trlyData,ppEvent_tr,ppEventCor_tr);

   // now do averaging 
   double B[3]     = {0,0,0}; 
   double B_err[3] = {0,0,0}; 
   // raw 
   CalculateAveragePP(method,fxprList,fxprData,ppEvent   ,B[0],B_err[0]);
   // drift corrected [fxpr] 
   CalculateAveragePP(method,fxprList,fxprData,ppEventCor,B[1],B_err[1]);
   // drift corrected [trly] 
   // CalculateAveragePP(true ,method,trlyList,trlyData,ppEventCor_tr,B[2],B_err[2]);

   // FIXME: while trolley drift gives garbage:
   B[2]     = B[0];
   B_err[2] = B_err[0];  
   // for(int i=0;i<3;i++) std::cout << Form("%.3lf +/- %.3lf Hz",B[i],B_err[i]) << std::endl; 

   double P2P[3]    = {0,0,0};
   double P2PErr[3] = {0,stdev_p2p,stdev_p2p}; 

   std::cout << "====================== RESULTS ======================" << std::endl;
   std::cout << "Raw results: " << std::endl;
   std::cout << Form("%.3lf +/- %.3lf Hz (%.3lf ppb)",B[0],B_err[0],B_err[0]/0.06179) << std::endl;
   std::cout << "-----------------------------------------------------" << std::endl; 
   std::cout << "Drift corrected [FXPR]: " << std::endl;
   std::cout << Form("%.3lf +/- %.3lf Hz (%.3lf ppb) +/- %.3lf (%.3lf ppb)",B[1],B_err[1],B_err[1]/0.06179,P2PErr[1],P2PErr[1]/0.06179) << std::endl;
   std::cout << "-----------------------------------------------------" << std::endl; 
   std::cout << "Drift corrected [TRLY]: " << std::endl;
   std::cout << Form("%.3lf +/- %.3lf Hz (%.3lf ppb)",B[2],B_err[2],B_err[2]/0.06179) << std::endl;

   char outpath[200],outdir[200];
   if(isBlind)  sprintf(outdir,"./output/blinded/%02d-%02d-%02d"  ,theDate.month,theDate.day,theDate.year-2000); 
   if(!isBlind) sprintf(outdir,"./output/unblinded/%02d-%02d-%02d",theDate.month,theDate.day,theDate.year-2000); 
   rc = MakeDirectory(outdir);
   sprintf(outpath,"%s/pp_shimmed-field_%s.csv",outdir,date.c_str()); 
   rc = PrintToFile(outpath,"shim",B,B_err,P2P,P2PErr); 

   // make some plots 

   // Plunging probe plots 
   TGraph *g1 = GetPPTGraph1("TimeStamp","freq",ppEvent);
   gm2fieldUtil::Graph::SetGraphParameters(g1,20,kBlack);

   TGraph *g2 = GetPPTGraph1("TimeStamp","freq",ppEventCor);
   gm2fieldUtil::Graph::SetGraphParameters(g2,20,kRed);

   TMultiGraph *mgPP = new TMultiGraph();
   mgPP->Add(g1,"lp");
   mgPP->Add(g2,"lp");

   TLegend *LP = new TLegend(0.6,0.6,0.8,0.8);
   LP->AddEntry(g1,"Raw","p");
   LP->AddEntry(g2,"FXPR Drift Corrected","p");

   double xMin = tStart/1E+9; 
   double xMax = tStop/1E+9; 

   TCanvas *c2 = new TCanvas("c2","PP & FXPR Data",1200,600);
   c2->Divide(1,2);

   c2->cd(1);
   mgPP->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mgPP,"PP Data","","Frequency (Hz)");
   gm2fieldUtil::Graph::UseTimeDisplay(mgPP);
   gm2fieldUtil::Graph::SetGraphLabelSizes(mgPP,0.05,0.06); 
   mgPP->GetXaxis()->SetLimits(xMin,xMax); 
   mgPP->Draw("a");
   LP->Draw("same");
   c2->Update();

   c2->cd(2);
   gFPAVG->Draw("alp");
   gm2fieldUtil::Graph::SetGraphLabels(gFPAVG,"Fixed Probe Average","","Frequency (Hz)");
   gm2fieldUtil::Graph::UseTimeDisplay(gFPAVG);
   gm2fieldUtil::Graph::SetGraphLabelSizes(gFPAVG,0.05,0.06);
   gFPAVG->GetXaxis()->SetLimits(xMin,xMax); 
   gFPAVG->Draw("alp");
   if(useFXPRFit) fxprFit->Draw("same"); 
   c2->Update();

   TString plotName = Form("%s/pp-shimmed_drift-cor_%s.png",plotDir,date.c_str());

   c2->cd(); 
   c2->Print(plotName);

   return 0;
}

