// Determine the gradient in the field using the PP (scan) data     

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

int GetImposedGradient_pp(std::string configFile){

   std::cout << "------------------------------------" << std::endl;
   std::cout << "GET IMPOSED GRADIENT (USING PP SCAN)" << std::endl;

   int rc=0;
   int method = gm2fieldUtil::Constants::kPhaseDerivative;

   InputManager *inputMgr = new InputManager();
   inputMgr->Load(configFile);
   inputMgr->Print();

   std::string date    = inputMgr->GetAnalysisDate();
   std::string fitFunc = inputMgr->GetFitFunction();
   bool isBlind        = inputMgr->IsBlind();
   bool isFullAnalysis = inputMgr->IsFullAnalysis();
   int probeNumber     = inputMgr->GetTrolleyProbe();
   int axis            = inputMgr->GetAxis();

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
   double blindValue = blind.value_pp;

   TString gradType;
   std::string Axis; 
   if(axis==0) Axis = "x"; 
   if(axis==1) Axis = "y"; 
   if(axis==2) Axis = "z"; 

   if( Axis.compare("x")==0 ) gradType = "norm-quad"; 
   if( Axis.compare("y")==0 ) gradType = "skew-quad"; 
   if( Axis.compare("z")==0 ) gradType = "azi";

   std::cout << Form("axis: %s (%s grad)",Axis.c_str(),gradType.Data()) << std::endl;  

   char inpath[200];
   sprintf(inpath,"./input/runlists/%s/pp-scan_%s-grad_%s.csv",date.c_str(),gradType.Data(),date.c_str());

   // std::vector<int> allRuns,run,bRun,gRun;
   // std::vector<double> sf;
   // std::vector<std::string> label;
   // ImportDeltaBFileList_csv(inpath,allRuns,label,sf);

   // const int NRUN = allRuns.size();
   // if(NRUN==0){
   //    std::cout << "No data!" << std::endl;
   //    return 1;
   // }

   std::vector<int> allRuns,run,bRun,gRun;
   std::vector<std::string> label;
   inputMgr->GetRunList(allRuns);
   inputMgr->GetRunLabels(label);

   const int NRUN = allRuns.size();
   if(NRUN==0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   for(int i=0;i<NRUN;i++){
      if(label[i].compare("skew-quad")==0) gRun.push_back(allRuns[i]);
      if(label[i].compare("norm-quad")==0) gRun.push_back(allRuns[i]);
      if(label[i].compare("azi")==0)       gRun.push_back(allRuns[i]);
      if(label[i].compare("bare")==0)      bRun.push_back(allRuns[i]);
      if(label[i].compare("bare-2")!=0){
         // save all runs that aren't the last bare run
         run.push_back(allRuns[i]);
      }
   }

   std::cout << "Gradient runs: " << std::endl;
   const int NG = gRun.size();
   for(int i=0;i<NG;i++) std::cout << gRun[i] << std::endl;

   const int NB = bRun.size();
   std::cout << "Bare field runs: " << std::endl;
   for(int i=0;i<NB;i++) std::cout << bRun[i] << std::endl;

   // PP data 
   std::vector<plungingProbeAnaEvent_t> ppEventBare,ppEventBareCor; 
   std::cout << "Getting bare run(s) " << std::endl; 
   for(int i=0;i<NB;i++) rc = GetPlungingProbeData(bRun[i],method,ppEventBare);
   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   std::vector<plungingProbeAnaEvent_t> ppEventGrad,ppEventGradCor; 
   std::cout << "Getting grad run(s) " << std::endl; 
   for(int i=0;i<NG;i++) rc = GetPlungingProbeData(gRun[i],method,ppEventGrad);
   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   const int NPPB = ppEventBare.size();
   for(int i=0;i<NPPB;i++) rc = ModifyPlungingProbeData(kLeastSquaresPhase,ppEventBare[i]);
   
   const int NPPG = ppEventGrad.size();
   for(int i=0;i<NPPG;i++) rc = ModifyPlungingProbeData(kLeastSquaresPhase,ppEventGrad[i]);

   if(isBlind){
      ApplyBlindingPP(blindValue,ppEventBare);
      ApplyBlindingPP(blindValue,ppEventGrad);
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

   const int NN = fxprData.size();
   const int NFP = fxprList.size();
   unsigned long long t0     = 0;
   unsigned long long tStart = fxprData[0].GpsTimeStamp[fxprList[0]];
   unsigned long long tStop  = fxprData[NN-1].GpsTimeStamp[fxprList[NFP-1]];
   unsigned long long tStep  = 1E+9;

   std::vector<fixedProbeEvent_t> fxprDataAvg;
   GetAverageFXPRVectorsNew(method,t0,tStart,tStop,tStep,fxprList,fxprData,fxprDataAvg);

   rc = CorrectPPForDriftDuringMeasurement(fxprDataAvg,ppEventBare,ppEventBareCor,true);
   rc = CorrectPPForDriftDuringMeasurement(fxprDataAvg,ppEventGrad,ppEventGradCor,true);

   // std::vector<double> fb,fa;
   // const int NB = ppEventBare[0].numTraces;  
   // for(int i=0;i<NB;i++) fb.push_back( ppEvent[0].freq[i] ); 
   // double mean_before  = gm2fieldUtil::Math::GetMean<double>(fb); 
   // double stdev_before = gm2fieldUtil::Math::GetStandardDeviation<double>(fb); 

   // const int NA = ppEventBareCor[0].numTraces;  
   // for(int i=0;i<NA;i++) fa.push_back( ppEventCor[0].freq[i] ); 
   // double mean_after  = gm2fieldUtil::Math::GetMean<double>(fa); 
   // double stdev_after = gm2fieldUtil::Math::GetStandardDeviation<double>(fa); 
  
   // std::cout << Form("mean (before) = %.3lf +/- %.3lf Hz",mean_before,stdev_before) << std::endl; 
   // std::cout << Form("mean (after)  = %.3lf +/- %.3lf Hz",mean_after ,stdev_after)  << std::endl;

   // Plunging probe plots 
   TGraphErrors *gb = GetPPTGraphErrors2(Axis.c_str(),"freq",ppEventBare);
   gm2fieldUtil::Graph::SetGraphParameters(gb,20,kBlack);

   TGraphErrors *gb_cor = GetPPTGraphErrors2(Axis.c_str(),"freq",ppEventBareCor);
   gm2fieldUtil::Graph::SetGraphParameters(gb_cor,20,kRed);

   TGraphErrors *gg = GetPPTGraphErrors2(Axis.c_str(),"freq",ppEventGrad);
   gm2fieldUtil::Graph::SetGraphParameters(gg,20,kBlack);

   TGraphErrors *gg_cor = GetPPTGraphErrors2(Axis.c_str(),"freq",ppEventGradCor);
   gm2fieldUtil::Graph::SetGraphParameters(gg_cor,20,kRed);

   TGraph *gDiff     = GetDiffPlot(gb,gg);
   gm2fieldUtil::Graph::SetGraphParameters(gDiff,21,kBlack); 

   TGraph *gDiff_cor = GetDiffPlot(gb_cor,gg_cor);
   gm2fieldUtil::Graph::SetGraphParameters(gDiff_cor,20,kRed); 

   TMultiGraph *mgb = new TMultiGraph();
   mgb->Add(gb,"lp"); 
   mgb->Add(gb_cor,"lp");

   TMultiGraph *mgg = new TMultiGraph();
   mgg->Add(gg,"lp"); 
   mgg->Add(gg_cor,"lp");

   TMultiGraph *mgd = new TMultiGraph();
   mgd->Add(gDiff,"p"); 
   mgd->Add(gDiff_cor,"p"); 

   // Fixed probe plot 
   TGraph *gFXPR = GetTGraphNew(fxprDataAvg);
   gm2fieldUtil::Graph::SetGraphParameters(gFXPR,20,kBlack);  

   TLegend *L = new TLegend(0.6,0.6,0.8,0.8); 
   L->AddEntry(gb,"Before Drift Correction","p");  
   L->AddEntry(gb_cor,"After Drift Correction","p");  

   TMultiGraph *mg2 = new TMultiGraph();
   mg2->Add(gFXPR,"lp"); 

   TString xAxisLabel = Form("%c (mm)",axis);

   TCanvas *c1 = new TCanvas("c1","PP & FXPR Data",1200,600);
   c1->Divide(1,3);
  
   c1->cd(1); 
   mgb->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mgb,"Bare Field",xAxisLabel,"Frequency (Hz)"); 
   gm2fieldUtil::Graph::SetGraphLabelSizes(mgb,0.05,0.06); 
   mgb->Draw("a");
   L->Draw("same"); 
   c1->Update(); 

   c1->cd(2); 
   mgg->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mgg,"Grad Applied",xAxisLabel,"Frequency (Hz)"); 
   gm2fieldUtil::Graph::SetGraphLabelSizes(mgg,0.05,0.06); 
   mgg->Draw("a");
   c1->Update(); 

   c1->cd(3); 
   mgd->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mgd,"Gradient",xAxisLabel,"Frequency (Hz)"); 
   gm2fieldUtil::Graph::SetGraphLabelSizes(mgd,0.05,0.06); 
   mgd->Draw("a");
   c1->Update(); 

   TString plotPath = Form("%s/pp-scan_%s-grad_bare-run-%d_grad-run-%d.png",plotDir,gradType.Data(),bRun[0],gRun[0]); 
   c1->cd();
   c1->Print(plotPath); 

   TCanvas *c2 = new TCanvas("c2","PP & FXPR Data",1200,600);
   c2->Divide(2,1);
  
   c2->cd(1); 
   gDiff->Draw("ap");
   gm2fieldUtil::Graph::SetGraphLabels(gDiff,"Gradient",xAxisLabel,"Frequency (Hz)"); 
   gm2fieldUtil::Graph::SetGraphLabelSizes(gDiff,0.04,0.04,0.8); 
   gDiff->Draw("ap");
   gDiff->Fit(fitFunc.c_str(),"Q"); 
   c2->Update(); 

   c2->cd(2); 
   gDiff_cor->Draw("ap");
   gm2fieldUtil::Graph::SetGraphLabels(gDiff_cor,"Gradient (FXPR)",xAxisLabel,"Frequency (Hz)"); 
   gm2fieldUtil::Graph::SetGraphLabelSizes(gDiff_cor,0.04,0.04,0.8); 
   gDiff_cor->Draw("ap");
   gDiff_cor->Fit(fitFunc.c_str(),"Q"); 
   c2->Update(); 

   plotPath = Form("%s/pp-scan_%s-grad_bare-run-%d_grad-run-%d_fits.png",plotDir,gradType.Data(),bRun[0],gRun[0]); 
   c2->cd();
   c2->Print(plotPath);

   TCanvas *c3 = new TCanvas("c3","FXPR Data",1200,600);

   c3->cd(2); 
   mg2->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mg2,"Fixed Probe Average","","Frequency (Hz)"); 
   gm2fieldUtil::Graph::UseTimeDisplay(mg2);  
   gm2fieldUtil::Graph::SetGraphLabelSizes(mg2,0.05,0.06); 
   mg2->Draw("a");
   c3->Update(); 

   plotPath = Form("%s/pp-scan_%s-grad_bare-run-%d_grad-run-%d_fxpr-avg.png",plotDir,gradType.Data(),bRun[0],gRun[0]); 
   c3->cd();
   c3->Print(plotPath);

   TF1 *myFit     = gDiff->GetFunction(fitFunc.c_str()); 
   TF1 *myFit_cor = gDiff_cor->GetFunction(fitFunc.c_str());

   // const int npar1 = myFit->GetNpar();  
   // const int npar2 = myFit_cor->GetNpar(); 

   // double par1[npar1],parErr1[npar1];  
   // double par2[npar2],parErr2[npar2];  
   // 
   // for(int i=0;i<npar1;i++){
   //    par1[i]    = myFit->GetParameter(i); 
   //    parErr1[i] = myFit->GetParError(i); 
   // }

   // for(int i=0;i<npar2;i++){
   //    par2[i]    = myFit_cor->GetParameter(i); 
   //    parErr2[i] = myFit_cor->GetParError(i); 
   // }
 
   // std::cout << "FIT RESULTS" << std::endl;
   // std::cout << "Raw" << std::endl;
   // for(int i=0;i<npar1;i++){
   //    std::cout << Form("p[%d] = %.3lf +/- %.3lf",i,par1[i],parErr1[i]) << std::endl;
   // } 
   // std::cout << "Drift Corrected (FXPR)" << std::endl;
   // for(int i=0;i<npar2;i++){
   //    std::cout << Form("p[%d] = %.3lf +/- %.3lf",i,par2[i],parErr1[i]) << std::endl;
   // } 

   // gather fit results

   int slopeIndex=0; 
   if( fitFunc.compare("pol1")==0 ) slopeIndex = 1; 
   if( fitFunc.compare("pol2")==0 ) slopeIndex = 1; 
   if( fitFunc.compare("pol3")==0 ) slopeIndex = 1; 
 
   double PR[3],ER[3];
   PR[0] = myFit->GetParameter(slopeIndex);
   PR[1] = myFit_cor->GetParameter(slopeIndex);
   PR[2] = 0;

   ER[0] = myFit->GetParError(slopeIndex);
   ER[1] = myFit_cor->GetParError(slopeIndex);
   ER[2] = 0;

   double drift[3]     = {0,0,0};
   double drift_err[3] = {0,0,0};

   // assume the difference as the drift error for now 
   drift_err[1] = TMath::Abs(PR[1]-PR[0]); 

   std::cout << "============================ RESULTS ============================" << std::endl;
   std::cout << Form("%s:        %.3lf +/- %.3lf Hz/mm",gradType.Data(),PR[0],ER[0]) << std::endl;
   std::cout << Form("%s [fxpr]: %.3lf +/- %.3lf +/- %.3lf Hz/mm",gradType.Data(),PR[1],ER[1],drift_err[1]) << std::endl;
   std::cout << Form("%s [trly]: %.3lf +/- %.3lf +/- %.3lf Hz/mm",gradType.Data(),PR[2],ER[2],drift_err[2]) << std::endl;

   char outdir[200];
   if(isBlind)  sprintf(outdir,"./output/blinded/%02d-%02d-%02d"  ,theDate.month,theDate.day,theDate.year-2000);
   if(!isBlind) sprintf(outdir,"./output/unblinded/%02d-%02d-%02d",theDate.month,theDate.day,theDate.year-2000);
   if(!isFullAnalysis){
      // we're doing a delta-b calculation during a measurement
      // store in a different location 
      sprintf(outdir,"./output/delta-b/%02d-%02d-%02d",theDate.month,theDate.day,theDate.year-2000);
   }
   rc = MakeDirectory(outdir);

   // print results 
   char tag[200],outpath[200]; 
   sprintf(tag,"%s-grad",gradType.Data());  
   sprintf(outpath,"%s/%s_pr-%02d_%s.csv" ,outdir,tag,probeNumber,date.c_str());

   std::string tagStr = tag; 
   PrintToFile(outpath,tagStr,PR,ER,drift,drift_err);

   delete inputMgr; 
 
   return 0;
}

