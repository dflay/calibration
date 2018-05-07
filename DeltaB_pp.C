// Compute DeltaB for PP data sets  

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
#include "./src/CustomUtilities.C"
#include "./src/DeltaBFuncs.C"

int DeltaB_pp(std::string date,int type,int finalPos,bool isBlind,bool isFullAnalysis,bool useFXPRFit){

   int rc=0;
   int method = gm2fieldUtil::Constants::kPhaseDerivative;

   date_t theDate; 
   GetDate(theDate);

   blind_t blind; 
   ImportBlinding(blind);
   double blindValue = blind.value_pp; 

   std::cout << "----------------------------------" << std::endl;
   std::cout << "PP DELTA B CALCULATION" << std::endl;

   std::string gradName;
   if(type==0) gradName = "rad"; 
   if(type==1) gradName = "vert"; 
   if(type==2) gradName = "azi"; 

   char tag[100];
   if(finalPos==1){
      sprintf(tag,"_final-location_"); 
   }else if(finalPos==0){
      sprintf(tag,"_"); 
   }else{
      std::cout << "Invalid type! finalPos = " << finalPos <<std::endl;
      return 1;
   }
   
   char inpath[200];   
   sprintf(inpath,"./input/runlists/%s/dB-pp%s%s-grad_%s.csv",date.c_str(),tag,gradName.c_str(),date.c_str());

   std::vector<int> allRuns;
   std::vector<double> sf;
   std::vector<std::string> label;
   rc = ImportDeltaBFileList_csv(inpath,allRuns,label,sf);
   if(rc!=0) return 1;

   const int NRUN = allRuns.size();

   std::vector<int> index,run,driftRun;
   SortRuns(label,allRuns,run,driftRun,index); 

   int bareIndex  = index[0]; 
   int gradIndex  = index[1]; 
   int bare2Index = index[2]; 

   // fixed probe list
   std::string fxprPath = "./input/probe-lists/fxpr-list.csv"; 
   std::vector<int> fxprList; 
   gm2fieldUtil::Import::ImportData1<int>(fxprPath,"csv",fxprList);

   // get graph of fxpr data for drift correction *across* runs 
   std::vector<double> stats;
   TGraph *gDrift = GetDriftTGraph(method,driftRun,fxprList,stats); 

   // now get FXPR data for drift corrections *during* runs
   // all runs
   std::vector<gm2field::fixedProbeFrequency_t> fxprData;
   for(int i=0;i<NRUN;i++) rc = gm2fieldUtil::RootHelper::GetFPFrequencies(allRuns[i],fxprData);
   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }
 
   // for the bare run
   std::cout << "Getting FXPR data for run " << run[bareIndex] << std::endl; 
   std::vector<gm2field::fixedProbeFrequency_t> fxprData_bare;
   rc = gm2fieldUtil::RootHelper::GetFPFrequencies(run[bareIndex],fxprData_bare);
   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   // for the grad run
   std::cout << "Getting FXPR data for run " << run[gradIndex] << std::endl; 
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

   // get an error estimate by subtracting off the fit, taking the stdev

   double mean_p2p=0; 
   double stdev_p2p_bare=0,stdev_p2p_grad=0;
   double STDEV_P2P_BARE=0,STDEV_P2P_GRAD=0;
   TGraph *gFPAVG_bare_trm = RemoveTrend(gFPAVG_bare,fxprFit_bare,mean_p2p,STDEV_P2P_BARE);
   TGraph *gFPAVG_grad_trm = RemoveTrend(gFPAVG_grad,fxprFit_grad,mean_p2p,STDEV_P2P_GRAD);

   if(useFXPRFit){
      stdev_p2p_bare = STDEV_P2P_BARE;
      stdev_p2p_grad = STDEV_P2P_GRAD;
   } 

   TString gLabel;

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
   // std::string trlyPath = "./input/probe-lists/trly-list.csv"; 
   // std::vector<int> trlyList;
   // gm2fieldUtil::Import::ImportData1<int>(trlyPath,"csv",trlyList);

   // std::vector<trolleyAnaEvent_t> trlyData;
   // for(int i=0;i<NRUN;i++) rc = GetTrolleyData(run[i],method,trlyData); 

   // TMultiGraph *mg2 = new TMultiGraph();
   // TLegend *L2 = new TLegend(0.6,0.6,0.8,0.8); 
   // const int NTR = trlyList.size();
   // TGraph **gTRLY = new TGraph*[NTR];
   // for(int i=0;i<NTR;i++){
   //    gTRLY[i] = GetTRLYTGraph(trlyList[i],"GpsTimeStamp","freq",trlyData);
   //    gm2fieldUtil::Graph::SetGraphParameters(gTRLY[i],20,kWhite+i+1);
   //    mg2->Add(gTRLY[i],"p");
   //    gLabel = Form("trly %03d",trlyList[i]); 
   //    L2->AddEntry(gTRLY[i],gLabel.Data(),"p"); 
   // }

   // c1->cd(2);
   // mg2->Draw("a");
   // gm2fieldUtil::Graph::SetGraphLabels(mg2,"Trolley Probe Data","","Frequency (Hz)");
   // gm2fieldUtil::Graph::UseTimeDisplay(mg2);  
   // mg2->Draw("a");
   // L2->Draw("same");  
   // c1->Update();

   // PP data 
   const int N3 = run.size(); 
   std::vector<plungingProbeAnaEvent_t> ppEvent,ppEvent_tr,ppEventCor,ppEventCor_tr; 
   for(int i=0;i<N3;i++){
      std::cout << "Getting run " << run[i] << std::endl; 
      rc = GetPlungingProbeData(run[i],method,ppEvent);
      if(rc!=0){
	 std::cout << "No data!" << std::endl;
	 return 1;
      }
   }

   const int NPP = ppEvent.size();
   for(int i=0;i<NPP;i++) rc = ModifyPlungingProbeData(kLeastSquaresPhase,ppEvent[i]);

   if(isBlind) ApplyBlindingPP(blindValue,ppEvent); 
  
   CopyPlungingProbe(ppEvent,ppEvent_tr);

   plungingProbeAnaEvent_t ppEventCor_bare,ppEventCor_grad,ppEventCor_new; 

   if(useFXPRFit){
      rc = CorrectPPForDriftDuringMeasurement(method,fxprFit_bare,ppEvent[bareIndex],ppEventCor_bare);
      rc = CorrectPPForDriftDuringMeasurement(method,fxprFit_grad,ppEvent[gradIndex],ppEventCor_grad);
      ppEventCor.push_back(ppEventCor_bare); 
      ppEventCor.push_back(ppEventCor_grad); 
   }else{ 
      rc = CorrectPPForDriftDuringMeasurement(fxprDataAvg_bare,ppEvent[bareIndex],ppEventCor_bare);
      rc = CorrectPPForDriftDuringMeasurement(fxprDataAvg_grad,ppEvent[gradIndex],ppEventCor_grad);
      ppEventCor.push_back(ppEventCor_bare); 
      ppEventCor.push_back(ppEventCor_grad);
      // for(int i=0;i<NN;i++){
      //    rc = CorrectPPForDriftDuringMeasurement(method,fxprList,fxprData,ppEvent[i],ppEventCor_new);
      //    ppEventCor.push_back(ppEventCor_new); 
      // }
   }
   // rc = CorrectPPForDriftDuringMeasurement(method,trlyList,trlyData,ppEvent_tr,ppEventCor_tr);

   // now do delta B calcs
   double theDrift=0,theDriftErr=0;
   double dB[3]        = {0,0,0}; 
   double dB_err[3]    = {0,0,0};
   double drift[3]     = {0,0,0};  
   double drift_err[3] = {0,0,0};  
   // raw 
   CalculatePPDeltaB(false,method,gDrift,ppEvent[bareIndex],ppEvent[gradIndex],dB[0],dB_err[0],drift[0],drift_err[0]);
   // drift corrected [fxpr] 
   CalculatePPDeltaB(true ,method,gDrift,ppEventCor[bareIndex],ppEventCor[gradIndex],dB[1],dB_err[1],drift[1],drift_err[1]);
   // drift corrected [trly] 
   // CalculatePPDeltaB(true,method,trlyList,trlyData,ppEventCor_tr[bareIndex],ppEventCor_tr[gradIndex],dB[2],dB_err[2],drift[2],drift_err[2]);

   // Plunging probe plots 
   TGraph *gPP_bare = GetPPTGraph3("TimeStamp","freq",ppEvent[bareIndex]);
   gm2fieldUtil::Graph::SetGraphParameters(gPP_bare,20,kBlack);

   TGraph *gPP_bare_fxpr = GetPPTGraph3("TimeStamp","freq",ppEventCor[bareIndex]);
   gm2fieldUtil::Graph::SetGraphParameters(gPP_bare_fxpr,20,kRed);

   TGraph *gPP_grad = GetPPTGraph3("TimeStamp","freq",ppEvent[gradIndex]);
   gm2fieldUtil::Graph::SetGraphParameters(gPP_grad,20,kBlack);

   TGraph *gPP_grad_fxpr = GetPPTGraph3("TimeStamp","freq",ppEventCor[gradIndex]);
   gm2fieldUtil::Graph::SetGraphParameters(gPP_grad_fxpr,20,kRed);

   TMultiGraph *mgPP_bare = new TMultiGraph();
   mgPP_bare->Add(gPP_bare     ,"lp");
   mgPP_bare->Add(gPP_bare_fxpr,"lp");

   TMultiGraph *mgPP_grad = new TMultiGraph();
   mgPP_grad->Add(gPP_grad     ,"lp");
   mgPP_grad->Add(gPP_grad_fxpr,"lp");

   TLegend *L = new TLegend(0.6,0.6,0.8,0.8); 
   L->AddEntry(gPP_bare,"Raw"           ,"p"); 
   L->AddEntry(gPP_bare_fxpr,"FXPR drift cor","p"); 

   // when trolley drift correction won't work... 
   dB[2]     = dB[0];
   dB_err[2] = dB_err[0];

   // add the stdev from the anchor points as the drift error
   // also add in stdev of P2P stdev after subtracting off the fit for both bare and grad runs
   double arg_sq=0; 
   for(int i=1;i<3;i++){
      arg_sq       = stats[1]*stats[1] + stats[3]*stats[3] + stdev_p2p_bare*stdev_p2p_bare + stdev_p2p_grad*stdev_p2p_grad; 
      drift_err[i] = TMath::Sqrt(arg_sq);
   }

   std::cout << "====================== RESULTS ======================" << std::endl;
   std::cout << "Raw results: " << std::endl;
   std::cout << Form("dB (%s): %.3lf +/- %.3lf Hz (%.3lf +/- %.3lf ppb)",label[gradIndex].c_str(),
                     dB[0],dB_err[0],dB[0]/0.06179,dB_err[0]/0.06179) << std::endl;
   std::cout << "-----------------------------------------------------" << std::endl; 
   std::cout << "Drift corrected [FXPR]: " << std::endl;
   std::cout << Form("dB (%s): %.3lf +/- %.3lf +/- %.3lf Hz (%.3lf +/- %.3lf +/- %.3lf ppb)",label[gradIndex].c_str(),
                     dB[1],dB_err[1],drift_err[1],dB[1]/0.06179,dB_err[1]/0.06179,drift_err[1]/0.06179) << std::endl;
   std::cout << "-----------------------------------------------------" << std::endl; 
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
   sprintf(outpath,"%s/dB-pp%s%s-grad_%s.csv"  ,outdir,tag,gradName.c_str(),date.c_str()); 
   rc = PrintToFile(outpath,label[gradIndex],dB,dB_err,drift,drift_err); 

   char plotDir[200];
   sprintf(plotDir,"./plots/%02d-%02d-%02d",theDate.month,theDate.day,theDate.year-2000); 
   rc = MakeDirectory(plotDir); 


   // draw some plots 
   TCanvas *c1 = new TCanvas("c1","Bare Run",1200,600);
   c1->Divide(1,2); 

   c1->cd(1);
   mgPP_bare->Draw("alp");
   gm2fieldUtil::Graph::SetGraphLabels(mgPP_bare,Form("PP Run %d",run[bareIndex]),"","Frequency (Hz)"); 
   gm2fieldUtil::Graph::UseTimeDisplay(mgPP_bare);
   mgPP_bare->GetXaxis()->SetLimits(xMin_bare,xMax_bare); 
   mgPP_bare->Draw("alp");
   L->Draw("same"); 
   c1->Update();

   c1->cd(2);
   gFPAVG_bare->Draw("alp");
   gm2fieldUtil::Graph::SetGraphLabels(gFPAVG_bare,Form("FXPR Run %d",run[bareIndex]),"","Frequency (Hz)"); 
   gm2fieldUtil::Graph::UseTimeDisplay(gFPAVG_bare); 
   gFPAVG_bare->GetXaxis()->SetLimits(xMin_bare,xMax_bare); 
   gFPAVG_bare->Draw("alp");
   if(useFXPRFit) fxprFit_bare->Draw("same"); 
   c1->Update(); 

   c1->cd();
   TString plotPath = Form("%s/pp_bare-run-%d.png",plotDir,run[bareIndex]); 
   c1->Print(plotPath);
   delete c1;  

   TCanvas *c2 = new TCanvas("c2","Grad Run",1200,600);
   c2->Divide(1,2);
 
   c2->cd(1);
   mgPP_grad->Draw("alp");
   gm2fieldUtil::Graph::SetGraphLabels(mgPP_grad,Form("PP Run %d",run[gradIndex]),"","Frequency (Hz)"); 
   gm2fieldUtil::Graph::UseTimeDisplay(mgPP_grad); 
   mgPP_grad->GetXaxis()->SetLimits(xMin_grad,xMax_grad); 
   mgPP_grad->Draw("alp");
   L->Draw("same"); 
   c2->Update();

   c2->cd(2);
   gFPAVG_grad->Draw("alp");
   gm2fieldUtil::Graph::SetGraphLabels(gFPAVG_grad,Form("FXPR Run %d",run[gradIndex]),"","Frequency (Hz)"); 
   gm2fieldUtil::Graph::UseTimeDisplay(gFPAVG_grad); 
   gFPAVG_grad->GetXaxis()->SetLimits(xMin_grad,xMax_grad); 
   gFPAVG_grad->Draw("alp");
   if(useFXPRFit) fxprFit_grad->Draw("same"); 
   c2->Update();

   c2->cd();
  
   plotPath = Form("%s/pp_grad-run-%d.png",plotDir,run[gradIndex]); 
   c2->Print(plotPath); 
   delete c2;  

   return 0;
}

