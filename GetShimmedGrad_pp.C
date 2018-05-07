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
#include "./include/fixedProbeEvent.h"

#include "./src/InputManager.C"
#include "./src/FXPRFuncs.C"
#include "./src/TRLYFuncs.C"
#include "./src/Consolidate.C"
#include "./src/CustomMath.C"
#include "./src/CustomGraph.C"
#include "./src/CustomImport.C"
#include "./src/CustomExport.C"
#include "./src/CustomAlgorithms.C"
#include "./src/CustomUtilities.C"
#include "./src/DeltaBFuncs.C"

int GetShimmedGrad_pp(){

   std::cout << "------------------------------------" << std::endl;
   std::cout << "GET SHIMMED GRADIENT (USING PP SCAN)" << std::endl;

   int rc=0;
   int method = gm2fieldUtil::Constants::kPhaseDerivative;

   std::string input_path = "./input/json/get-shim-grad-pp.json";
   InputManager *inputMgr = new InputManager(); 
   inputMgr->Load(input_path);
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

   if( Axis.compare("x")==0 ) gradType = "rad"; 
   if( Axis.compare("y")==0 ) gradType = "vert"; 
   if( Axis.compare("z")==0 ) gradType = "azi";

   // char inpath[200];
   // sprintf(inpath,"./input/runlists/%s/pp-scan_shimmed_%s.csv",date.c_str(),date.c_str());
   // std::vector<double> sf;
   // ImportDeltaBFileList_csv(inpath,run,label,sf);

   std::vector<int> run,subRun;
   std::vector<std::string> label;
   inputMgr->GetRunList(run); 
   inputMgr->GetRunLabels(label); 

   const int NRUN = run.size();
   if(NRUN==0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   // find the trolley run
   int trlyRun=0,ppRun=0; 
   for(int i=0;i<NRUN;i++){
      if(label[i].compare("trly")==0){
	 trlyRun = run[i]; 
      }else if(label[i].compare("pp")==0){
	 ppRun = run[i];
      }else{
	 subRun.push_back(run[i]);
      }
   }

   // PP data 
   std::vector<plungingProbeAnaEvent_t> ppData,ppEvent,ppEventCor; 
   std::cout << "Getting run(s) " << std::endl; 
   rc = GetPlungingProbeData(ppRun,method,ppData);
   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   const int NEV = ppData.size(); 
   std::cout << "PP events: " << ppData.size() << std::endl;
   for(int i=0;i<NEV;i++) rc = ModifyPlungingProbeData(kLeastSquaresPhase,ppData[i]);
  
   if(isBlind){
      ApplyBlindingPP(blindValue,ppData);
   }

   // now parse PP data using subrun list 
   rc = FilterPlungingProbeData(subRun,ppData,ppEvent); 
   std::cout << "Filtered PP data to use NMR-DAQ runs: " << std::endl;
   const int NPP = ppEvent.size();
   for(int i=0;i<NPP;i++){
      std::cout << Form("%d: x = %.3lf mm, y = %.3lf mm, z = %.3lf mm",
                        ppEvent[i].run,ppEvent[i].r[0],ppEvent[i].y[0],ppEvent[i].phi[0]) << std::endl;
      if(ppEvent[i].run==312){  // FIXME: dumb hack...
	 for(int j=0;j<25;j++) ppEvent[i].y[j] = 0.; 
      }
   } 

   // fixed probe data
   std::string fxprPath = "./input/probe-lists/fxpr-list.csv";
   std::vector<int> fxprList;
   gm2fieldUtil::Import::ImportData1<int>(fxprPath,"csv",fxprList);

   std::vector<gm2field::fixedProbeFrequency_t> fxprData;
   rc = gm2fieldUtil::RootHelper::GetFPFrequencies(ppRun,fxprData);
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

   rc = CorrectPPForDriftDuringMeasurement(fxprDataAvg,ppEvent,ppEventCor,true);

   // Plunging probe plots 
   TGraphErrors *g = GetPPTGraphErrors2(Axis.c_str(),"freq",ppEvent);
   gm2fieldUtil::Graph::SetGraphParameters(g,20,kBlack);

   TGraphErrors *g_cor = GetPPTGraphErrors2(Axis.c_str(),"freq",ppEventCor);
   gm2fieldUtil::Graph::SetGraphParameters(g_cor,20,kRed);

   // Fixed probe plot 
   TGraph *gFXPR = GetTGraphNew(fxprDataAvg);
   gm2fieldUtil::Graph::SetGraphParameters(gFXPR,20,kBlack);  

   TString xAxisLabel = Form("%s (mm)",Axis.c_str());

   TCanvas *c1 = new TCanvas("c1","PP Data",1200,600);
   c1->Divide(1,2);
  
   c1->cd(1); 
   g->Draw("ap");
   gm2fieldUtil::Graph::SetGraphLabels(g,"Shimmed Field",xAxisLabel,"Frequency (Hz)"); 
   gm2fieldUtil::Graph::SetGraphLabelSizes(g,0.05,0.06); 
   g->Draw("ap");
   g->Fit(fitFunc.c_str(),"Q"); 
   c1->Update(); 

   c1->cd(2); 
   g_cor->Draw("ap");
   gm2fieldUtil::Graph::SetGraphLabels(g_cor,"Drift Corrected (FXPR)",xAxisLabel,"Frequency (Hz)"); 
   gm2fieldUtil::Graph::SetGraphLabelSizes(g_cor,0.05,0.06); 
   g_cor->Draw("ap");
   g_cor->Fit(fitFunc.c_str(),"Q"); 
   c1->Update(); 

   TString plotPath = Form("%s/pp-shimmed-scan-%s_run-%d.png",plotDir,Axis.c_str(),ppRun); 
   c1->cd();
   c1->Print(plotPath); 

   TCanvas *c2 = new TCanvas("c2","FXPR Data",1200,600);

   c2->cd(); 
   gFXPR->Draw("ap");
   gm2fieldUtil::Graph::SetGraphLabels(gFXPR,"Fixed Probe Average","","Frequency (Hz)"); 
   gm2fieldUtil::Graph::UseTimeDisplay(gFXPR);  
   gm2fieldUtil::Graph::SetGraphLabelSizes(gFXPR,0.04,0.04); 
   gFXPR->Draw("ap");
   c2->Update(); 

   plotPath = Form("%s/pp-shimmed-scan_fxpr-avg_run-%d.png",plotDir,ppRun); 
   c2->cd();
   c2->Print(plotPath);

   // get fit info 
   TF1 *myFit     = g->GetFunction(fitFunc.c_str()); 
   TF1 *myFit_cor = g_cor->GetFunction(fitFunc.c_str());

   const int NPAR = myFit->GetNpar();
   double par[NPAR]   ,parErr[NPAR]; 
   double parCor[NPAR],parCorErr[NPAR]; 
   for(int i=0;i<NPAR;i++){
      par[i]       = myFit->GetParameter(i); 
      parErr[i]    = myFit->GetParError(i); 
      parCor[i]    = myFit_cor->GetParameter(i); 
      parCorErr[i] = myFit_cor->GetParError(i); 
   }

   // FIXME: finish implementing trolley probe coordinates 
   double pos[2] = {0,0}; 
   rc = GetTrolleyProbePosition(probeNumber,pos);

   // FIXME: What to do for axis = 2 (z)?
   double trly_pr_pos=0;
   if(axis!=2) trly_pr_pos = pos[axis];  

   // get gradient based on fit type and which trolley probe we're looking at  
   double dBdr=0,dBdr_err=0,dBdr_cor=0,dBdr_cor_err=0;
   double PR[3],ER[3];
   int slopeIndex=0; 
   if( fitFunc.compare("pol1")==0 ){
      dBdr         = par[1];  
      dBdr_err     = parErr[1]; 
      dBdr_cor     = parCor[1];  
      dBdr_cor_err = parCorErr[1]; 
   }else if( fitFunc.compare("pol2")==0 ){
      dBdr         = par[1]    + 2.*par[2]*trly_pr_pos; 
      dBdr_err     = TMath::Sqrt( parErr[1]*parErr[1] + parErr[2]*parErr[2] );  
      dBdr_cor     = parCor[1] + 2.*parCor[2]*trly_pr_pos;  
      dBdr_cor_err = TMath::Sqrt( parCorErr[1]*parCorErr[1] + parCorErr[2]*parErr[2] );  
   }else if( fitFunc.compare("pol3")==0 ){
      dBdr         = par[1]    + 2.*par[2]*trly_pr_pos    + 3.*par[3]*TMath::Power(trly_pr_pos,2.);  
      dBdr_err     = TMath::Sqrt( parErr[1]*parErr[1] + parErr[2]*parErr[2] + parErr[3]*parErr[3] );  
      dBdr_cor     = parCor[1] + 2.*parCor[2]*trly_pr_pos + 3.*parCor[3]*TMath::Power(trly_pr_pos,2.);  
      dBdr_cor_err = TMath::Sqrt( parCorErr[1]*parCorErr[1] + parCorErr[2]*parCorErr[2] + parCorErr[3]*parCorErr[3] );  
   }
 
   PR[0] = dBdr; 
   PR[1] = dBdr_cor;  
   PR[2] = 0;

   ER[0] = dBdr_err; 
   ER[1] = dBdr_cor_err; 
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
   char outpath[200],tag[10];
   sprintf(outpath,"%s/%s-grad_pr-%02d_%s.csv" ,outdir,gradType.Data(),probeNumber,date.c_str());
   sprintf(tag    ,"%s-grad",gradType.Data()); 
   std::string strTag = tag;
   PrintToFile(outpath,strTag,PR,ER,drift,drift_err);

   delete inputMgr; 
   
   return 0;
}

