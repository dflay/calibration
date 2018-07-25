// Compute DeltaB for TRLY data sets 
// Uses an ABA approach 
// For use with PRODUCTION calibration data  

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

int DeltaB_trly_prod(std::string configFile){

   std::cout << "----------------------------------" << std::endl;
   std::cout << "TRLY DELTA B CALCULATION (ABA Method)" << std::endl;

   int rc=0;
   int method = gm2fieldUtil::Constants::kPhaseDerivative;

   InputManager *inputMgr = new InputManager();
   inputMgr->Load(configFile);
   inputMgr->Print();

   std::string date = inputMgr->GetAnalysisDate();
   bool isBlind     = inputMgr->IsBlind();
   int probeNumber  = inputMgr->GetTrolleyProbe();
   int axis         = inputMgr->GetAxis();

   date_t theDate; 
   GetDate(theDate);

   blind_t blind; 
   ImportBlinding(blind);
   double blindValue = blind.value_tr; 

   std::string gradName;
   if(axis==0) gradName = "rad"; 
   if(axis==1) gradName = "vert"; 
   if(axis==2) gradName = "azi"; 

   std::vector<int> run;
   std::vector<std::string> label;
   inputMgr->GetRunList(run);
   inputMgr->GetRunLabels(label);

   const int NRUN = run.size();
   if(NRUN==0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   // determine the correct ordering of the SCC on/off cycles 
   std::vector<gm2field::surfaceCoils_t> sccData;
   for(int i=0;i<NRUN;i++) rc = gm2fieldUtil::RootHelper::GetSCCData(run[i],sccData);
   if(rc!=0) return 1; 

   int coilSet  = 0;
   if(axis==2) coilSet = -1;
 
   double thr   = 10E-3;         // in A 
   double delta = 2.*60. + 40.;  // 2 min 40 sec   
   std::vector<double> sccOff,sccOn;
   rc = FindTransitionTimes(coilSet,thr,delta,sccData,sccOff,sccOn);

   bool sccStartOn = false; 
   if(rc==1) sccStartOn = true;

   if(sccStartOn){
      std::cout << "SCC was ON to start the sequence" << std::endl;
   }else{
      std::cout << "SCC was OFF to start the sequence" << std::endl;
   } 

   // TRLY data
   std::vector<trolleyAnaEvent_t> trlyData;
   for(int i=0;i<NRUN;i++){
      rc = GetTrolleyData("",run[i],method,trlyData);
      if(rc!=0){
	 std::cout << "No data!" << std::endl;
	 return 1;
      }
   }

   if(isBlind) ApplyBlindingTRLY(blindValue,trlyData);

   // get the mean field with SCC off and on 
   int nev = 20;  // keep 20 events in the analysis
   std::vector<double> bareTime,bare,bareErr,sccTime,scc,sccErr;  
   rc = GetTRLYStats_sccToggle(probeNumber-1,nev,sccOff,trlyData,bareTime,bare,bareErr); 
   rc = GetTRLYStats_sccToggle(probeNumber-1,nev,sccOn ,trlyData,sccTime ,scc ,sccErr ); 

   // do delta B calcs
   // raw difference 
   std::vector<double> diff,diffErr;
   rc = GetDifference(scc,sccErr,bare,bareErr,diff,diffErr);  
   double mean  = gm2fieldUtil::Math::GetMean<double>(diff); 
   double stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(diff); 

   // ABA difference (bare first) 
   std::vector<double> diff_aba,diffErr_aba;
   if(sccStartOn){
      rc = GetDifference_ABA_sccFirst(scc,sccErr,bare,bareErr,diff_aba,diffErr_aba);  
   }else{
      rc = GetDifference_ABA(scc,sccErr,bare,bareErr,diff_aba,diffErr_aba);  
   }

   double mean_aba  = gm2fieldUtil::Math::GetMean<double>(diff_aba); 
   double stdev_aba = gm2fieldUtil::Math::GetStandardDeviation<double>(diff_aba); 
 
   double dB[3]        = {0,0,0}; 
   double dB_err[3]    = {0,0,0};
   double drift[3]     = {0,0,0};  
   double drift_err[3] = {0,0,0};  

   dB[0]     = mean;
   dB_err[0] = stdev;

   dB[1]     = mean_aba;
   dB_err[1] = stdev_aba;  

   const int ND = diff.size();
   std::vector<double> trial; 
   for(int i=0;i<ND;i++) trial.push_back(i+1);
   
   const int NDA = diff_aba.size();
   std::vector<double> trial_aba;
   for(int i=0;i<NDA;i++) trial_aba.push_back(i+1);
   
   // Plots
 
   TGraphErrors *gTRLY_bare  = gm2fieldUtil::Graph::GetTGraphErrors(bareTime ,bare      ,bareErr      );
   TGraphErrors *gTRLY_scc   = gm2fieldUtil::Graph::GetTGraphErrors(sccTime  ,scc       ,sccErr       );
   TGraphErrors *gDiff       = gm2fieldUtil::Graph::GetTGraphErrors(trial    ,diff      ,diffErr      );  
   TGraphErrors *gDiff_aba   = gm2fieldUtil::Graph::GetTGraphErrors(trial_aba,diff_aba  ,diffErr_aba  );  

   gm2fieldUtil::Graph::SetGraphParameters(gTRLY_bare ,20,kBlack  );
   gm2fieldUtil::Graph::SetGraphParameters(gTRLY_scc  ,20,kRed    );
   gm2fieldUtil::Graph::SetGraphParameters(gDiff      ,20,kBlue   );
   gm2fieldUtil::Graph::SetGraphParameters(gDiff_aba  ,20,kGreen+2);

   TMultiGraph *mgTRLY = new TMultiGraph();
   mgTRLY->Add(gTRLY_bare,"lp");
   mgTRLY->Add(gTRLY_scc ,"lp");

   TMultiGraph *mgDiff = new TMultiGraph();
   mgDiff->Add(gDiff      ,"lp");
   mgDiff->Add(gDiff_aba  ,"lp");

   TLegend *L = new TLegend(0.6,0.6,0.8,0.8); 
   L->AddEntry(gTRLY_bare   ,"Bare","p"); 
   L->AddEntry(gTRLY_scc    ,"SCC" ,"p");
   
   TLegend *LD = new TLegend(0.6,0.6,0.8,0.8); 
   LD->AddEntry(gDiff      ,"Raw"  ,"p");  
   LD->AddEntry(gDiff_aba  ,"ABA"  ,"p");  

   std::cout << Form("====================== RESULTS FOR PROBE %02d ======================",probeNumber) << std::endl;
   std::cout << "Raw results: " << std::endl;
   std::cout << Form("dB (%s): %.3lf +/- %.3lf Hz (%.3lf +/- %.3lf ppb)",gradName.c_str(),
                     dB[0],dB_err[0],dB[0]/0.06179,dB_err[0]/0.06179) << std::endl;
   std::cout << "-----------------------------------------------------" << std::endl; 
   std::cout << "Drift corrected [ABA]: " << std::endl;
   std::cout << Form("dB (%s): %.3lf +/- %.3lf Hz (%.3lf +/- %.3lf ppb)",gradName.c_str(),
                     dB[1],dB_err[1],dB[1]/0.06179,dB_err[1]/0.06179) << std::endl;

   char outpath[200],outdir[200],datedir[50];
   sprintf(datedir,"%02d-%02d-%02d",theDate.month,theDate.day,theDate.year-2000); 

   if(isBlind)  sprintf(outdir,"./output/blinded/%s"  ,datedir); 
   if(!isBlind) sprintf(outdir,"./output/unblinded/%s",datedir); 

   rc = MakeDirectory(outdir); 
   sprintf(outpath,"%s/dB-trly_final-location_%s-grad_pr-%02d_%s.csv"  ,outdir,gradName.c_str(),probeNumber,date.c_str());
   rc = PrintToFile(outpath,gradName,dB,dB_err,drift,drift_err); 

   char plotDir[200];
   sprintf(plotDir,"./plots/%s",datedir); 
   rc = MakeDirectory(plotDir); 

   // draw some plots 
   TCanvas *c1 = new TCanvas("c1","Data",1200,600);
   c1->Divide(1,2); 

   c1->cd(1);
   mgTRLY->Draw("alp");
   gm2fieldUtil::Graph::SetGraphLabels(mgTRLY,"PP Data","","Frequency (Hz)"); 
   gm2fieldUtil::Graph::UseTimeDisplay(mgTRLY);
   // mgTRLY_bare->GetXaxis()->SetLimits(xMin_bare,xMax_bare); 
   mgTRLY->Draw("alp");
   L->Draw("same"); 
   c1->Update();

   c1->cd(2);
   mgDiff->Draw("alp");
   gm2fieldUtil::Graph::SetGraphLabels(mgDiff,"SCC-Bare","Trial","Frequency Difference (Hz)"); 
   // mgTRLY_bare->GetXaxis()->SetLimits(xMin_bare,xMax_bare); 
   mgDiff->Draw("alp");
   LD->Draw("same"); 
   c1->Update();

   c1->cd();
   TString plotPath = Form("%s/trly_dB_%s-grad_run-%d_pr-%02d.png",plotDir,gradName.c_str(),run[0],probeNumber); 
   c1->Print(plotPath);
   delete c1;  

   delete inputMgr; 

   return 0;
}