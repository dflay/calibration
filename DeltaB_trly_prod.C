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
#include "TLine.h"

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
#include "./src/TRLYFuncs.C"
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
   inputMgr->UseAxis(); 
   inputMgr->Load(configFile);
   inputMgr->Print();

   std::string date   = inputMgr->GetAnalysisDate();
   bool isBlind       = inputMgr->IsBlind();
   bool useTimeWeight = inputMgr->GetTimeWeightStatus();  
   bool loadTimes     = inputMgr->GetSCCTimeStatus(); 
   int probeNumber    = inputMgr->GetTrolleyProbe();
   int axis           = inputMgr->GetAxis();

   date_t theDate; 
   GetDate(theDate);

   char outpath[200],outdir[200];

   if(isBlind)  sprintf(outdir,"./output/blinded/%s"  ,theDate.getDateString().c_str()); 
   if(!isBlind) sprintf(outdir,"./output/unblinded/%s",theDate.getDateString().c_str());
 
   std::string gradName;
   if(axis==0) gradName = "rad"; 
   if(axis==1) gradName = "vert"; 
   if(axis==2) gradName = "azi";

   rc = MakeDirectory(outdir); 
   sprintf(outpath,"%s/dB-trly_final-location_%s-grad_pr-%02d_%s.csv"  ,outdir,gradName.c_str(),probeNumber,date.c_str());

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

   bool loadOnline = false; 
   for(int i=0;i<NRUN;i++){
      if(run[i]==-1){
	 loadOnline = true;
         std::cout << "Will load online Delta-B value" << std::endl;
      }
   }

   int coilSet  = 0;
   double thr   = 10E-3;         // in A 
   double delta = 2.*60. + 40.;  // 2 min 40 sec   // for x, y 
   
   if(axis==2){
      coilSet = -1;
      delta = 1.*60.;
   }

   // determine the correct ordering of the SCC on/off cycles 
   std::vector<gm2field::surfaceCoils_t> sccData;
   for(int i=0;i<NRUN;i++) rc = gm2fieldUtil::RootHelper::GetSCCData(run[i],sccData);
   if(rc!=0) return 1; 
   
   bool sccStartOn = false; 
   std::vector<double> sccOff,sccOn;
  
   std::vector<deltab_t> trly_dB; 

   char inpath_db[200];
   sprintf(inpath_db,"./input/delta-b/trly_xyz_07-18.csv");

   double dB[3]        = {0,0,0}; 
   double dB_err[3]    = {0,0,0};
   double drift[3]     = {0,0,0};  
   double drift_err[3] = {0,0,0};  

   if(loadOnline){
      // use online results since we don't have a run to work with 
      rc = LoadDeltaBData_trlyXYZ(inpath_db,probeNumber,trly_dB);
      dB[0] = trly_dB[axis].dB; 
      dB[1] = trly_dB[axis].dB_fxpr; 
      dB[2] = 0.; 
      dB_err[0] = trly_dB[axis].dB_err; 
      dB_err[1] = trly_dB[axis].dB_fxpr_err; 
      dB_err[2] = 0.; 
      rc = PrintToFile(outpath,gradName,dB,dB_err,drift,drift_err); 
   }else{
      if(loadTimes){
	 // better to use pre-defined transition times   
	 rc = LoadTRLYSCCTimes(probeNumber,sccOff,sccOn);
	 if(sccOn[0]<sccOff[0]) sccStartOn = true;
      }else{
	 rc = FindTransitionTimes(coilSet,thr,delta,sccData,sccOff,sccOn);
	 if(rc<0){
	    std::cout << "No SCC transitions!" << std::endl;
	    return 1;
	 }
	 if(rc==1) sccStartOn = true;
      }
   }

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

   std::vector<double> X; 
   int NEV = trlyData.size();
   for(int i=0;i<NEV;i++) X.push_back( trlyData[i].freq[probeNumber-1] );
   double MEAN  = gm2fieldUtil::Math::GetMean<double>(X); 
   double STDEV = gm2fieldUtil::Math::GetStandardDeviation<double>(X); 

   double yMin = MEAN - 3.*STDEV;   
   double yMax = MEAN + 3.*STDEV;   

   const int NL = sccOn.size();
   TLine **tON  = new TLine*[NL];
   TLine **tOFF = new TLine*[NL];
   for(int i=0;i<NL;i++){
      tON[i] = new TLine(sccOn[i],yMin,sccOn[i],yMax);
      tON[i]->SetLineColor(kGreen+1);
      tON[i]->SetLineWidth(2);
      tON[i]->SetLineStyle(2);
      tOFF[i] = new TLine(sccOff[i],yMin,sccOff[i],yMax);
      tOFF[i]->SetLineColor(kRed);
      tOFF[i]->SetLineWidth(2);
      tOFF[i]->SetLineStyle(2);

   }

   if(isBlind) ApplyBlindingTRLY(blindValue,trlyData);

   // get the mean field with SCC off and on 
   int nev = 30;  // keep 30 events in the analysis
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
      rc = GetDifference_ABA_sccFirst(useTimeWeight,sccTime,scc,sccErr,bareTime,bare,bareErr,diff_aba,diffErr_aba);  
   }else{
      rc = GetDifference_ABA(useTimeWeight,sccTime,scc,sccErr,bareTime,bare,bareErr,diff_aba,diffErr_aba);  
   }

   double mean_aba  = gm2fieldUtil::Math::GetMean<double>(diff_aba); 
   double stdev_aba = gm2fieldUtil::Math::GetStandardDeviation<double>(diff_aba); 
 
   const int ND = diff.size();
   std::vector<double> trial; 
   for(int i=0;i<ND;i++) trial.push_back(i+1);
   
   const int NDA = diff_aba.size();
   std::vector<double> trial_aba;
   for(int i=0;i<NDA;i++) trial_aba.push_back(i+1);

   if(ND==1){
      stdev     = diffErr[0]; 
      mean_aba  = 0;
      stdev_aba = 0;
   }

   dB[0]     = mean;
   dB_err[0] = stdev;

   dB[1]     = mean_aba;
   dB_err[1] = stdev_aba;  
   
   // Plots

   TGraph *gTR               = GetTRLYTGraph(probeNumber-1,"GpsTimeStamp","freq",trlyData);
   TGraphErrors *gTRLY_bare  = gm2fieldUtil::Graph::GetTGraphErrors(bareTime ,bare      ,bareErr      );
   TGraphErrors *gTRLY_scc   = gm2fieldUtil::Graph::GetTGraphErrors(sccTime  ,scc       ,sccErr       );
   TGraphErrors *gDiff       = gm2fieldUtil::Graph::GetTGraphErrors(trial    ,diff      ,diffErr      );  
   TGraphErrors *gDiff_aba   = gm2fieldUtil::Graph::GetTGraphErrors(trial_aba,diff_aba  ,diffErr_aba  );  

   gm2fieldUtil::Graph::SetGraphParameters(gTR       ,20,kBlack  );
   gm2fieldUtil::Graph::SetGraphParameters(gTRLY_bare,20,kBlack  );
   gm2fieldUtil::Graph::SetGraphParameters(gTRLY_scc ,20,kRed    );
   gm2fieldUtil::Graph::SetGraphParameters(gDiff     ,20,kBlue   );
   gm2fieldUtil::Graph::SetGraphParameters(gDiff_aba ,20,kGreen+2);

   TMultiGraph *mgTRLY = new TMultiGraph();
   mgTRLY->Add(gTRLY_bare,"lp");
   mgTRLY->Add(gTRLY_scc ,"lp");

   TMultiGraph *mgDiff = new TMultiGraph();
   mgDiff->Add(gDiff      ,"lp");
   mgDiff->Add(gDiff_aba  ,"lp");

   std::cout << Form("====================== RESULTS FOR PROBE %02d ======================",probeNumber) << std::endl;
   std::cout << "Raw results: " << std::endl;
   std::cout << Form("dB (%s): %.3lf +/- %.3lf Hz (%.3lf +/- %.3lf ppb)",gradName.c_str(),
                     dB[0],dB_err[0],dB[0]/0.06179,dB_err[0]/0.06179) << std::endl;
   std::cout << "-----------------------------------------------------" << std::endl; 
   std::cout << "Drift corrected [ABA]: " << std::endl;
   std::cout << Form("dB (%s): %.3lf +/- %.3lf Hz (%.3lf +/- %.3lf ppb)",gradName.c_str(),
                     dB[1],dB_err[1],dB[1]/0.06179,dB_err[1]/0.06179) << std::endl;

   rc = PrintToFile(outpath,gradName,dB,dB_err,drift,drift_err); 

   char plotDir[200];
   sprintf(plotDir,"./plots/%s",theDate.getDateString().c_str()); 
   rc = MakeDirectory(plotDir); 

   // draw some plots 
   TCanvas *c1 = new TCanvas("c1","Data",1200,600);
   c1->Divide(1,2); 

   c1->cd(1);
   mgTRLY->Draw("alp");
   gm2fieldUtil::Graph::SetGraphLabels(mgTRLY,"TRLY Data (Bare = Black, Red = SCC)","","Frequency (Hz)"); 
   gm2fieldUtil::Graph::UseTimeDisplay(mgTRLY);
   gm2fieldUtil::Graph::SetGraphLabelSizes(mgTRLY,0.05,0.06); 
   // mgTRLY_bare->GetXaxis()->SetLimits(xMin_bare,xMax_bare); 
   mgTRLY->Draw("alp");
   c1->Update();

   c1->cd(2);
   mgDiff->Draw("alp");
   gm2fieldUtil::Graph::SetGraphLabels(mgDiff,"SCC-Bare (Blue = Raw, Green = ABA)","Trial","Frequency Difference (Hz)"); 
   gm2fieldUtil::Graph::SetGraphLabelSizes(mgDiff,0.05,0.06); 
   // mgTRLY_bare->GetXaxis()->SetLimits(xMin_bare,xMax_bare); 
   mgDiff->Draw("alp");
   c1->Update();

   c1->cd();
   TString plotPath = Form("%s/trly_dB_%s-grad_run-%d_pr-%02d.png",plotDir,gradName.c_str(),run[0],probeNumber); 
   c1->Print(plotPath);
   delete c1;  

   TCanvas *c2 = new TCanvas("c1","Data",1200,600);

   c2->cd();
   gTR->Draw("alp");
   gm2fieldUtil::Graph::SetGraphLabels(gTR,Form("TRLY %02d Data",probeNumber),"","Frequency (Hz)"); 
   gm2fieldUtil::Graph::UseTimeDisplay(gTR);
   gTR->Draw("alp");
   for(int i=0;i<NL;i++){
      tON[i]->Draw("same"); 
      tOFF[i]->Draw("same"); 
   }
   c2->Update();

   c2->cd();
   plotPath = Form("%s/trly_dB_%s-grad_run-%d_pr-%02d_all-data.png",plotDir,gradName.c_str(),run[0],probeNumber); 
   c2->Print(plotPath);
   delete c2;  

   delete inputMgr; 

   return 0;
}
