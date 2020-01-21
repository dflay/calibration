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
#include "MovingAverage.h"
#include "Blinder.h"

#include "./include/date.h"
#include "./include/Constants.h"
#include "./include/drift.h"
#include "./include/fixedProbeEvent.h"
#include "./include/plungingProbeAnaEvent.h"
#include "./include/trolleyAnaEvent.h"

#include "./src/CustomUtilities.C"
#include "./src/CustomMath.C"
#include "./src/CustomImport.C"
#include "./src/CustomExport.C"
#include "./src/CustomGraph.C"
#include "./src/CustomAlgorithms.C"
#include "./src/OscFuncs.C"
#include "./src/BlindFuncs.C"
#include "./src/InputManager.C"
#include "./src/FitFuncs.C"
#include "./src/TRLYFuncs.C"
#include "./src/SystFuncs.C"

int DeltaB_trly_prod(std::string configFile){

   std::cout << "----------------------------------" << std::endl;
   std::cout << "TRLY DELTA B CALCULATION (ABA Method)" << std::endl;

   int rc=0;
   int method     = gm2fieldUtil::Constants::kPhaseDerivative;
   int fxprMethod = gm2fieldUtil::Constants::kPhaseDerivative; 

   InputManager *inputMgr = new InputManager();
   inputMgr->UseAxis(); 
   inputMgr->Load(configFile);
   inputMgr->Print();

   std::string prodVersion = inputMgr->GetProductionTag(); 
   std::string date        = inputMgr->GetAnalysisDate();
   std::string blindLabel  = inputMgr->GetBlindLabel();

   bool isBlind            = inputMgr->IsBlind();
   bool useTimeWeight      = inputMgr->GetTimeWeightStatus(); 
   bool useOscCor          = false; // inputMgr->GetOscCorStatus(); // never use oscillation corrections here   
   bool loadTimes          = inputMgr->GetSCCTimeStatus(); 
   int probeNumber         = inputMgr->GetTrolleyProbe();
   int axis                = inputMgr->GetAxis();
   int runPeriod           = inputMgr->GetRunPeriod();
   int nev                 = inputMgr->GetNumEventsToAvg();  
   // systematics 
   bool isSyst             = inputMgr->GetSystStatus(); 
   bool varyDB_time        = inputMgr->GetVaryTimeStatus("tr","db");
   double dB_delta         = inputMgr->GetDeltaTime("tr","db");  

   // update the analysis method according to Ran's guidance 
   if(prodVersion.compare("v9_21_01")==0) method = gm2fieldUtil::Constants::kHilbertPhaseLinear;

   date_t theDate; 
   GetDate(theDate);

   std::string plotDir = GetPath("plots" ,isBlind,blindLabel,theDate.getDateString());
   std::string outDir  = GetPath("output",isBlind,blindLabel,theDate.getDateString());
 
   std::string gradName;
   if(axis==0) gradName = "rad"; 
   if(axis==1) gradName = "vert"; 
   if(axis==2) gradName = "azi";

   std::string trLabel; 
   if(axis==0) trLabel = "trx"; 
   if(axis==1) trLabel = "try"; 
   if(axis==2) trLabel = "trz";

   char outpath[200];
   sprintf(outpath,"%s/dB-trly_final-location_%s-grad_pr-%02d.csv",outDir.c_str(),gradName.c_str(),probeNumber);

   int blindUnits  = inputMgr->GetBlindUnits(); 
   double blindMag = inputMgr->GetBlindScale(); 
   gm2fieldUtil::Blinder *myBlind = new gm2fieldUtil::Blinder(blindLabel,blindMag,blindUnits);
   double blindValue = myBlind->GetBlinding(2); // in Hz

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

   int coilSet  = 1;             // 0 = bottom, 1 = top, -1 = azi  
   double thr   = 100E-3;        // in A 
   double delta = 2.*60. + 40.;  // 2 min 40 sec   // for x, y 
   
   if(axis==2){
      coilSet = -1;
      delta = 1.*60.;
   }

   if(runPeriod==2){
      delta = 2.*50.;
      // if(axis==2) delta = 2.*60. + 12.; 
   }

   // determine the correct ordering of the SCC on/off cycles 
   std::vector<surfaceCoilEvent_t> sccData;
   if(!loadOnline){
      for(int i=0;i<NRUN;i++) rc = GetSurfaceCoilData(run[i],sccData,prodVersion);
      if(rc!=0) return 1;
   } 
   
   bool sccStartOn = false; 
   std::vector<double> sccOff,sccOn;
  
   std::vector<deltab_t> trly_dB; 

   char inpath_db[200];
   sprintf(inpath_db,"./input/delta-b/trly_xyz_run-%d.csv",runPeriod);

   double dB[3]        = {0,0,0}; 
   double dB_err[3]    = {0,0,0};
   double drift[3]     = {0,0,0};  
   double drift_err[3] = {0,0,0};  

   if(loadOnline){
      // use online results since we don't have a run to work with 
      rc    = LoadDeltaBData_trlyXYZ(inpath_db,probeNumber,trly_dB);
      dB[0] = trly_dB[axis].dB; 
      dB[1] = trly_dB[axis].dB_fxpr; 
      dB[2] = 0.; 
      dB_err[0] = trly_dB[axis].dB_err; 
      dB_err[1] = trly_dB[axis].dB_fxpr_err; 
      dB_err[2] = 0.; 
      rc = PrintToFile(outpath,gradName,dB,dB_err,drift,drift_err); 
      return 0;
   }else{
      if(loadTimes){
	 // better to use pre-defined transition times   
	 rc = LoadSCCTimes(probeNumber,runPeriod,prodVersion,trLabel,sccOff,sccOn);
	 if(sccOn[0]<sccOff[0]) sccStartOn = true;
      }else{
	 rc = FindTransitionTimes(coilSet,axis,thr,delta,sccData,sccOff,sccOn);
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

   if(isSyst && varyDB_time){
      std::cout << "[DeltaB_trly_prod]: SYSTEMATIC VARIATION! Delta-B times will be randomized by up to " << dB_delta << " sec" <<std::endl;
      rc = systFunc::RandomizeTimeValues(dB_delta,sccOff); 
      rc = systFunc::RandomizeTimeValues(dB_delta,sccOn ); 
   }

   // save the times -- DO NOT NEED AFTER FIRST ANA PASS
   // char scc_time_path[200];                                                              
   // sprintf(scc_time_path,"./input/scc-times/run-%d/%s-%02d.txt",runPeriod,trLabel.c_str(),probeNumber);
   // if(prodVersion.compare("nearline")==0) rc = PrintToFile_sccTimes(scc_time_path,sccOff,sccOn);

   // TRLY data
   std::vector<trolleyAnaEvent_t> trlyData;
   for(int i=0;i<NRUN;i++){
      rc = GetTrolleyData(run[i],method,trlyData,prodVersion);
      if(rc!=0){
	 std::cout << "No data!" << std::endl;
	 return 1;
      }
   }
   if(isBlind) ApplyBlindingTRLY(blindValue,trlyData);

   std::vector<int> fxprList;
   inputMgr->GetFXPRList(fxprList);

   std::vector<averageFixedProbeEvent_t> fxprData;
   bool subtractDrift = inputMgr->GetFXPRDriftStatus();  
   int period = inputMgr->GetNumEventsTimeWindow();
   for(int i=0;i<NRUN;i++){
      rc = GetFixedProbeData_avg(run[i],fxprMethod,fxprList,fxprData,prodVersion,subtractDrift,period,0);
      if(rc!=0){
	 std::cout << "No data!" << std::endl;
	 return 1;
      }
   }

   TGraphErrors *gFXPR = GetFXPRTGraph_avg("GpsTimeStamp","freq","NONE",fxprData);
   gm2fieldUtil::Graph::SetGraphParameters(gFXPR,20,kBlack); 

   std::vector<double> X; 
   int NTR = trlyData.size();
   for(int i=0;i<NTR;i++) X.push_back( trlyData[i].freq[probeNumber-1] );
   double MEAN  = gm2fieldUtil::Math::GetMean<double>(X); 
   double STDEV = gm2fieldUtil::Math::GetStandardDeviation<double>(X); 

   double yMin = MEAN - 3.*STDEV;   
   double yMax = MEAN + 3.*STDEV;   

   double yMin_scc = -3; 
   double yMax_scc =  3; 

   const int NL     = sccOn.size();
   TLine **tON      = new TLine*[NL];
   TLine **tOFF     = new TLine*[NL];
   TLine **tON_scc  = new TLine*[NL];
   TLine **tOFF_scc = new TLine*[NL];
   for(int i=0;i<NL;i++){
      tON[i] = new TLine(sccOn[i],yMin,sccOn[i],yMax);
      tON[i]->SetLineColor(kGreen+1);
      tON[i]->SetLineWidth(2);
      tON[i]->SetLineStyle(2);
      tOFF[i] = new TLine(sccOff[i],yMin,sccOff[i],yMax);
      tOFF[i]->SetLineColor(kRed);
      tOFF[i]->SetLineWidth(2);
      tOFF[i]->SetLineStyle(2);
      // for the SCC plot 
      tON_scc[i] = new TLine(sccOn[i],yMin_scc,sccOn[i],yMax_scc);
      tON_scc[i]->SetLineColor(kGreen+1);
      tON_scc[i]->SetLineWidth(2);
      tON_scc[i]->SetLineStyle(2);
      tOFF_scc[i] = new TLine(sccOff[i],yMin_scc,sccOff[i],yMax_scc);
      tOFF_scc[i]->SetLineColor(kRed);
      tOFF_scc[i]->SetLineWidth(2);
      tOFF_scc[i]->SetLineStyle(2);
   }

   // find the reference time t0 for oscillation corrections
   // WARNING: Can't do this here because we're CHANGING THE FIELD using SCC!  
   double t0=-1;

   // get the mean field with SCC off and on 
   std::vector<double> bareTime,bare,bareErr,sccTime,scc,sccErr;  
   if(useOscCor) std::cout << "[DeltaB_trly_prod]: Bare field data" << std::endl;
   rc = GetTRLYStats_sccToggle(useOscCor,probeNumber-1,nev,sccOff,fxprData,trlyData,bareTime,bare,bareErr,t0); 
   if(useOscCor) std::cout << "[DeltaB_trly_prod]: Grad field data" << std::endl;
   rc = GetTRLYStats_sccToggle(useOscCor,probeNumber-1,nev,sccOn ,fxprData,trlyData,sccTime ,scc ,sccErr ,t0); 

   // do delta B calcs
   // raw difference 
   std::vector<double> diff,diffErr;
   rc = GetDifference(scc,sccErr,bare,bareErr,diff,diffErr); 

   double mean=0,stdev=0,err=0; 
   rc = GetWeightedAverageStats(diff,diffErr,mean,err,stdev); 

   int NN = bare.size(); 

   // ABA difference 
   double mean_aba=0,stdev_aba=0; 
   std::vector<double> diff_aba,diffErr_aba;
   if(NN>1){
      if(sccStartOn){
	 // SCC was first! A = SCC, B = baseline  
	 rc = GetDifference_ABA_final(useTimeWeight,sccTime,scc,sccErr,bareTime,bare,bareErr,diff_aba,diffErr_aba);  
      }else{
	 // Baseline was first! A = baseline, B = SCC  
	 rc = GetDifference_ABA_final(useTimeWeight,bareTime,bare,bareErr,sccTime,scc,sccErr,diff_aba,diffErr_aba); 
	 // need to invert results since we want SCC - baseline 
	 for(int i=0;i<NN;i++) diff_aba[i] *= -1.;  
      }
      rc = GetWeightedAverageStats(diff_aba,diffErr_aba,mean_aba,err,stdev_aba); 
   }else{
      // not enough events!
      mean_aba  = 0;
      stdev_aba = 0; 
   }
 
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
  
   // if we have a single ABA trial, use statistical uncertainty as the error 
   if(NDA==1) dB_err[1] = diffErr_aba[0];
 
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
   if(NN>1) mgDiff->Add(gDiff_aba  ,"lp");

   std::cout << Form("====================== RESULTS FOR PROBE %02d ======================",probeNumber) << std::endl;
   std::cout << "Raw results: " << std::endl;
   std::cout << Form("dB (%s): %.3lf +/- %.3lf Hz",gradName.c_str(),dB[0],dB_err[0]) << std::endl;
   std::cout << "-----------------------------------------------------" << std::endl; 
   std::cout << "Drift corrected [ABA]: " << std::endl;
   std::cout << Form("dB (%s): %.3lf +/- %.3lf Hz",gradName.c_str(),dB[1],dB_err[1]) << std::endl;

   rc = PrintToFile(outpath,gradName,dB,dB_err,drift,drift_err); 

   char msg[200];
   if(dB[1]==0){
      sprintf(msg,"[DeltaB_trly_prod]: No ABA data for probe %02d, axis %d!",probeNumber,axis);
      Logger::PrintMessage(Logger::kERR,"default",msg,'a');
   }

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
   TString plotPath = Form("%s/trly_dB_%s-grad_pr-%02d.png",plotDir.c_str(),gradName.c_str(),probeNumber); 
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
   plotPath = Form("%s/trly_dB_%s-grad_pr-%02d_all-data.png",plotDir.c_str(),gradName.c_str(),probeNumber); 
   c2->Print(plotPath);
   delete c2;  

   TGraph *gSCCb = GetSCCPlot( 0,sccData);
   TGraph *gSCCt = GetSCCPlot( 1,sccData);
   TGraph *gSCCa = GetSCCPlot(-1,sccData);

   gm2fieldUtil::Graph::SetGraphParameters(gSCCb,20,kBlack);
   gm2fieldUtil::Graph::SetGraphParameters(gSCCt,20,kRed  );
   gm2fieldUtil::Graph::SetGraphParameters(gSCCa,20,kBlue );

   TMultiGraph *mgSCC = new TMultiGraph();
   mgSCC->Add(gSCCb,"lp");
   mgSCC->Add(gSCCt,"lp");
   mgSCC->Add(gSCCa,"lp");

   TString Title_scc = Form("SCC Data (black = bot, red = top, blue = azi)");

   TCanvas *c3 = new TCanvas("c3","SCC Data",1200,600);
   c3->cd();

   mgSCC->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mgSCC,Title_scc,"","Current (A)");
   gm2fieldUtil::Graph::UseTimeDisplay(mgSCC);
   mgSCC->GetYaxis()->SetRangeUser(yMin_scc,yMax_scc); 
   mgSCC->Draw("a");
   for(int i=0;i<NL;i++){
      tON_scc[i]->Draw("same");
      tOFF_scc[i]->Draw("same");
   } 
   c3->Update(); 

   c3->cd();
   plotPath = Form("%s/trly_dB_scc-currents_%s-grad_pr-%02d.png",plotDir.c_str(),gradName.c_str(),probeNumber); 
   c3->Print(plotPath);
   delete c3;

   // for the sake of making plots
   std::vector<double> trTime_bare,trFreq_bare,trFreq_cor_bare;
   std::vector<double> trTime_scc ,trFreq_scc ,trFreq_cor_scc;
   rc = CorrectOscillation_trly(probeNumber-1,nev,bareTime,fxprData,trlyData,trTime_bare,trFreq_bare,trFreq_cor_bare);
   rc = CorrectOscillation_trly(probeNumber-1,nev,sccTime ,fxprData,trlyData,trTime_scc,trFreq_scc,trFreq_cor_scc);
 
   // to make it easier to see, subtract off the linear trend
   double f0_bare = trFreq_bare[0];
   double f0_scc  = trFreq_scc[0];  
   int SIZE = trTime_bare.size(); 
   double intercept=0,slope=0,r=0;
   rc = gm2fieldUtil::Math::LeastSquaresFitting(trTime_bare,trFreq_bare,intercept,slope,r);
   for(int i=0;i<SIZE;i++) trFreq_bare[i] -= -f0_bare + (intercept + slope*trTime_bare[i]); 
   rc = gm2fieldUtil::Math::LeastSquaresFitting(trTime_bare,trFreq_cor_bare,intercept,slope,r);
   for(int i=0;i<SIZE;i++) trFreq_cor_bare[i] -= -f0_bare + (intercept + slope*trTime_bare[i]); 

   SIZE = trTime_scc.size(); 
   rc = gm2fieldUtil::Math::LeastSquaresFitting(trTime_scc,trFreq_scc,intercept,slope,r);
   for(int i=0;i<SIZE;i++) trFreq_scc[i] -= -f0_scc + (intercept + slope*trTime_scc[i]); 
   rc = gm2fieldUtil::Math::LeastSquaresFitting(trTime_scc,trFreq_cor_scc,intercept,slope,r);
   for(int i=0;i<SIZE;i++) trFreq_cor_scc[i] -= -f0_scc + (intercept + slope*trTime_scc[i]); 

   TGraph *gTR_bare_raw = gm2fieldUtil::Graph::GetTGraph(trTime_bare,trFreq_bare); 
   TGraph *gTR_bare_cor = gm2fieldUtil::Graph::GetTGraph(trTime_bare,trFreq_cor_bare); 
   TGraph *gTR_scc_raw  = gm2fieldUtil::Graph::GetTGraph(trTime_scc ,trFreq_scc); 
   TGraph *gTR_scc_cor  = gm2fieldUtil::Graph::GetTGraph(trTime_scc ,trFreq_cor_scc);

   gm2fieldUtil::Graph::SetGraphParameters(gTR_bare_raw,21,kBlack); 
   gm2fieldUtil::Graph::SetGraphParameters(gTR_bare_cor,20,kRed); 
   gm2fieldUtil::Graph::SetGraphParameters(gTR_scc_raw ,21,kBlack); 
   gm2fieldUtil::Graph::SetGraphParameters(gTR_scc_cor ,20,kRed);

   TMultiGraph *mgTR_bare = new TMultiGraph(); 
   mgTR_bare->Add(gTR_bare_raw,"lp");  
   if(useOscCor) mgTR_bare->Add(gTR_bare_cor,"lp");  
   
   TMultiGraph *mgTR_scc = new TMultiGraph(); 
   mgTR_scc->Add(gTR_scc_raw,"lp");  
   if(useOscCor) mgTR_scc->Add(gTR_scc_cor,"lp");  

   TString Title_tr = Form("TRLY %02d SCC OFF",probeNumber);
   if(useOscCor) Title_tr += Form(" (black = raw, red = osc cor)");
   TString yAxisTitle_tr = Form("Frequency (Hz)");

   TCanvas *c4 = new TCanvas("c4","TRLY Data",1200,600);
   c4->Divide(1,2);
 
   c4->cd(1);
   mgTR_bare->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mgTR_bare,Title_tr,"",yAxisTitle_tr);
   gm2fieldUtil::Graph::UseTimeDisplay(mgTR_bare);
   gm2fieldUtil::Graph::SetGraphLabelSizes(mgTR_bare,0.05,0.06); 
   mgTR_bare->Draw("a");
   c4->Update();

   c4->cd(2);
   mgTR_scc->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mgTR_scc,"SCC ON","",yAxisTitle_tr);
   gm2fieldUtil::Graph::UseTimeDisplay(mgTR_scc);
   gm2fieldUtil::Graph::SetGraphLabelSizes(mgTR_scc,0.05,0.06); 
   mgTR_scc->Draw("a");
   c4->Update();

   c4->cd(); 
   plotPath = Form("%s/trly_dB_%s-grad_all-events_pr-%02d.png",plotDir.c_str(),gradName.c_str(),probeNumber); 
   c4->Print(plotPath);
   delete c4;

   TCanvas *c5 = new TCanvas("c5","FXPR Data",1200,600);
   c5->cd();

   gFXPR->Draw("alp");
   gm2fieldUtil::Graph::SetGraphLabels(gFXPR,"FXPR Data","","Frequency (Hz)");
   gm2fieldUtil::Graph::UseTimeDisplay(gFXPR);
   gFXPR->Draw("alp");
   c5->Update(); 

   plotPath = Form("%s/trly_dB_%s-grad_fxpr-data_pr-%02d.png",plotDir.c_str(),gradName.c_str(),probeNumber); 
   c5->Print(plotPath);

   delete inputMgr; 

   return 0;
}
