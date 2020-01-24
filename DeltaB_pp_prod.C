// Compute DeltaB for PP data sets 
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
#include "./src/CustomImport.C"
#include "./src/CustomExport.C"
#include "./src/CustomMath.C"
#include "./src/CustomGraph.C"
#include "./src/CustomAlgorithms.C"
#include "./src/BlindFuncs.C"
#include "./src/InputManager.C"
#include "./src/FitFuncs.C"
#include "./src/OscFuncs.C"

int SortPPDAQRuns_alt(std::vector<plungingProbeAnaEvent_t> data,
                      std::vector<double> &tOdd ,std::vector<double> &fOdd,
                      std::vector<double> &tEven,std::vector<double> &fEven); 

int SortPPDAQRuns(std::vector<plungingProbeAnaEvent_t> data,bool sccStartOn,
                  std::vector<double> &sccTime ,std::vector<double> &scc ,std::vector<double> &sccErr,
                  std::vector<double> &bareTime,std::vector<double> &bare,std::vector<double> &bareErr);

int DeltaB_pp_prod(std::string configFile){

   std::cout << "----------------------------------" << std::endl;
   std::cout << "PP DELTA B CALCULATION (ABA Method)" << std::endl;

   int rc=0;
   int prMethod = gm2fieldUtil::Constants::kPhaseDerivative;
   int ppMethod = plungingProbeAnalysis::kLeastSquaresPhase;

   InputManager *inputMgr = new InputManager();
   inputMgr->UseAxis(); 
   inputMgr->Load(configFile);
   inputMgr->Print();

   std::string prodVersion   = inputMgr->GetProductionTag();
   std::string nmrAnaVersion = inputMgr->GetNMRANATag();  
   std::string blindLabel    = inputMgr->GetBlindLabel();
   std::string cutFile       = inputMgr->GetCutFile();

   bool isBlind              = inputMgr->IsBlind();
   bool useTimeWeight        = inputMgr->GetTimeWeightStatus();
   bool useOscCor            = false; // inputMgr->GetOscCorStatus(); // never use oscillation corrections here  
   int probeNumber           = inputMgr->GetTrolleyProbe();
   int axis                  = inputMgr->GetAxis();
   int runPeriod             = inputMgr->GetRunPeriod();
   // systematics
   bool isSyst               = inputMgr->GetSystStatus();
   int systDirNum            = inputMgr->GetSystDirNum(); 

   double tempCorValue = 0;
   bool useTempCor_pp  = inputMgr->GetTempCorStatus_pp();
   if(useTempCor_pp) tempCorValue = inputMgr->GetTempCor_pp(); 

   char cutPath[200];
   sprintf(cutPath,"./input/json/run-%d/%s",runPeriod,cutFile.c_str());
   std::string cutpath = cutPath;

   // date_t theDate; 
   // GetDate(theDate);
   std::string theDate = inputMgr->GetAnaDate();

   std::string plotDir = GetPath("plots" ,isBlind,blindLabel,theDate,isSyst,systDirNum);
   std::string outDir  = GetPath("output",isBlind,blindLabel,theDate,isSyst,systDirNum);

   int blindUnits  = inputMgr->GetBlindUnits(); 
   double blindMag = inputMgr->GetBlindScale(); 
   gm2fieldUtil::Blinder *myBlind = new gm2fieldUtil::Blinder(blindLabel,blindMag,blindUnits);
   double blindValue = myBlind->GetBlinding(1); // in Hz

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
   std::vector<surfaceCoilEvent_t> sccData;
   for(int i=0;i<NRUN;i++) rc = GetSurfaceCoilData(run[i],sccData,prodVersion);
   if(rc!=0) return 1; 

   // load SCC times or not? 
   std::string dB_key = "load-pp-scc-times"; 
   bool loadTimes = inputMgr->GetValueFromKey<bool>(dB_key); 
 
   int coilSet  = 1;
   double thr   = 100E-3;        // in A 
   double delta = 2.*60. + 40.;  // 2 min 40 sec   
   std::vector<double> sccOff,sccOn;
   bool sccStartOn = false;

   std::string ppLabel; 
   if(axis==0) ppLabel = "ppx"; 
   if(axis==1) ppLabel = "ppy"; 
   if(axis==2) ppLabel = "ppz"; 
   
   if(axis==2){
      coilSet = -1;
      delta   = 100.;  // 1 min 40 sec  
   }
 
   if(loadTimes){
      // better to use pre-defined transition times   
      rc = LoadSCCTimes(probeNumber,runPeriod,prodVersion,ppLabel,sccOff,sccOn);
      if(sccOn[0]<sccOff[0]) sccStartOn = true;
   }else{
      rc = FindTransitionTimes(coilSet,axis,thr,delta,sccData,sccOff,sccOn);
      if(rc<0){
	 std::cout << "No SCC transitions!" << std::endl;
	 return 1;
      }
      if(rc==1) sccStartOn = true;
   }

   if(sccStartOn){
      std::cout << "SCC was ON to start the sequence" << std::endl;
   }else{
      std::cout << "SCC was OFF to start the sequence" << std::endl;
   } 

   // save the times -- DO NOT NEED AFTER FIRST ANA PASS 
   // char scc_time_path[200]; 
   // sprintf(scc_time_path,"./input/scc-times/run-%d/%s-%02d.txt",runPeriod,ppLabel.c_str(),probeNumber);
   // if(prodVersion.compare("nearline")==0) rc = PrintToFile_sccTimes(scc_time_path,sccOff,sccOn);

   const int NSCC     = sccOn.size();

   // FXPR data 
   std::vector<int> fxprList;
   inputMgr->GetFXPRList(fxprList);

   std::vector<averageFixedProbeEvent_t> fxprData;
   bool subtractDrift = inputMgr->GetFXPRDriftStatus();  
   int period = inputMgr->GetNumEventsTimeWindow(); // 10;
   for(int i=0;i<NRUN;i++){
      rc = GetFixedProbeData_avg(run[i],prMethod,fxprList,fxprData,prodVersion,subtractDrift,period,0);
      if(rc!=0){
         std::cout << "No data!" << std::endl;
         return 1;
      }
   }

   TGraphErrors *gFXPR = GetFXPRTGraph_avg("GpsTimeStamp","freq","NONE",fxprData);
   gm2fieldUtil::Graph::SetGraphParameters(gFXPR,20,kBlack);

   // PP data
   bool useNMRANA = true; 
   const int N3 = run.size(); 
   std::vector<plungingProbeAnaEvent_t> ppInput,ppEvent; 
   for(int i=0;i<N3;i++){
      std::cout << "Getting PP data for run " << run[i] << "..." << std::endl; 
      rc = GetPlungingProbeData(run[i],prMethod,ppMethod,ppInput,prodVersion,nmrAnaVersion,cutpath,useNMRANA,tempCorValue);
      if(rc!=0){
	 std::cout << "No data!" << std::endl;
	 return 1;
      }
   }

   if(isBlind) ApplyBlindingPP(blindValue,ppInput);

   // for diagnostics 

   // reference time for oscillation correction
   // WARNING: Can't do this here because we're CHANGING THE FIELD using SCC!  
   double t0 = -1; 

   // oscillation correction 
   if(useOscCor){
      std::cout << "DOING OSCILLATION CORRECTION" << std::endl;
      rc = CorrectOscillation_pp(fxprData,ppInput,ppEvent,t0);
   }else{
      CopyPlungingProbe(ppInput,ppEvent);
   }
   
   const int NPP = ppEvent.size(); 

   // gather the PP DAQ runs to analyze
   std::vector<double> sccTime,scc,sccErr,bareTime,bare,bareErr; 
   rc = SortPPDAQRuns(ppEvent,sccStartOn,sccTime,scc,sccErr,bareTime,bare,bareErr);

   // for diagnostic plots 
   std::vector<double> tOdd_raw,fOdd_raw,tEven_raw,fEven_raw; 
   std::vector<double> tOdd_cor,fOdd_cor,tEven_cor,fEven_cor; 
   rc = SortPPDAQRuns_alt(ppInput,tOdd_raw,fOdd_raw,tEven_raw,fEven_raw);
   rc = SortPPDAQRuns_alt(ppEvent,tOdd_cor,fOdd_cor,tEven_cor,fEven_cor);

   TGraph *gOdd_raw  = gm2fieldUtil::Graph::GetTGraph(tOdd_raw ,fOdd_raw);
   TGraph *gOdd_cor  = gm2fieldUtil::Graph::GetTGraph(tOdd_cor ,fOdd_cor);
   TGraph *gEven_raw = gm2fieldUtil::Graph::GetTGraph(tEven_raw,fEven_raw);
   TGraph *gEven_cor = gm2fieldUtil::Graph::GetTGraph(tEven_cor,fEven_cor);

   gm2fieldUtil::Graph::SetGraphParameters(gOdd_raw,21,kBlack);
   gm2fieldUtil::Graph::SetGraphParameters(gOdd_cor,20,kRed);
   gm2fieldUtil::Graph::SetGraphParameters(gEven_raw,21,kBlack);
   gm2fieldUtil::Graph::SetGraphParameters(gEven_cor,20,kRed);

   TMultiGraph *mgo = new TMultiGraph();
   mgo->Add(gOdd_raw,"lp"); 
   mgo->Add(gOdd_cor,"lp"); 

   TMultiGraph *mge = new TMultiGraph();
   mge->Add(gEven_raw,"lp"); 
   mge->Add(gEven_cor,"lp"); 

   // do delta B calcs
   // raw difference 
   std::vector<double> diff,diffErr;
   rc = GetDifference(scc,sccErr,bare,bareErr,diff,diffErr);
  
   rc = CheckDifference(diff,diffErr);
   if(rc!=-1){
      SwapEntries(rc,scc,bare); 
      diff.clear();
      diffErr.clear();
      rc = GetDifference(scc,sccErr,bare,bareErr,diff,diffErr);
   }

   double mean=0,stdev=0,err=0;
   rc = GetWeightedAverageStats(diff,diffErr,mean,err,stdev); 

   double mean_aba=0,stdev_aba=0; 
   std::vector<double> diff_aba,diffErr_aba;
   if(NPP>1){
      // ABA difference  
      if(sccStartOn){
	 // SCC is first; A = SCC, B = baseline 
	 rc = GetDifference_ABA_final(useTimeWeight,sccTime,scc,sccErr,bareTime,bare,bareErr,diff_aba,diffErr_aba);  
      }else{
	 // baseline is first; A = baseline, B = SCC 
	 rc = GetDifference_ABA_final(useTimeWeight,bareTime,bare,bareErr,sccTime,scc,sccErr,diff_aba,diffErr_aba); 
	 // need to flip the sign; we want SCC - baseline 
	 for(int i=0;i<NPP;i++) diff_aba[i] *= -1.;  
      }
      rc = GetWeightedAverageStats(diff_aba,diffErr_aba,mean_aba,err,stdev_aba); 
   }else{
      // not enough events! 
      mean_aba  = 0;
      stdev_aba = 0;
   }
 
   double dB[3]        = {0,0,0}; 
   double dB_err[3]    = {0,0,0};
   double drift[3]     = {0,0,0};  
   double drift_err[3] = {0,0,0};  

   dB[0]     = mean;
   dB_err[0] = stdev;

   const int ND = diff.size();
   std::vector<double> trial; 
   for(int i=0;i<ND;i++){
      trial.push_back(i+1);
      std::cout << Form("RAW trial %02d: %.3lf +/- %.3lf Hz",i+1,diff[i],diffErr[i]) << std::endl;
   }

   const int NDA = diff_aba.size();
   std::vector<double> trial_aba;
   for(int i=0;i<NDA;i++){
      trial_aba.push_back(i+1);
      std::cout << Form("ABA trial %02d: %.3lf +/- %.3lf Hz",i+1,diff_aba[i],diffErr_aba[i]) << std::endl;
   }
  
   // single trial, use shot uncertainty 
   if(ND==1){
      dB_err[0] = diffErr[0]; 
   }

   dB[1]     = mean_aba;
   dB_err[1] = stdev_aba; 

   // if we have a single ABA trial, use statistical uncertainty as the error 
   if(NDA==1) dB_err[1] = diffErr_aba[0];  

   // clean up junk results
   for(int i=0;i<3;i++) if( gm2fieldUtil::Math::IsInfOrNaN<double>(dB[i]) )     dB[i] = 0.;
   for(int i=0;i<3;i++) if( gm2fieldUtil::Math::IsInfOrNaN<double>(dB_err[i]) ) dB_err[i] = 0.;

   // Plots

   TGraph *gPP             = GetPPTGraph1("TimeStamp","freq",ppInput); 
   TGraph *gPP_oscCor      = GetPPTGraph1("TimeStamp","freq",ppEvent); 

   TGraphErrors *gPP_bare  = gm2fieldUtil::Graph::GetTGraphErrors(bareTime ,bare      ,bareErr      );
   TGraphErrors *gPP_scc   = gm2fieldUtil::Graph::GetTGraphErrors(sccTime  ,scc       ,sccErr       );
   TGraphErrors *gDiff     = gm2fieldUtil::Graph::GetTGraphErrors(trial    ,diff      ,diffErr      );  
   TGraphErrors *gDiff_aba = gm2fieldUtil::Graph::GetTGraphErrors(trial_aba,diff_aba  ,diffErr_aba  );  

   gm2fieldUtil::Graph::SetGraphParameters(gPP,20    ,kBlack); 
   gm2fieldUtil::Graph::SetGraphParameters(gPP_oscCor,20,kRed); 
  
   gm2fieldUtil::Graph::SetGraphParameters(gPP_bare   ,20,kBlack  );
   gm2fieldUtil::Graph::SetGraphParameters(gPP_scc    ,20,kRed    );
   gm2fieldUtil::Graph::SetGraphParameters(gDiff      ,20,kBlue   );
   gm2fieldUtil::Graph::SetGraphParameters(gDiff_aba  ,20,kGreen+2);

   TMultiGraph *mg = new TMultiGraph();
   mg->Add(gPP       ,"lp"); 
   mg->Add(gPP_oscCor,"lp"); 

   TMultiGraph *mgPP = new TMultiGraph();
   mgPP->Add(gPP_bare,"lp");
   mgPP->Add(gPP_scc ,"lp");

   TMultiGraph *mgDiff = new TMultiGraph();
   mgDiff->Add(gDiff      ,"lp");
   if(NPP>1) mgDiff->Add(gDiff_aba  ,"lp");

   TLegend *L = new TLegend(0.6,0.6,0.8,0.8); 
   L->AddEntry(gPP_bare   ,"Bare","p"); 
   L->AddEntry(gPP_scc    ,"SCC" ,"p");
   
   TLegend *LD = new TLegend(0.6,0.6,0.8,0.8); 
   LD->AddEntry(gDiff      ,"Raw"  ,"p");  
   if(NPP>1) LD->AddEntry(gDiff_aba  ,"ABA"  ,"p");  

   std::cout << Form("====================== RESULTS FOR PROBE %02d ======================",probeNumber) << std::endl;
   std::cout << "Raw results: " << std::endl;
   std::cout << Form("dB (%s): %.3lf +/- %.3lf Hz",gradName.c_str(),dB[0],dB_err[0]) << std::endl;
   std::cout << "-----------------------------------------------------" << std::endl; 
   std::cout << "Drift corrected [ABA]: " << std::endl;
   std::cout << Form("dB (%s): %.3lf +/- %.3lf Hz",gradName.c_str(),dB[1],dB_err[1]) << std::endl;

   char msg[200]; 
   if(dB[1]==0){
      sprintf(msg,"[DeltaB_pp_prod]: No ABA data for probe %02d, axis %d!",probeNumber,axis);
      Logger::PrintMessage(Logger::kERR,"default",msg,'a');  
   }

   char outpath[200];
   sprintf(outpath,"%s/dB-pp_final-location_%s-grad_pr-%02d.csv",outDir.c_str(),gradName.c_str(),probeNumber);
   rc = PrintToFile(outpath,gradName,dB,dB_err,drift,drift_err); 

   // draw some plots 
   TCanvas *c1 = new TCanvas("c1","Data",1200,600);
   c1->Divide(1,2); 

   c1->cd(1);
   mgPP->Draw("alp");
   gm2fieldUtil::Graph::SetGraphLabels(mgPP,"PP Data (Bare = Black, Red = SCC)","","Frequency (Hz)"); 
   gm2fieldUtil::Graph::UseTimeDisplay(mgPP);
   gm2fieldUtil::Graph::SetGraphLabelSizes(mgPP,0.05,0.06);  
   // mgPP_bare->GetXaxis()->SetLimits(xMin_bare,xMax_bare); 
   mgPP->Draw("alp");
   // L->Draw("same"); 
   c1->Update();

   c1->cd(2);
   mgDiff->Draw("alp");
   gm2fieldUtil::Graph::SetGraphLabels(mgDiff,"SCC-Bare (Raw = Blue, ABA = Green)","Trial","Frequency Difference (Hz)");
   gm2fieldUtil::Graph::SetGraphLabelSizes(mgDiff,0.05,0.06);  
   // mgPP_bare->GetXaxis()->SetLimits(xMin_bare,xMax_bare); 
   mgDiff->Draw("alp");
   // LD->Draw("same"); 
   c1->Update();

   c1->cd();
   TString plotPath = Form("%s/pp_dB_%s-grad_pr-%02d.png",plotDir.c_str(),gradName.c_str(),probeNumber); 
   c1->Print(plotPath);
   delete c1;  

   double yMin_scc = -50; 
   double yMax_scc =  50; 

   TLine **tON  = new TLine*[NSCC];
   TLine **tOFF = new TLine*[NSCC];
   for(int i=0;i<NSCC;i++){
      tON[i] = new TLine(sccOn[i],yMin_scc,sccOn[i],yMax_scc);
      tON[i]->SetLineColor(kGreen+1);
      tON[i]->SetLineWidth(2);
      tON[i]->SetLineStyle(2);
      tOFF[i] = new TLine(sccOff[i],yMax_scc,sccOff[i],yMin_scc);
      tOFF[i]->SetLineColor(kRed);
      tOFF[i]->SetLineWidth(2);
      tOFF[i]->SetLineStyle(2);
   }

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

   TCanvas *c2 = new TCanvas("c2","SCC Data",1200,600);
   c2->cd();

   mgSCC->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mgSCC,Title_scc,"","Current (A)");
   gm2fieldUtil::Graph::UseTimeDisplay(mgSCC);
   mgSCC->GetYaxis()->SetRangeUser(yMin_scc,yMax_scc);
   mgSCC->Draw("a");
   for(int i=0;i<NSCC;i++) tON[i]->Draw("same"); 
   for(int i=0;i<NSCC;i++) tOFF[i]->Draw("same"); 
   c2->Update(); 

   c2->cd();
   plotPath = Form("%s/pp_dB_scc-currents_%s-grad_pr-%02d.png",plotDir.c_str(),gradName.c_str(),probeNumber); 
   c2->Print(plotPath);
   delete c2; 

   TString Title_o;
   TString Title_e = Form("All PP Data");
   if(sccStartOn){
      Title_e += Form("SCC ON");
      Title_o = Form("SCC OFF");
   }else{
      Title_e += Form("SCC OFF");
      Title_o = Form("SCC ON");
   }
   if(useOscCor) Title_e += Form(" (black = raw, red = osc cor)");
   TString yAxisTitle_pp = Form("Frequency (Hz)");

   TCanvas *c3 = new TCanvas("c3","PP Data",1200,600);
   c3->Divide(1,2); 
 
   c3->cd(1);
   mge->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mge,Title_e,"",yAxisTitle_pp);
   gm2fieldUtil::Graph::UseTimeDisplay(mge);
   mge->Draw("a");
   c3->Update(); 

   c3->cd(2);
   mgo->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mgo,Title_o,"",yAxisTitle_pp);
   gm2fieldUtil::Graph::UseTimeDisplay(mgo);
   mgo->Draw("a");
   c3->Update(); 
   
   
   plotPath = Form("%s/pp_dB_%s-grad_all-events_pr-%02d.png",plotDir.c_str(),gradName.c_str(),probeNumber); 
   c3->Print(plotPath);
   delete c3;  

   TCanvas *c4 = new TCanvas("c4","FXPR Data",1200,600);
   c4->cd();

   gFXPR->Draw("alp");
   gm2fieldUtil::Graph::SetGraphLabels(gFXPR,"FXPR Data","","Frequency (Hz)");
   gm2fieldUtil::Graph::UseTimeDisplay(gFXPR); 
   gFXPR->Draw("alp");
   c4->Update();

   plotPath = Form("%s/pp_dB_%s-grad_fxpr-data_pr-%02d.png",plotDir.c_str(),gradName.c_str(),probeNumber);
   c4->Print(plotPath);

   delete inputMgr;
 
   return 0;
}
//______________________________________________________________________________
int SortPPDAQRuns_alt(std::vector<plungingProbeAnaEvent_t> data,
                      std::vector<double> &tOdd ,std::vector<double> &fOdd,
                      std::vector<double> &tEven,std::vector<double> &fEven){

   int M=0;
   const int N = data.size();  
   for(int i=0;i<N;i++){
      // compute the average for the PP DAQ run 
      std::cout << "Processing PP DAQ run " << data[i].run << std::endl; 
      M = data[i].numTraces; 
      for(int j=0;j<M;j++){
	 if(i%2!=0){
	    // odd run 
	    tOdd.push_back(data[i].time[j]/1E+9);  
	    fOdd.push_back(data[i].freq[j]); 
	 }else{
	    // even run   
	    tEven.push_back(data[i].time[j]/1E+9);  
	    fEven.push_back(data[i].freq[j]);
	 }
      } 
   }
   return 0;
}
//______________________________________________________________________________
int SortPPDAQRuns(std::vector<plungingProbeAnaEvent_t> data,bool sccStartOn,
                  std::vector<double> &sccTime ,std::vector<double> &scc ,std::vector<double> &sccErr,
                  std::vector<double> &bareTime,std::vector<double> &bare,std::vector<double> &bareErr){

   int M=0;
   double mean=0,stdev=0;
   std::vector<double> x;
   std::vector<double> tOdd,Odd,OddErr; 
   std::vector<double> tEven,Even,EvenErr; 
   const int N = data.size();  
   for(int i=0;i<N;i++){
      // compute the average for the PP DAQ run 
      std::cout << "Processing PP DAQ run " << data[i].run << std::endl; 
      M = data[i].numTraces; 
      for(int j=0;j<M;j++) x.push_back(data[i].freq[j]); 
      mean  = gm2fieldUtil::Math::GetMean<double>(x); 
      stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(x); 
      if(i%2!=0){
	 // odd run 
	 tOdd.push_back(data[i].time[0]/1E+9);  
         Odd.push_back(mean); 
	 OddErr.push_back(stdev);
      }else{
	 // even run   
	 tEven.push_back(data[i].time[0]/1E+9);  
	 Even.push_back(mean);
	 EvenErr.push_back(stdev);
      }
      // clean up for next PP DAQ run 
      x.clear(); 
   }

   const int NN = tOdd.size(); 

   for(int i=0;i<NN;i++){
      if(sccStartOn){
	 // if SCC is ON, evens are SCC 
	 sccTime.push_back(tEven[i]);  
	 scc.push_back(Even[i]);  
	 sccErr.push_back(EvenErr[i]);  
	 bareTime.push_back(tOdd[i]);  
	 bare.push_back(Odd[i]);  
	 bareErr.push_back(OddErr[i]);  
      }else{
	 // if SCC is ON, odds are SCC 
	 bareTime.push_back(tEven[i]);  
	 bare.push_back(Even[i]);  
	 bareErr.push_back(EvenErr[i]);  
	 sccTime.push_back(tOdd[i]);  
	 scc.push_back(Odd[i]);  
	 sccErr.push_back(OddErr[i]);  
      }
   }

   return 0;
}

