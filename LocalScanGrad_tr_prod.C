// Get shimmed field gradient along the z direction
// NOTE: We use the TRLY here, not the PP! 
// for the special cases where we don't have PP data 

#include <cstdlib> 
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>  
#include <cmath> 

#include "TStyle.h"
#include "TLine.h"
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
#include "MovingAverage.h"
#include "Blinder.h"

#include "./include/Constants.h"
#include "./include/drift.h"
#include "./include/plungingProbeAnaEvent.h"
#include "./include/trolleyAnaEvent.h"
#include "./include/fixedProbeEvent.h"

#include "./src/CustomMath.C"
#include "./src/CustomImport.C"
#include "./src/CustomExport.C"
#include "./src/CustomGraph.C"
#include "./src/CustomAlgorithms.C"
#include "./src/CustomUtilities.C"
#include "./src/InputManager.C" 
#include "./src/BlindFuncs.C"
#include "./src/OscFuncs.C"
#include "./src/FitFuncs.C"
#include "./src/MyFits.C"
#include "./src/FitErr.C"
#include "./src/TRLYFuncs.C"
#include "./src/SystFuncs.C"

double fitFunc(double *x,double *p); 

int LocalScanGrad_tr_prod(std::string configFile){

   int rc=0;
   int method = gm2fieldUtil::Constants::kPhaseDerivative;

   std::cout << "------------------------------------" << std::endl;
   std::cout << "GET IMPOSED GRADIENT (USING TRLY Z SCANS)" << std::endl;

   InputManager *inputMgr = new InputManager();
   inputMgr->UseAxis();         // need to grab the axis data in the JSON file 
   inputMgr->Load(configFile);
   inputMgr->Print();

   std::string blindLabel    = inputMgr->GetBlindLabel();
   std::string prodVersion   = inputMgr->GetProductionTag();

   bool isBlind              = inputMgr->IsBlind();
   bool useOscCor            = inputMgr->GetOscCorStatus();
   int runPeriod             = inputMgr->GetRunPeriod();
   int nev                   = 20; // inputMgr->GetNumEventsToAvg(); 
   int probeNumber           = inputMgr->GetTrolleyProbe();

   // systematics 
   bool isSyst               = inputMgr->GetSystStatus();
   bool varyFit              = inputMgr->GetSystFitStatus("shim");  
   int systDirNum            = inputMgr->GetSystDirNum(); 

   std::string theDate = inputMgr->GetAnaDate();

   std::string plotDir = GetPath("plots" ,isBlind,blindLabel,theDate,isSyst,systDirNum);
   std::string outDir  = GetPath("output",isBlind,blindLabel,theDate,isSyst,systDirNum);

   std::vector<int> run;
   std::vector<std::string> label;
   inputMgr->GetRunList(run);
   inputMgr->GetRunLabels(label);

   // Trolley data 
   std::vector<trolleyAnaEvent_t> trlyData1;   
   rc = GetTrolleyData(run[0],method,trlyData1,prodVersion);
   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   // blinding 
   int blindUnits  = inputMgr->GetBlindUnits();
   double blindMag = inputMgr->GetBlindScale();
   gm2fieldUtil::Blinder *myBlind = new gm2fieldUtil::Blinder(blindLabel,blindMag,blindUnits);
   double blindValue = myBlind->GetBlinding(2); // in Hz

   if(isBlind) ApplyBlindingTRLY(blindValue,trlyData1);

   // fixed probe data 
   std::vector<int> fxprList;
   inputMgr->GetFXPRList(fxprList);

   std::vector<averageFixedProbeEvent_t> fxprData1;
   bool subtractDrift = inputMgr->GetFXPRDriftStatus();  
   int period = inputMgr->GetNumEventsToAvg(); 
   rc = GetFixedProbeData_avg(run[0],method,fxprList,fxprData1,prodVersion,subtractDrift,period,0);
   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   // reference time for oscillation correction
   // need first time of first TRLY event for the probe of interest  
   double t0 = trlyData1[0].time[probeNumber-1]/1E+9;
   std::cout << Form("REFERENCE TIME is %s",gm2fieldUtil::GetStringTimeStampFromUTC(t0).c_str()) << std::endl;

   std::vector<double> time1;
   time1.push_back(t0); 

   // apply oscillation correction 
   int NTR1 = trlyData1.size(); 
   // double lastTime = trlyData1[NTR1-1].time[probeNumber-1]/1E+9; 
   std::vector<double> trTime1,trPhi1,trFreq1,trFreq1_cor;
   // time.push_back(lastTime); // want all events so we use the last trolley event   
   rc = CorrectOscillation_trly(probeNumber-1,nev,time1,fxprData1,trlyData1,trTime1,trFreq1,trFreq1_cor);
   std::cout << "-------------" << std::endl;
   trFreq1.clear();
   trFreq1_cor.clear();
   rc = CorrectOscillation_trly(probeNumber-1,nev,time1,fxprData1,trlyData1,trPhi1 ,trFreq1,trFreq1_cor,"phi");
   std::cout << "-------------" << std::endl;

   TString xAxis = Form("phi");
   TString yAxis = Form("freq");

   std::vector<double> FREQ1,FREQ; 
   NTR1 = trFreq1.size();
   for(int i=0;i<NTR1;i++){
      if(useOscCor){
	 FREQ1.push_back(trFreq1_cor[i]); 
      }else{
	 FREQ1.push_back(trFreq1[i]); 
      }
      FREQ.push_back(trFreq1[i]); 
      std::cout << Form("%s: %.3lf, %.3lf",gm2fieldUtil::GetStringTimeStampFromUTC(trTime1[i]).c_str(),trPhi1[i],FREQ1[i]) << std::endl;
   }

   TGraph *g1 = gm2fieldUtil::Graph::GetTGraph(trPhi1,FREQ1);
   gm2fieldUtil::Graph::SetGraphParameters(g1,21,kBlack);

   TGraph *g1raw = GetTRLYTGraph(probeNumber-1,"phi","freq",trlyData1,1.0);
   gm2fieldUtil::Graph::SetGraphParameters(g1raw,20,kBlack);
   
   double angle   = inputMgr->GetTrolleyAngle(); // run 1 = 189.153, run 2 = 189.345  
   // std::cout << Form("angle = %.3lf deg, DeltaB = %.3lf Hz",angle,dB) << std::endl;  

   TMultiGraph *mg = new TMultiGraph();
   mg->Add(g1,"lp"); 

   TString fitName = "myFit";
   const int npar  = 4;  
   double delta    = 1.0; // roughly 124 mm 
   double min      = angle - delta;   
   double max      = angle + delta;
   std::cout << Form("Fit range: %.3lf to %.3lf",min,max) << std::endl;

   double par[npar]    = {1.,750.,1.,angle};
   double parErr[npar] = {0,0,0,0};  
   TF1 *myFit = new TF1(fitName,fitFunc,min,max,npar);
   for(int i=0;i<npar;i++) myFit->SetParameter(i,par[i]);  
   myFit->FixParameter(3,par[3]);   
   myFit->SetLineColor(kRed);

   TString Title    = Form("TRLY %02d Data, Run %d",probeNumber,run[0]);
   TString plotPath = Form("%s/tr-%02d_local-scan-z.png",plotDir.c_str(),probeNumber);

   TCanvas *c1 = new TCanvas("c1","TRLY Data",1200,600);
   c1->Divide(1,2);
   
   c1->cd(1);
   g1->Draw("ap");
   gm2fieldUtil::Graph::SetGraphLabels(g1,Title,"","Frequency (Hz)");
   gm2fieldUtil::Graph::SetGraphLabelSizes(g1,0.05,0.06); 
   g1->Draw("ap");
   TFitResultPtr fitResult = g1->Fit(fitName,"QS");
   c1->Update(); 

   c1->cd(2);
   g1raw->Draw("ap");
   gm2fieldUtil::Graph::SetGraphLabels(g1raw,"No Oscillation Correction","","Frequency (Hz)");
   gm2fieldUtil::Graph::SetGraphLabelSizes(g1raw,0.05,0.06); 
   g1raw->Draw("ap");
   c1->Update();

   c1->cd();
   c1->Print(plotPath); 

   for(int i=0;i<npar;i++){
      par[i]    = myFit->GetParameter(i); 
      parErr[i] = myFit->GetParError(i); 
      std::cout << Form("p[%d] = %.3lf ",i,par[i]) << std::endl;
   }

   std::cout << "FITTED RESULTS" << std::endl;
   double XX[3]    = {angle,0.,0.}; 
   double sf       = 1./124.; // NOTE: Not dividing by current anymore!  
   double dBdz     = sf*myFit->Derivative(XX[0]);  // evaluate at the trolley probe position of interest 
   double dBdz_err = sf*GetFitError(myFit,fitResult,MyPolyFitFuncDerivative_impGradZ,XX);

   if(isSyst && varyFit){
      std::cout << "[LocalScanGrad_tr_prod]: SYSTEMATIC VARIATION! Will vary fit result within (Gaussian) uncertainties" << std::endl;
      rc = systFunc::RandomizeFitValue(dBdz,dBdz_err);
   }
   std::cout << Form("dB/dz = %.3lf +/- %.3lf Hz/mm",dBdz,dBdz_err) << std::endl; 

   // print to file -- use same format as PP local scan code! 
   double PR[3] = {dBdz    ,dBdz    ,0}; 
   double ER[3] = {dBdz_err,dBdz_err,0};
   double DR[3] = {0,0,0};  

   std::string strTag = "azi-grad"; 
   char outpath[200];
   sprintf(outpath,"%s/%s_pr-%02d.csv",outDir.c_str(),strTag.c_str(),probeNumber); 
   PrintToFile(outpath,strTag,PR,ER,DR,DR);

   // reference plots 
   TGraph *g1r = gm2fieldUtil::Graph::GetTGraph(trTime1,trFreq1);
   gm2fieldUtil::Graph::SetGraphParameters(g1r,21,kBlack);

   TGraphErrors *gFXPR1 = GetFXPRTGraph_avg("GpsTimeStamp","freq","NONE",fxprData1);
   gm2fieldUtil::Graph::SetGraphParameters(gFXPR1,20,kBlack);

   TCanvas *c2 = new TCanvas("c2","Reference Plots",1200,600);
   c2->Divide(1,2);  

   c2->cd(1);
   g1r->Draw("alp");
   gm2fieldUtil::Graph::SetGraphLabels(g1r,Form("TRLY Data, Probe %02d, Run %d",probeNumber,runPeriod),"","Frequency (Hz)");
   gm2fieldUtil::Graph::SetGraphLabelSizes(g1r,0.05,0.06); 
   gm2fieldUtil::Graph::UseTimeDisplay(g1r);  
   g1r->Draw("alp");
   c2->Update(); 

   c2->cd(2);
   gFXPR1->Draw("alp");
   gm2fieldUtil::Graph::SetGraphLabels(gFXPR1,Form("FXPR Data, Run %d",run[0]),"","Frequency (Hz)");
   gm2fieldUtil::Graph::SetGraphLabelSizes(gFXPR1,0.05,0.06); 
   gm2fieldUtil::Graph::UseTimeDisplay(gFXPR1);  
   gFXPR1->Draw("alp");
   c2->Update(); 

   c2->cd();
   plotPath = Form("%s/azi-grad_pr-%02d.png",plotDir.c_str(),probeNumber);
   c2->Print(plotPath);  

   return 0;
}
//______________________________________________________________________________
double fitFunc(double *x,double *p){
   double X = x[0] - p[3]; 
   double f = p[0] + p[1]*X + p[2]*X*X;
   return f; 
}
