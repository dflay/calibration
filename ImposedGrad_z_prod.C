// Get imposed gradient along the z direction 

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
int PrintToFile_trlyDBz(const char *outpath,double dB,double dB_err); 

int ImposedGrad_z_prod(std::string configFile){

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
   int nev                   = 10; // inputMgr->GetNumEventsToAvg(); 
   int probeNumber           = inputMgr->GetTrolleyProbe();

   // systematics 
   bool isSyst               = inputMgr->GetSystStatus();
   bool varyFit              = inputMgr->GetSystFitStatus("imp-grad");  
   int systDirNum            = inputMgr->GetSystDirNum(); 

   // date_t theDate;
   // GetDate(theDate);
   std::string theDate = inputMgr->GetAnaDate();

   std::string plotDir = GetPath("plots" ,isBlind,blindLabel,theDate,isSyst,systDirNum);
   std::string outDir  = GetPath("output",isBlind,blindLabel,theDate,isSyst,systDirNum);

   // std::vector<int> run;
   // run.push_back(5455);  
   // run.push_back(5457);  

   std::vector<int> run;
   std::vector<std::string> label;
   inputMgr->GetRunList(run);
   inputMgr->GetRunLabels(label);

   // Trolley data 
   std::vector<trolleyAnaEvent_t> trlyData1,trlyData2;   
   rc = GetTrolleyData(run[0],method,trlyData1,prodVersion);
   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   // method = gm2fieldUtil::Constants::kHilbertPhaseLinear;  // this is a gradient field, so we use Hilbert  
   rc = GetTrolleyData(run[1],method,trlyData2,prodVersion);
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
   if(isBlind) ApplyBlindingTRLY(blindValue,trlyData2);

   // fixed probe data 
   std::vector<int> fxprList;
   inputMgr->GetFXPRList(fxprList);

   std::vector<averageFixedProbeEvent_t> fxprData1,fxprData2;
   bool subtractDrift = inputMgr->GetFXPRDriftStatus();  
   int period = inputMgr->GetNumEventsToAvg(); 
   rc = GetFixedProbeData_avg(run[0],method,fxprList,fxprData1,prodVersion,subtractDrift,period,0);
   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   rc = GetFixedProbeData_avg(run[1],method,fxprList,fxprData2,prodVersion,subtractDrift,period,0);
   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   // load time stamps 
   std::vector<double> time,time1,time2; 
   rc = LoadScanTimes(run[0],runPeriod,prodVersion,time1);
   rc = LoadScanTimes(run[1],runPeriod,prodVersion,time2);

   // std::cout << Form("run %05d TIMES",run[0]) << std::endl;
   // for(int i=0;i<time1.size();i++) std::cout << gm2fieldUtil::GetStringTimeStampFromUTC(time1[i]) << std::endl; 
   // std::cout << Form("run %05d TIMES",run[1]) << std::endl;
   // for(int i=0;i<time2.size();i++) std::cout << gm2fieldUtil::GetStringTimeStampFromUTC(time2[i]) << std::endl; 

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

   int NTR2 = trlyData2.size(); 
   // lastTime = trlyData2[NTR2-1].time[probeNumber-1]/1E+9;
   // time.clear(); 
   // time.push_back(lastTime); // want all events so we use the last trolley event   
   std::vector<double> trTime2,trPhi2,trFreq2,trFreq2_cor;
   rc = CorrectOscillation_trly(probeNumber-1,nev,time2,fxprData2,trlyData2,trTime2,trFreq2,trFreq2_cor);
   std::cout << "-------------" << std::endl;
   trFreq2.clear();
   trFreq2_cor.clear();
   rc = CorrectOscillation_trly(probeNumber-1,nev,time2,fxprData2,trlyData2,trPhi2 ,trFreq2,trFreq2_cor,"phi");
   std::cout << "-------------" << std::endl;

   TString xAxis = Form("phi");
   TString yAxis = Form("freq");

   std::vector<double> FREQ1,FREQ2; 
   NTR1 = trFreq1.size();
   for(int i=0;i<NTR1;i++){
      if(useOscCor){
	 FREQ1.push_back(trFreq1_cor[i]); 
      }else{
	 FREQ1.push_back(trFreq1[i]); 
      }
      std::cout << Form("%s: %.3lf, %.3lf",gm2fieldUtil::GetStringTimeStampFromUTC(trTime1[i]).c_str(),trPhi1[i],FREQ1[i]) << std::endl;
   }

   NTR2 = trFreq2.size();
   for(int i=0;i<NTR2;i++){
      // if(i%10==0) std::cout << gm2fieldUtil::GetStringTimeStampFromUTC(trTime2[i]) << std::endl;
      if(useOscCor){
	 FREQ2.push_back(trFreq2_cor[i]); 
      }else{
	 FREQ2.push_back(trFreq2[i]); 
      }
   }

   // this is the baseline 
   TGraph *g1 = gm2fieldUtil::Graph::GetTGraph(trPhi1,FREQ1);
   gm2fieldUtil::Graph::SetGraphParameters(g1,20,kBlack);

   // this is the gradient-on data
   TGraph *g2 = gm2fieldUtil::Graph::GetTGraph(trPhi2,FREQ2); 
   gm2fieldUtil::Graph::SetGraphParameters(g2,20,kRed);

   const int NPTS = 20;
   TGraph *gDiff = gm2fieldUtil::Graph::GetTGraphDifference(NPTS,g1,g2);
   gm2fieldUtil::Graph::SetGraphParameters(gDiff,20,kBlue);
   
   double angle   = inputMgr->GetTrolleyAngle(); // run 1 = 189.153, run 2 = 189.345  
   double dB      = gDiff->Eval(angle);
   double current = inputMgr->GetDBZCurrent(); // usually 0.82; 

   std::cout << Form("angle = %.3lf deg, DeltaB = %.3lf Hz",angle,dB) << std::endl;  

   TMultiGraph *mg = new TMultiGraph();
   mg->Add(g1,"lp"); 
   mg->Add(g2,"lp");

   TLegend *L = new TLegend(0.6,0.6,0.8,0.8);
   L->AddEntry(g1,Form("Run %d",run[0]),"lp"); 
   L->AddEntry(g2,Form("Run %d",run[1]),"lp");

   double yMin = -200; 
   double yMax = -15; 
   TLine *aLine = new TLine(angle,yMin,angle,yMax); 
   aLine->SetLineWidth(2);
   aLine->SetLineStyle(2);
   aLine->SetLineColor(kGreen+2);  

   TString fitName = "myFit";
   const int npar = 4;  
   double delta = 0.6; // 0.02; // roughly 2.5 mm 
   double min = angle - delta;   
   double max = angle + delta;
   std::cout << Form("Fit range: %.3lf to %.3lf",min,max) << std::endl;

   double par[npar] = {1.,750.,1.,angle};
   double parErr[npar] = {0,0,0,0};  
   TF1 *myFit = new TF1(fitName,fitFunc,min,max,npar);
   for(int i=0;i<npar;i++) myFit->SetParameter(i,par[i]);  
   myFit->FixParameter(3,par[3]);   
   myFit->SetLineColor(kRed);

   TString Title = Form("TRLY %02d Data: Black = Bare (Run %d), Red = Azi %.2lf A (Run %d)",probeNumber,run[0],current,run[1]);
   TString plotPath = Form("%s/tr-%02d_imposed-grad-z.png",plotDir.c_str(),probeNumber);

   TCanvas *c1 = new TCanvas("c1","TRLY Data",1200,600);
   c1->Divide(1,2);  
   
   c1->cd(1);
   mg->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mg,Title,"","Frequency (Hz)");
   gm2fieldUtil::Graph::SetGraphLabelSizes(mg,0.05,0.06); 
   mg->Draw("a");
   c1->Update();

   c1->cd(2);
   gStyle->SetOptFit(111); 
   gDiff->Draw("alp");
   gm2fieldUtil::Graph::SetGraphLabels(gDiff,"Azi - Bare","","Frequency (Hz)");
   gm2fieldUtil::Graph::SetGraphLabelSizes(gDiff,0.05,0.06); 
   gDiff->Draw("alp");
   TFitResultPtr fitResult = gDiff->Fit(fitName,"RQS");
   aLine->Draw("same"); 
   c1->Update();

   c1->cd();
   c1->Print(plotPath); 

   TF1 *theFit = gDiff->GetFunction(fitName); 
   for(int i=0;i<npar;i++){
      par[i]    = theFit->GetParameter(i); 
      parErr[i] = theFit->GetParError(i); 
      // std::cout << Form("p[%d] = %.3lf ",i,par[i]) << std::endl;
   }

   std::cout << "FITTED RESULTS" << std::endl;
   double XX[3] = {angle,0.,0.}; 
   double dB_fitted    = theFit->Eval(XX[0]); 
   double dB_fittedErr = parErr[0]; // GetFitError(myFit,fitResult,MyPolyFitFuncDerivative,XX); 
   std::cout << Form("dB = %.3lf +/- %.3lf Hz",dB_fitted,dB_fittedErr) << std::endl; 

   std::cout << "IMPOSED GRADIENT" << std::endl;
   // divide by current as well!  
   // double dBdz     = par[1]/124./current;     // convert from Hz/deg -> Hz/mm 
   // double dBdz_err = parErr[1]/124./current;  // convert from Hz/deg -> Hz/mm 
   // std::cout << Form("dB/dz = %.3lf +/- %.3lf Hz/mm",dBdz,dBdz_err) << std::endl; 

   double sf       = 1./124.; // NOTE: Not dividing by current anymore!  
   double dBdz     = sf*theFit->Derivative(XX[0]);  
   double dBdz_err = sf*GetFitError(theFit,fitResult,MyPolyFitFuncDerivative_impGradZ,XX); 

   if(isSyst && varyFit){
      std::cout << "[ImposedGrad_z_prod]: SYSTEMATIC VARIATION! Will vary fit result within (Gaussian) uncertainties" << std::endl;
      rc = systFunc::RandomizeFitValue(dBdz,dBdz_err);
   }
   std::cout << Form("dB/dz = %.3lf +/- %.3lf Hz/mm",dBdz,dBdz_err) << std::endl; 

   char outpath[200];
   sprintf(outpath,"%s/imposed-grad_z_pr-%02d.csv",outDir.c_str(),probeNumber); 
   rc = PrintToFile_trlyDBz(outpath,dBdz,dBdz_err);

   // reference plots 
   TGraph *g1r = gm2fieldUtil::Graph::GetTGraph(trTime1,trFreq1);
   gm2fieldUtil::Graph::SetGraphParameters(g1r,21,kBlack);

   TGraph *g1c = gm2fieldUtil::Graph::GetTGraph(trTime1,trFreq1_cor);
   gm2fieldUtil::Graph::SetGraphParameters(g1c,20,kRed);

   TGraph *g2r = gm2fieldUtil::Graph::GetTGraph(trTime2,trFreq2);
   gm2fieldUtil::Graph::SetGraphParameters(g2r,21,kBlack);

   TGraph *g2c = gm2fieldUtil::Graph::GetTGraph(trTime2,trFreq2_cor);
   gm2fieldUtil::Graph::SetGraphParameters(g2c,20,kRed);

   TMultiGraph *mg1 = new TMultiGraph(); 
   mg1->Add(g1r,"lp"); 
   mg1->Add(g1c,"lp"); 

   TMultiGraph *mg2 = new TMultiGraph(); 
   mg2->Add(g2r,"lp"); 
   mg2->Add(g2c,"lp"); 

   TString Title1 = Form("Run %d: black = raw, red = osc cor",run[0]);
   TString Title2 = Form("Run %d: black = raw, red = osc cor",run[1]);

   TCanvas *c2 = new TCanvas("c2","TRLY Reference Plots",1200,600);
   c2->Divide(1,2);  

   c2->cd(1);
   mg1->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mg1,Title1,"","Frequency (Hz)");
   gm2fieldUtil::Graph::SetGraphLabelSizes(mg1,0.05,0.06); 
   gm2fieldUtil::Graph::UseTimeDisplay(mg1);  
   mg1->Draw("a");
   c2->Update(); 

   c2->cd(2);
   mg2->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mg2,Title2,"","Frequency (Hz)");
   gm2fieldUtil::Graph::SetGraphLabelSizes(mg2,0.05,0.06);
   gm2fieldUtil::Graph::UseTimeDisplay(mg2);  
   mg2->Draw("a");
   c2->Update(); 

   c2->cd();
   plotPath = Form("%s/azi-grad_pr-%02d.png",plotDir.c_str(),probeNumber);
   c2->Print(plotPath);  

   TGraphErrors *gFXPR1 = GetFXPRTGraph_avg("GpsTimeStamp","freq","NONE",fxprData1);
   gm2fieldUtil::Graph::SetGraphParameters(gFXPR1,20,kBlack);

   TGraphErrors *gFXPR2 = GetFXPRTGraph_avg("GpsTimeStamp","freq","NONE",fxprData2);
   gm2fieldUtil::Graph::SetGraphParameters(gFXPR2,20,kBlack);

   TCanvas *c3 = new TCanvas("c3","FXPR Data",1200,600); 
   c3->Divide(1,2);
 
   c3->cd(1); 
   gFXPR1->Draw("alp");
   gm2fieldUtil::Graph::SetGraphLabels(gFXPR1,Form("FXPR Data, Run %d",run[0]),"","Frequency (Hz)");
   gm2fieldUtil::Graph::SetGraphLabelSizes(gFXPR1,0.05,0.06); 
   gm2fieldUtil::Graph::UseTimeDisplay(gFXPR1);  
   gFXPR1->Draw("alp");
   c3->Update(); 

   c3->cd(2); 
   gFXPR2->Draw("alp");
   gm2fieldUtil::Graph::SetGraphLabels(gFXPR2,Form("FXPR Data, Run %d",run[1]),"","Frequency (Hz)");
   gm2fieldUtil::Graph::SetGraphLabelSizes(gFXPR2,0.05,0.06); 
   gm2fieldUtil::Graph::UseTimeDisplay(gFXPR2);  
   gFXPR2->Draw("alp");
   c3->Update(); 

   c3->cd();
   plotPath = Form("%s/azi-grad_fxpr-data.png",plotDir.c_str());
   c3->Print(plotPath);  

   // now determine the gradient in the shimmed field -- will use as an estimate for some probes 
   double parBL[npar]    = {1.,750.,1.,angle};
   double parBLErr[npar] = {0,0,0,0};
   TF1 *myFitBL = new TF1("myFitBL",fitFunc,min,max,npar);
   for(int i=0;i<npar;i++) myFit->SetParameter(i,parBL[i]);  
   myFitBL->FixParameter(3,par[3]);   
   myFitBL->SetLineColor(kRed);

   TCanvas *c4 = new TCanvas("c4","Shimmed Field Plot",1200,600);

   c4->cd();
   g2->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(g2,"Shimmed Field","Azimuthal Position (deg)","Frequency (Hz)");
   g2->Draw("a");
   TFitResultPtr fitResultBL = g2->Fit("myFitBL","RQS");
   c4->Update(); 

   plotPath = Form("%s/azi-grad_est_pr-%02d.png",plotDir.c_str(),probeNumber);
   c4->Print(plotPath);  

   TF1 *theFitBL = g2->GetFunction("myFitBL"); 
   for(int i=0;i<npar;i++){
      parBL[i]    = theFitBL->GetParameter(i); 
      parBLErr[i] = theFitBL->GetParError(i); 
      // std::cout << Form("p[%d] = %.3lf ",i,par[i]) << std::endl;
   }

   double dBdz_bl     = sf*theFitBL->Derivative(XX[0]);  
   double dBdz_bl_err = sf*GetFitError(theFitBL,fitResultBL,MyPolyFitFuncDerivative_impGradZ,XX);  

   if(isSyst && varyFit){
      // std::cout << "[ImposedGrad_z_prod]: SYSTEMATIC VARIATION! Will vary fit result within (Gaussian) uncertainties" << std::endl;
      rc = systFunc::RandomizeFitValue(dBdz_bl,dBdz_bl_err);
   }

   // assume 100% uncertainty on probes 1, 3, 9, 11 
   if(probeNumber==1||probeNumber==3||probeNumber==9||probeNumber==11) dBdz_bl_err = TMath::Abs(dBdz_bl); 

   sprintf(outpath,"%s/azi-grad_est_pr-%02d.csv",outDir.c_str(),probeNumber); 
   rc = PrintToFile_trlyDBz(outpath,dBdz_bl,dBdz_bl_err);

   return 0;
}
//______________________________________________________________________________
double fitFunc(double *x,double *p){
   double X = x[0] - p[3]; 
   double f = p[0] + p[1]*X + p[2]*X*X;
   return f; 
}
//______________________________________________________________________________
int PrintToFile_trlyDBz(const char *outpath,double dB,double dB_err){
   char outStr[200]; 

   std::ofstream outfile;
   outfile.open(outpath);
   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      outfile << "#grad(Hz/mm),grad-err(Hz/mm)"  << std::endl;
      sprintf(outStr,"%.3lf,%.3lf",dB,dB_err);
      outfile << outStr << std::endl;
      outfile.close();
      sprintf(outStr,"Data printed to file: %s",outpath); 
      PrintMessage("ImposedGrad_z_prod",outStr); 
   }
   return 0;
}

