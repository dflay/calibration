// Compare PP and TRLY measurements across multiple swaps back and forth  
// TODO 
// - Implement scheme to bring in imposed gradients (from input files).  
//    Fit to function (pol2) and evaluate at the correct (r,y,z) of trolley probe.  
//    Don't have azimuthal done yet, so hold off there 
// - Combine shot errors and misalignment errors for final uncertainty  

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
#include "TText.h"

#include "RootTreeStructs.h"
#include "gm2fieldMath.h"
#include "gm2fieldGraph.h"
#include "gm2fieldRootHelper.h"
#include "gm2fieldUnits.h"
#include "TemperatureSensor.h"
#include "MovingAverage.h"

#include "./include/date.h"
#include "./include/Constants.h"
#include "./include/fixedProbeEvent.h"
#include "./include/trolleyAnaEvent.h"
#include "./include/nmr_meas.h"
#include "./include/perturbation.h" 

#include "./src/InputManager.C"
#include "./src/FitFuncs.C"
#include "./src/FXPRFuncs.C"
#include "./src/TRLYFuncs.C"
#include "./src/CalibFuncs.C"
#include "./src/Consolidate.C"
#include "./src/CustomUtilities.C"
#include "./src/CustomMath.C"
#include "./src/CustomGraph.C"
#include "./src/CustomImport.C"
#include "./src/CustomExport.C"
#include "./src/CustomAlgorithms.C"

double gMarkerSize = 0.8;
double gXSize      = 0.05;  
double gYSize      = 0.06;  

int GetStats(std::vector<double> x,std::vector<double> dx,double &mean,double &err,double &stdev);
int PrintResults(std::string outpath,std::vector<double> trial,
                 std::vector<double> x1,std::vector<double> dx1,
                 std::vector<double> x2,std::vector<double> dx2,
                 std::vector<double> x3,std::vector<double> dx3,
                 std::vector<double> x4,std::vector<double> dx4);

int LoadCalibData(const char *inpath,std::vector<int> &run,
                  std::vector<double> &time,std::vector<double> &temp,
                  std::vector<double> &x,std::vector<double> &dx,
                  std::vector<double> &y,std::vector<double> &dy);

int CalibrateMultiSwap(std::string configFile){

   int rc=0;
   int method = gm2fieldUtil::Constants::kPhaseDerivative;

   InputManager *inputMgr = new InputManager();
   inputMgr->Load(configFile);
   inputMgr->Print();

   std::string anaDate = inputMgr->GetAnalysisDate();
   bool isBlind        = inputMgr->IsBlind();
   int probe           = inputMgr->GetTrolleyProbe(); 

   date_t theDate;
   GetDate(theDate);

   char plotDir[200];
   sprintf(plotDir,"./plots/%02d-%02d-%02d",theDate.month,theDate.day,theDate.year-2000);
   rc = MakeDirectory(plotDir);

   char outDir[200];
   sprintf(outDir,"./output"); 
   if(isBlind)  sprintf(outDir,"%s/blinded"  ,outDir);
   if(!isBlind) sprintf(outDir,"%s/unblinded",outDir);
   sprintf(outDir,"%s/%02d-%02d-%02d",outDir,theDate.month,theDate.day,theDate.year-2000); 
   rc = MakeDirectory(outDir);

   char outPath[500]; 
   sprintf(outPath,"%s/results_%s.csv",outDir,anaDate.c_str());
   std::string outpath = outPath;  

   // load the PP data 
   std::vector<int> ppRun;
   std::vector<double> ppTime,ppTemp,ppFreq,ppFreqErr,ppFreqCor,ppFreqCorErr;
   char dataPath[500]; 
   sprintf(dataPath,"%s/pp-data_%s.csv",outDir,anaDate.c_str());
   rc = LoadCalibData(dataPath,ppRun,ppTime,ppTemp,ppFreq,ppFreqErr,ppFreqCor,ppFreqCorErr); 
   if(rc!=0){
      std::cout << "No plunging probe data!" << std::endl;
      return 1;
   } 

   // load the trolley data
   std::vector<int> trRun;
   std::vector<double> trTime,trTemp,trFreq,trFreqErr,trFreqCor,trFreqCorErr;
   sprintf(dataPath,"%s/trly-data_pr-%02d_%s.csv",outDir,probe,anaDate.c_str());
   rc = LoadCalibData(dataPath,trRun,trTime,trTemp,trFreq,trFreqErr,trFreqCor,trFreqCorErr);
   if(rc!=0){
      std::cout << "No trolley data!" << std::endl;
      return 1;
   } 

   int NPP = ppRun.size();
   int NTR = trRun.size();
   if(NPP!=NTR){
      std::cout << "Number of runs for PP and TRLY do not match!  Exiting..." << std::endl;
      std::cout << "PP runs:   " << NPP << std::endl;
      std::cout << "TRLY runs: " << NTR << std::endl;
      return 1;
   }

   // load in perturbation data
   perturbation_t ppPert;
   char inpath[500]; 
   sprintf(inpath,"./input/perturbation/pp-pert.csv");
   LoadPerturbationData(inpath,ppPert);

   // put PP data into nmr_meas_t struct 
   nmr_meas_t data;
   std::vector<nmr_meas_t> ppData; 
   for(int i=0;i<NPP;i++){
      data.freq          = ppFreq[i];
      data.freq_err      = ppFreqErr[i];
      data.freq_fxpr     = ppFreqCor[i];
      data.freq_fxpr_err = ppFreqCorErr[i];
      data.freq_trly     = 0; 
      data.freq_trly_err = 0;
      ppData.push_back(data);  
   }

   double freqFree[3],freqFreeErr[3];
   std::vector<double> ppFreqFree,ppFreqFreeErr,ppFreqFreeCor,ppFreqFreeCorErr;

   // compute omega_p free (PP) 
   for(int i=0;i<NPP;i++){
      GetOmegaP_free(ppData[i],ppPert,freqFree,freqFreeErr);
      ppFreqFree.push_back(freqFree[0]); 
      ppFreqFreeErr.push_back(freqFreeErr[0]); 
      ppFreqFreeCor.push_back(freqFree[1]); 
      ppFreqFreeCorErr.push_back(freqFreeErr[1]);
      std::cout << Form("%d: %.3lf, %.3lf",i,ppFreq[i],ppFreqFree[i]) << std::endl; 
   }

   double diff=0,err=0,diff_prev=0,err_prev=0,diff_cor=0,err_cor=0,diff_cor_w=0,err_cor_w=0;
   double diff_free=0,err_free=0,diff_free_prev=0,err_free_prev=0,diff_free_cor=0,err_free_cor=0,diff_free_cor_w=0,err_free_cor_w=0;
   std::vector<double> D,PPR,TRR;
   std::vector<double> trial,DIFF,ERR,DIFF_cor,ERR_cor,DIFF_cor_w,ERR_cor_w; 
   std::vector<double> DIFF_FREE,ERR_FREE,DIFF_FREE_cor,ERR_FREE_cor,DIFF_FREE_cor_w,ERR_FREE_cor_w;

   for(int i=0;i<NPP;i++) PPR.push_back(ppRun[i]); 
   for(int i=0;i<NTR;i++) TRR.push_back(trRun[i]); 

   std::cout << "Processing " << NPP << " trials..." << std::endl;

   // get the raw difference 
   double arg=0,argErr=0;
   std::vector<double> rT,rD,rD_err; 
   for(int i=0;i<NPP;i++){
      arg = ppFreq[i] - trFreq[i];  
      argErr = TMath::Sqrt( ppFreqErr[i]*ppFreqErr[i] + trFreqErr[i]*trFreqErr[i] ); 
      rT.push_back(i+1); 
      rD.push_back(arg); 
      rD_err.push_back(argErr); 
   }

   // create a weight for each trial based on time of the measurement
   double arg_time=0,arg_time_prev=0; 
   std::vector<double> timeDiffPrev,timeDiff,timeWeight,timeAllPairs;  

   for(int i=1;i<NPP;i++){
      // straight difference on each PP, TRLY data point 
      arg           = TMath::Abs(ppTime[i-1]-trTime[i-1]);  
      timeDiffPrev.push_back(arg);
      timeAllPairs.push_back(arg);
      // look at i PP vs i-1 TRLY 
      arg           = TMath::Abs(ppTime[i]-trTime[i-1]);  
      timeDiff.push_back(arg);
      timeAllPairs.push_back(arg);
   }

   // straight difference on each PP, TRLY data point 
   arg           = TMath::Abs(ppTime[NPP-1]-trTime[NPP-1]);  
   timeDiffPrev.push_back(arg);

   double w_prev=0,w=0;
   double dt_tot  = 0;
   double dt_mean = gm2fieldUtil::Math::GetMean<double>(timeAllPairs); 
   double dt_sig  = gm2fieldUtil::Math::GetStandardDeviation<double>(timeAllPairs); 

   double dt=0,dt_prev=0;

   // now do the ABA analysis 
   for(int i=1;i<NPP;i++){
      std::cout << "========================= TRIAL " << i << " =========================" << std::endl;
      // compute the time weight
      arg_time_prev = TMath::Abs(ppTime[i-1] - trTime[i-1]); 
      arg_time      = TMath::Abs(ppTime[i] - trTime[i-1]);
      dt_tot        = TMath::Abs(ppTime[i]-ppTime[i-1]); 
      dt_prev       = TMath::Abs(arg_time_prev-dt_mean);  
      dt            = TMath::Abs(arg_time-dt_mean);  
      // if(dt_prev>dt){
	 // prev time is bad 
	 w_prev = arg_time/dt_tot;
	 w      = arg_time_prev/dt_tot;
      // }else{
      //    w_prev = arg_time_prev/dt_tot;
      //    w      = arg_time/dt_tot;
      // } 
      // compute difference and error 
      diff_prev = ppFreq[i-1] - trFreq[i-1]; 
      err_prev  = TMath::Sqrt( ppFreqErr[i-1]*ppFreqErr[i-1] + trFreqErr[i-1]*trFreqErr[i-1] ); 
      // order of measurements was PP, then TRLY.  So now consider the next difference, looking at PP(i) vs TRLY(i-1) 
      diff      = ppFreq[i] - trFreq[i-1]; 
      err       = TMath::Sqrt( ppFreqErr[i]*ppFreqErr[i] + trFreqErr[i-1]*trFreqErr[i-1] ); 
      // take average of these differences 
      diff_cor  = 0.5*(diff_prev + diff);  
      err_cor   = 0.5*(err_prev + err);
      // take average of these differences -- time weighted  
      diff_cor_w  = w_prev*diff_prev + w*diff;  
      err_cor_w   = w_prev*err_prev  + w*err;
      // now for the free proton version  
      // compute difference and error 
      diff_free_prev = ppFreqFree[i-1] - trFreq[i-1]; 
      err_free_prev  = TMath::Sqrt( ppFreqFreeErr[i-1]*ppFreqFreeErr[i-1] + trFreqErr[i-1]*trFreqErr[i-1] ); 
      // order of measurements was PP, then TRLY.  So now consider the next difference, looking at PP(i) vs TRLY(i-1) 
      diff_free      = ppFreqFree[i] - trFreq[i-1]; 
      err_free       = TMath::Sqrt( ppFreqFreeErr[i]*ppFreqFreeErr[i] + trFreqErr[i-1]*trFreqErr[i-1] ); 
      // take average of these differences 
      diff_free_cor  = 0.5*(diff_free_prev + diff_free);  
      err_free_cor   = 0.5*(err_free_prev + err_free); 
      // time-weighted  
      diff_free_cor_w = w_prev*diff_free_prev + w*diff_free;  
      err_free_cor_w  = w_prev*err_free_prev  + w*err_free; 
      // now store the result (no correction) 
      trial.push_back(i); 
      DIFF.push_back(diff_prev);
      ERR.push_back(err_prev);  
      DIFF_cor.push_back(diff_cor);
      ERR_cor.push_back(err_cor); 
      DIFF_cor_w.push_back(diff_cor_w);
      ERR_cor_w.push_back(err_cor_w); 
      DIFF_FREE.push_back(diff_free_prev);
      ERR_FREE.push_back(err_free_prev);  
      DIFF_FREE_cor.push_back(diff_free_cor);
      ERR_FREE_cor.push_back(err_free_cor); 
      DIFF_FREE_cor_w.push_back(diff_free_cor_w);
      ERR_FREE_cor_w.push_back(err_free_cor_w); 
      // print the results 
      std::cout << "Frequencies: " << std::endl;
      std::cout << Form("PP [A]: %.3lf +/- %.3lf, TRLY [B]: %.3lf +/- %.3lf, PP [A]: %.3lf +/- %.3lf",
                        ppFreq[i-1],ppFreqErr[i-1],trFreq[i-1],trFreqErr[i-1],ppFreq[i],ppFreqErr[i]) << std::endl;
      std::cout << "Differences: " << std::endl;
      std::cout << Form("diff1 = %.3lf +/- %.3lf, diff2 = %.3lf +/- %.3lf",diff_prev,err_prev,diff,err) << std::endl; 
      std::cout << Form("result: PP - TRLY = %.3lf +/- %.3lf, drift-cor = %.3lf +/- %.3lf, drift-cor-w = %.3lf +/- %.3lf",
                        diff_prev,err_prev,diff_cor,err_cor,diff_cor_w,err_cor_w) << std::endl;
      std::cout << "------------------------------" << std::endl; 
      std::cout << "Using the free proton PP data: " << std::endl;
      std::cout << Form("PP [A]: %.3lf +/- %.3lf, TRLY [B]: %.3lf +/- %.3lf, PP [A]: %.3lf +/- %.3lf",
                        ppFreqFree[i-1],ppFreqFreeErr[i-1],trFreq[i-1],trFreqErr[i-1],ppFreqFree[i],ppFreqFreeErr[i]) << std::endl;
      std::cout << "Differences: " << std::endl;
      std::cout << Form("diff1 = %.3lf +/- %.3lf, diff2 = %.3lf +/- %.3lf",diff_free_prev,err_free_prev,diff_free,err_free) << std::endl; 
      std::cout << Form("result: PP - TRLY = %.3lf +/- %.3lf, drift-cor = %.3lf +/- %.3lf, drift-cor-w = %.3lf +/- %.3lf",
                        diff_free_prev,err_free_prev,diff_free_cor,err_free_cor,diff_free_cor_w,err_free_cor_w) << std::endl;
   }

   std::cout << "=====================================" << std::endl; 

   // get mean and stdev of the runs
   double mean_diff=0,stdev_diff=0,err_diff=0; 
   double mean_diff_cor=0,stdev_diff_cor=0,err_diff_cor=0; 
   double mean_diff_cor_w=0,stdev_diff_cor_w=0,err_diff_cor_w=0; 
   rc = GetStats(DIFF    ,ERR    ,mean_diff    ,err_diff    ,stdev_diff); 
   rc = GetStats(DIFF_cor,ERR_cor,mean_diff_cor,err_diff_cor,stdev_diff_cor); 
   rc = GetStats(DIFF_cor_w,ERR_cor_w,mean_diff_cor_w,err_diff_cor_w,stdev_diff_cor_w); 

   double yMin = 0.5*(mean_diff + mean_diff_cor) - 7.;
   double yMax = 0.5*(mean_diff + mean_diff_cor) + 5.;

   // get mean and stdev of the runs (free)
   double mean_diff_free=0    ,stdev_diff_free=0    ,err_diff_free=0; 
   double mean_diff_free_cor=0,stdev_diff_free_cor=0,err_diff_free_cor=0; 
   double mean_diff_free_cor_w=0,stdev_diff_free_cor_w=0,err_diff_free_cor_w=0; 
   rc = GetStats(DIFF_FREE    ,ERR_FREE    ,mean_diff_free    ,err_diff_free    ,stdev_diff_free); 
   rc = GetStats(DIFF_FREE_cor,ERR_FREE_cor,mean_diff_free_cor,err_diff_free_cor,stdev_diff_free_cor); 
   rc = GetStats(DIFF_FREE_cor_w,ERR_FREE_cor_w,mean_diff_free_cor_w,err_diff_free_cor_w,stdev_diff_free_cor_w); 

   double yMinFree = 0.5*(mean_diff_free + mean_diff_free_cor) - 7.;
   double yMaxFree = 0.5*(mean_diff_free + mean_diff_free_cor) + 5.;

   std::cout << "MEAN" << std::endl;
   std::cout << Form("[raw]:     %.3lf +/- %.3lf +/- %.3lf Hz",mean_diff    ,err_diff,stdev_diff)     << std::endl; 
   std::cout << Form("[ABA]:     %.3lf +/- %.3lf +/- %.3lf Hz",mean_diff_cor,err_diff_cor,stdev_diff_cor) << std::endl; 
   std::cout << Form("[ABA] [w]: %.3lf +/- %.3lf +/- %.3lf Hz",mean_diff_cor_w,err_diff_cor_w,stdev_diff_cor_w) << std::endl; 

   std::cout << "MEAN (free)" << std::endl;
   std::cout << Form("[raw]:     %.3lf +/- %.3lf +/- %.3lf Hz",mean_diff_free    ,err_diff_free,stdev_diff_free)     << std::endl; 
   std::cout << Form("[ABA]:     %.3lf +/- %.3lf +/- %.3lf Hz",mean_diff_free_cor,err_diff_free_cor,stdev_diff_free_cor) << std::endl; 
   std::cout << Form("[ABA] [w]: %.3lf +/- %.3lf +/- %.3lf Hz",mean_diff_free_cor_w,err_diff_free_cor_w,stdev_diff_free_cor_w) << std::endl; 

   // make a special verison of the PP, TRLY plots with a better scale 
   std::vector<double> F,FPP,FTR; 
   for(int i=0;i<NPP;i++) F.push_back(ppFreq[i]);
   for(int i=0;i<NPP;i++) F.push_back(trFreq[i]);
   double mean_freq = gm2fieldUtil::Math::GetMean<double>(F);
   for(int i=0;i<NPP;i++){
      FPP.push_back( ppFreq[i]-mean_freq ); 
      FTR.push_back( trFreq[i]-mean_freq ); 
   } 

   // make plots 
   TGraphErrors *gPP             = gm2fieldUtil::Graph::GetTGraphErrors(ppTime,FPP,ppFreqErr);
   TGraphErrors *gTR             = gm2fieldUtil::Graph::GetTGraphErrors(trTime,FTR,trFreqErr);
   TGraphErrors *gDiff           = gm2fieldUtil::Graph::GetTGraphErrors(trial,DIFF,ERR);
   TGraphErrors *gDiff_cor       = gm2fieldUtil::Graph::GetTGraphErrors(trial,DIFF_cor,ERR_cor);
   TGraphErrors *gDiff_cor_w     = gm2fieldUtil::Graph::GetTGraphErrors(trial,DIFF_cor_w,ERR_cor_w);
   TGraphErrors *gDiffFree       = gm2fieldUtil::Graph::GetTGraphErrors(trial,DIFF_FREE,ERR_FREE);
   TGraphErrors *gDiffFree_cor   = gm2fieldUtil::Graph::GetTGraphErrors(trial,DIFF_FREE_cor,ERR_FREE_cor);
   TGraphErrors *gDiffFree_cor_w = gm2fieldUtil::Graph::GetTGraphErrors(trial,DIFF_FREE_cor_w,ERR_FREE_cor_w);
   TGraphErrors *gRaw            = gm2fieldUtil::Graph::GetTGraphErrors(rT,rD,rD_err);
   TGraph       *gTimeDiffPrev   = gm2fieldUtil::Graph::GetTGraph(rT,timeDiffPrev);  
   TGraph       *gTimeDiff       = gm2fieldUtil::Graph::GetTGraph(trial,timeDiff);  

   gm2fieldUtil::Graph::SetGraphParameters(gPP          ,20,kBlue );
   gm2fieldUtil::Graph::SetGraphParameters(gTR          ,20,kRed  );
   gm2fieldUtil::Graph::SetGraphParameters(gDiff        ,20,kBlack);
   gm2fieldUtil::Graph::SetGraphParameters(gDiff_cor    ,20,kRed  );
   gm2fieldUtil::Graph::SetGraphParameters(gDiff_cor_w  ,20,kBlue );
   gm2fieldUtil::Graph::SetGraphParameters(gDiffFree    ,20,kBlack);
   gm2fieldUtil::Graph::SetGraphParameters(gDiffFree_cor,20,kRed  );
   gm2fieldUtil::Graph::SetGraphParameters(gDiffFree_cor_w,20,kBlue );
   gm2fieldUtil::Graph::SetGraphParameters(gRaw         ,20,kBlack);
   gm2fieldUtil::Graph::SetGraphParameters(gTimeDiffPrev,20,kMagenta);
   gm2fieldUtil::Graph::SetGraphParameters(gTimeDiff    ,20,kBlack);

   // now append the averaged results to the individual data sets 
   trial.push_back(99); // not a real trial 
   DIFF.push_back(mean_diff); 
   ERR.push_back(stdev_diff); 
   DIFF_cor.push_back(mean_diff_cor); 
   ERR_cor.push_back(stdev_diff_cor); 
   DIFF_cor_w.push_back(mean_diff_cor_w); 
   ERR_cor_w.push_back(stdev_diff_cor_w); 
   DIFF_FREE.push_back(mean_diff_free); 
   ERR_FREE.push_back(stdev_diff_free); 
   DIFF_FREE_cor.push_back(mean_diff_free_cor); 
   ERR_FREE_cor.push_back(stdev_diff_free_cor); 
   DIFF_FREE_cor_w.push_back(mean_diff_free_cor_w); 
   ERR_FREE_cor_w.push_back(stdev_diff_free_cor_w); 

   // print to file 
   rc = PrintResults(outpath,trial,DIFF,ERR,DIFF_cor,ERR_cor,DIFF_FREE,ERR_FREE,DIFF_FREE_cor,ERR_FREE_cor); 

   TMultiGraph *mgp = new TMultiGraph(); 
   mgp->Add(gPP,"lp"); 
   mgp->Add(gTR,"lp"); 

   TMultiGraph *mg = new TMultiGraph();
   mg->Add(gDiff_cor,"lp");
   mg->Add(gDiff_cor_w,"lp");

   TMultiGraph *mgf = new TMultiGraph();
   mgf->Add(gDiffFree_cor,"lp");
   mgf->Add(gDiffFree_cor_w,"lp");

   TMultiGraph *mgt = new TMultiGraph();
   mgt->Add(gTimeDiffPrev,"lp"); 
   mgt->Add(gTimeDiff    ,"lp"); 

   TLegend *L = new TLegend(0.6,0.6,0.8,0.8);
   L->AddEntry(gDiff    ,"Raw"       ,"p");    
   L->AddEntry(gDiff_cor,"ABA Method","p");   
 
   TLegend *LP = new TLegend(0.6,0.6,0.8,0.8);
   LP->AddEntry(gPP,"Plunging Probe","p");    
   LP->AddEntry(gTR,"Trolley"       ,"p");

   TCanvas *c1 = new TCanvas("c1","PP, TRLY plots",1200,800);
   c1->Divide(1,3);

   c1->cd(1);
   mgp->Draw("a"); 
   gm2fieldUtil::Graph::SetGraphLabels(mgp,"PP (Blue), TRLY (Red) Data","",Form("Frequency-%.3lf MHz (Hz)",mean_freq/1E+6));
   gm2fieldUtil::Graph::UseTimeDisplay(mgp);
   gm2fieldUtil::Graph::SetGraphLabelSizes(mgp,gXSize,gYSize);  
   mgp->Draw("a");
   // LP->Draw("same"); 
   c1->Update(); 

   c1->cd(2);
   gRaw->Draw("ap"); 
   gm2fieldUtil::Graph::SetGraphLabels(gRaw,"PP-TRLY Data","Trial","Frequency Difference (Hz)");
   gm2fieldUtil::Graph::SetGraphLabelSizes(gRaw,gXSize,gYSize);  
   gRaw->Draw("ap");
   c1->Update(); 

   c1->cd(3);
   mgt->Draw("ap"); 
   gm2fieldUtil::Graph::SetGraphLabels(mgt,"Cycle Time","Trial","Cycle Time (sec)");
   gm2fieldUtil::Graph::SetGraphLabelSizes(gTimeDiff,gXSize,gYSize);  
   mgt->Draw("ap");
   c1->Update(); 

   // save the plot  
   c1->cd(); 
   TString calibPlots = Form("%s/pp-trly-data_%s.png",plotDir,anaDate.c_str());
   c1->Print(calibPlots);

   TString statsText = Form("#mu = %.3lf Hz, #sigma = %.3lf Hz",mean_diff_cor,stdev_diff_cor);
   TText *myText = new TText(0.6,0.6,statsText); 
   myText->SetTextFont(43); 
   myText->SetTextSize(40); 

   TCanvas *c2 = new TCanvas("c2","Calib Data",1200,600);
   c2->cd();

   mg->Draw("a"); 
   gm2fieldUtil::Graph::SetGraphLabels(mg,"PP-TRLY (red = ABA Method, blue = time-weighted ABA)","Trial","Frequency Difference (Hz)"); 
   mg->GetYaxis()->SetRangeUser(yMin,yMax);
   mg->Draw("a");
   myText->Draw("same"); 
   c2->Update(); 

   // save the plot  
   c2->cd(); 
   calibPlots = Form("%s/calib_%s.png",plotDir,anaDate.c_str());
   c2->Print(calibPlots);

   TCanvas *c3 = new TCanvas("c3","Calib Data (Free Proton)",1200,600);
   c3->cd();

   mgf->Draw("a"); 
   gm2fieldUtil::Graph::SetGraphLabels(mgf,"PP (free) - TRLY (red = ABA Method, blue = time-weighted ABA)","Trial","Frequency Difference (Hz)"); 
   mgf->GetYaxis()->SetRangeUser(yMinFree,yMaxFree);
   mgf->Draw("a");
   c3->Update(); 

   // save the plot  
   c3->cd(); 
   calibPlots = Form("%s/calib-free_%s.png",plotDir,anaDate.c_str());
   c3->Print(calibPlots);

   return 0;
}
//______________________________________________________________________________
int GetStats(std::vector<double> x,std::vector<double> dx,double &mean,double &err,double &stdev){
   // get weighted mean, stdev of a data set x +/- dx
   double arg=0;
   const int N = x.size();
   std::vector<double> weight; 
   for(int i=0;i<N;i++){
      arg = 1./( dx[i]*dx[i] );
      weight.push_back(arg); 
   }

   int rc = gm2fieldUtil::Math::GetWeightedMean<double>(x,weight,mean,err); 
   stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(x);

   return 0; 
}
//______________________________________________________________________________
int LoadCalibData(const char *inpath,std::vector<int> &run,std::vector<double> &time,std::vector<double> &temp,
                  std::vector<double> &x,std::vector<double> &dx,
                  std::vector<double> &y,std::vector<double> &dy){

   int ir; 
   double itime,itemp,ix,idx,iy,idy;
   std::string sr,stime,stemp,sx,sdx,sy,sdy;

   std::ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      std::cout << "Cannot open the file: " << inpath << std::endl;
      return 1;
   }else{ 
      while( !infile.eof() ){
         std::getline(infile,sr   ,',');
         std::getline(infile,stime,',');
         std::getline(infile,stemp,',');
         std::getline(infile,sx   ,',');
         std::getline(infile,sdx  ,',');
         std::getline(infile,sy   ,',');
         std::getline(infile,sdy      );
         ir    = std::atof( sr.c_str()  ); 
         itime = std::atof( stime.c_str()  ); 
         itemp = std::atof( stemp.c_str()  ); 
         ix    = std::atof( sx.c_str()  );
         idx   = std::atof( sdx.c_str() );
         iy    = std::atof( sy.c_str()  );
         idy   = std::atof( sdy.c_str() );
         run.push_back(ir);
         time.push_back(itime);  
         temp.push_back(itemp);  
         x.push_back(ix);
         dx.push_back(idx);
         y.push_back(iy);
         dy.push_back(idy);
      }
      infile.close();
      run.pop_back(); 
      x.pop_back();  
      time.pop_back();  
      temp.pop_back();  
      dx.pop_back(); 
      y.pop_back(); 
      dy.pop_back(); 
   }
   return 0;
}
//______________________________________________________________________________
int PrintResults(std::string outpath,std::vector<double> trial,
                 std::vector<double> x1,std::vector<double> dx1,
                 std::vector<double> x2,std::vector<double> dx2,
                 std::vector<double> x3,std::vector<double> dx3,
                 std::vector<double> x4,std::vector<double> dx4){

   const int N = trial.size(); 
   char myStr[1000],header[500]; 
   // prepare the header 
   sprintf(header,"#diff,err,diff_cor,err_cor,diff_free,err_free,diff_free_cor,err_free_cor");  

   std::ofstream outfile;
   outfile.open(outpath.c_str());
   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      outfile << header << std::endl; 
      for(int i=0;i<N;i++){
	 sprintf(myStr,"%02d,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf",
                 (int)trial[i],x1[i],dx1[i],x2[i],dx2[i],x3[i],dx3[i],x4[i],dx4[i]);
	 outfile << myStr << std::endl;
      }
      std::cout << "The data has been written to file: " << outpath << std::endl;
   }  
   return 0;
}
