// Determine the Delta-B values when toggling SCC currents on and off   

#include <cstdlib> 
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>  
#include <cmath> 

#include "TH2D.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPad.h"
#include "TSpline.h"

#include "RootTreeStructs.h"
#include "gm2fieldMath.h"
#include "gm2fieldGraph.h"
#include "gm2fieldRootHelper.h"
#include "gm2fieldUnits.h"
#include "TemperatureSensor.h"

#include "./include/Constants.h"
#include "./include/fixedProbeEvent.h"
#include "./include/trolleyAnaEvent.h"
#include "./include/sccEvent.h" 

#include "./src/FitFuncs.C"
#include "./src/FXPRFuncs.C"
#include "./src/TRLYFuncs.C"
#include "./src/Consolidate.C"
#include "./src/CustomMath.C"
#include "./src/CustomGraph.C"
#include "./src/CustomImport.C"
#include "./src/CustomExport.C"
#include "./src/CustomAlgorithms.C"

TGraph *GetSCCPlot(int type,std::vector<gm2field::surfaceCoils_t> data); 

int PrintToFile_avg(const char *outpath,int probe,double scc_on,double scc_on_err,double bare,double bare_err);
int PrintToFile_avg(const char *outpath,int probe,double diff,double diff_err,double diff_aba,double diff_aba_err);
int PrintToFile_scc(const char *outpath,int probe,int trial,double scc_on,double scc_on_err,double bare,double bare_err); 

int CombineResults(std::vector<double> lo,std::vector<double> hi,std::vector<double> &u);
int GetStats(int probe,int nev,std::vector<double> time,std::vector<trolleyAnaEvent_t> Data,
             std::vector<double> &MEAN,std::vector<double> &STDEV);

int CalculateDiff(std::vector<double> mu1,std::vector<double> sig1,
                  std::vector<double> mu2,std::vector<double> sig2,
                  std::vector<double> &DIFF,std::vector<double> &ERR);

int SCCToggle(){

   gStyle->SetPalette(1); // nice colors  
   
   int M=0;  
   int nev=20;  // number of events to use for each plateau region 

   int rc=0;
   int method = gm2fieldUtil::Constants::kPhaseDerivative;

   std::vector<int> run;
   gm2fieldUtil::Import::GetRunList(run);
   const int NRUNS = run.size();

   int coilSet = -9; 
   std::cout << "Use bottom (0), top (1), or azi coils (-1)? ";
   std::cin  >> coilSet;

   // Surface coil data 
   std::vector<gm2field::surfaceCoils_t> sccData; 
   for(int i=0;i<NRUNS;i++) rc = gm2fieldUtil::RootHelper::GetSCCData(run[i],sccData);

   TGraph *gSCCb = GetSCCPlot(0,sccData); 
   TGraph *gSCCt = GetSCCPlot(1,sccData);

   double thr   = 10E-3;         // in A 
   double delta = 2.*60. + 40.;  // 2 min 40 sec   
   std::vector<double> sccOff,sccOn;
   rc = FindTransitionTimes(coilSet,thr,delta,sccData,sccOff,sccOn);

   bool sccStartOn = false;
   if(rc==1) sccStartOn = true;

   if(sccStartOn){
      std::cout << "SCC was ON at start of analysis" << std::endl;
   }else{
      std::cout << "SCC was OFF at start analysis" << std::endl;
   }

   double yMin_tr = 50E+3; 
   double yMax_tr = 60E+3;

   const int NToff = sccOff.size(); 
   const int NTon  = sccOn.size(); 
   TLine **tOff = new TLine*[NToff]; 
   TLine **tOn  = new TLine*[NTon]; 
   for(int i=0;i<NToff;i++){
      tOff[i] = new TLine(sccOff[i],yMin_tr,sccOff[i],yMax_tr); 
      tOff[i]->SetLineColor(kRed); 
   }
   for(int i=0;i<NTon;i++){
      tOn[i]  = new TLine(sccOn[i] ,yMin_tr,sccOn[i] ,yMax_tr); 
      tOn[i]->SetLineColor(kGreen+1); 
   }

   gm2fieldUtil::Graph::SetGraphParameters(gSCCb,20,kBlack);  
   gm2fieldUtil::Graph::SetGraphParameters(gSCCt,20,kRed  ); 

   TMultiGraph *mgSCC = new TMultiGraph();
   mgSCC->Add(gSCCb,"lp");  
   mgSCC->Add(gSCCt,"lp");  

   // Trolley data
   int trlyMethod = gm2fieldUtil::Constants::kPhaseDerivative;  
   std::vector<trolleyAnaEvent_t> trlyData;     
   for(int i=0;i<NRUNS;i++) rc = GetTrolleyData("",run[i],trlyMethod,trlyData);
   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   const int NTRLY = 17; 
   TGraph **gTRLY = new TGraph*[NTRLY]; 
   for(int i=0;i<NTRLY;i++){
      gTRLY[i] = GetTRLYTGraph(i,"GpsTimeStamp","freq",trlyData);
      gm2fieldUtil::Graph::SetGraphParameters(gTRLY[i],21,kBlack);
      gTRLY[i]->SetMarkerSize(0.5); 
   }

   // analysis
 
   double mean=0,stdev=0;
   std::vector<double> trial,diff,err; 
   std::vector<double> offMean_tr,offStdev_tr,onMean_tr,onStdev_tr; 
   std::vector<double> diff_tr,err_tr; 

   sccTrlyEvent_t scc_trial;                             // single trial  
   std::vector<sccTrlyEvent_t> scc_probe;                // single probe, a collection of trials 
   std::vector< std::vector<sccTrlyEvent_t> > sccEvent;  // all probes 
  
   char outpath[500];
   sprintf(outpath,"./scc-toggle_all-data_run-%05d.txt",run[0]); 
 
   for(int i=0;i<NTRLY;i++){
      // get statistics 
      rc = GetStats(i,nev,sccOff,trlyData,offMean_tr,offStdev_tr);
      rc = GetStats(i,nev,sccOn ,trlyData,onMean_tr ,onStdev_tr );
      rc = CalculateDiff(onMean_tr,onStdev_tr,offMean_tr,offStdev_tr,diff_tr,err_tr);
      if(rc!=0) return 1;
      // fill vectors 
      M = offMean_tr.size();  
      for(int j=0;j<M;j++){
	 scc_trial.probeID        = i+1;
         scc_trial.freq_bare      = offMean_tr[j];  
         scc_trial.freq_bare_err  = offStdev_tr[j];  
         scc_trial.freq_scc       = onMean_tr[j];  
         scc_trial.freq_scc_err   = onStdev_tr[j];  
         scc_trial.freq_diff      = diff_tr[j];  
         scc_trial.freq_diff_err  = err_tr[j]; 
	 scc_probe.push_back(scc_trial); 
	 // print to file 
	 PrintToFile_scc(outpath,i+1,j+1,scc_trial.freq_scc,scc_trial.freq_scc_err,scc_trial.freq_bare,scc_trial.freq_bare_err);
      }
      sccEvent.push_back(scc_probe); 
      // clear vectors 
      scc_probe.clear();
      offMean_tr.clear(); 
      offStdev_tr.clear(); 
      onMean_tr.clear(); 
      onStdev_tr.clear(); 
      diff_tr.clear(); 
      err_tr.clear(); 
   }   

   std::vector<double> md,md_err,md_aba,md_aba_err; 
   std::vector<double> x,y,z,dy,dz,aba,aba_err;
  
   double mean_scc=0,stdev_scc=0; 
   double mean_bare=0,stdev_bare=0; 
   double mean_aba=0,stdev_aba=0;

   double zMin_h=10E+6; 
   double zMax_h=-10E+6; 
 
   // print results 
   sprintf(outpath,"./scc-toggle_avg_run-%05d.txt",run[0]); 
   for(int i=0;i<NTRLY;i++){
      M = sccEvent[i].size();
      for(int j=0;j<M;j++){
	 x.push_back(sccEvent[i][j].freq_diff); 
	 y.push_back(sccEvent[i][j].freq_scc); 
	 z.push_back(sccEvent[i][j].freq_bare); 
	 dy.push_back(sccEvent[i][j].freq_scc_err); 
	 dz.push_back(sccEvent[i][j].freq_bare_err); 
	 // std::cout << Form("trial %02d: scc on = %.3lf +/- %.3lf, scc off = %.3lf +/- %.3lf, diff = %.3lf +/- %.3lf",
         //                   j+1,sccEvent[i][j].freq_scc,sccEvent[i][j].freq_scc_err,sccEvent[i][j].freq_bare,sccEvent[i][j].freq_bare_err,
	 //       	           sccEvent[i][j].freq_diff,sccEvent[i][j].freq_diff_err) << std::endl;
      }

      mean_scc   = gm2fieldUtil::Math::GetMean<double>(y); 
      stdev_scc  = gm2fieldUtil::Math::GetStandardDeviation<double>(y); 
      mean_bare  = gm2fieldUtil::Math::GetMean<double>(z); 
      stdev_bare = gm2fieldUtil::Math::GetStandardDeviation<double>(z);
      mean       = gm2fieldUtil::Math::GetMean<double>(x); 
      stdev      = gm2fieldUtil::Math::GetStandardDeviation<double>(x); 
      if(zMin_h>mean) zMin_h = mean;  
      if(zMax_h<mean) zMax_h = mean;  
      // get the ABA result
      if(sccStartOn){
	 GetDifference_ABA_sccFirst(y,dy,z,dz,aba,aba_err);
      }else{
	 GetDifference_ABA(y,dy,z,dz,aba,aba_err);
      }
      mean_aba  = gm2fieldUtil::Math::GetMean<double>(aba); 
      stdev_aba = gm2fieldUtil::Math::GetStandardDeviation<double>(aba); 
      std::cout << Form("PROBE %02d, mean = %.3lf +/- %.3lf Hz (%.3lf ppb), mean_aba = %.3lf +/- %.3lf Hz (%.3lf ppb)",
                        i+1,mean,stdev,stdev/0.06179,mean_aba,stdev_aba,stdev_aba/0.06179) << std::endl; 
      std::cout << "--------------------------------------------" << std::endl;
      // print to file
      PrintToFile_avg(outpath,i+1,mean,stdev,mean_aba,stdev_aba);
      // store results 
      md.push_back(mean); 
      md_err.push_back(stdev);
      md_aba.push_back(mean_aba); 
      md_aba_err.push_back(stdev_aba);
      // get ready for next probe
      x.clear();
      y.clear();
      z.clear();
      dy.clear();
      dz.clear(); 
      aba.clear();
      aba_err.clear();
   }
 
   // add 50% to the limits 
   zMin_h *= 1.5; 
   zMax_h *= 1.5; 

   // now make some plots  

   trolleyProbePosition_t trlyPos;
   rc = GetTrolleyProbePositions(trlyPos);  

   TGraph *gTrolleyPos = GetTRLYPositionsGraph();   

   TGraph2D *g2D = new TGraph2D();
   for(int i=0;i<17;i++) g2D->SetPoint(i,trlyPos.r[i],trlyPos.y[i],md_aba[i]); 

   const int xBin = 100; 
   double xMin_h  = -45; 
   double xMax_h  =  45; 
   const int yBin = 100; 
   double yMin_h  = -45; 
   double yMax_h  =  45;
   TH2D *h2D = new TH2D("h","h",xBin,xMin_h,xMax_h,yBin,yMin_h,yMax_h);

   TString Title;
   TCanvas *c1 = new TCanvas("c1","TRLY Data",1000,800);
   c1->cd();

   TPad *tpPad = new TPad("tpPad","",0,0,1,1); 
   tpPad->SetFillStyle(4000);  // transparent 
   tpPad->SetFrameFillStyle(0); 

   tpPad->Draw();
   tpPad->cd();
   gStyle->SetOptStat(0);
   g2D->SetHistogram(h2D);
   g2D->Draw("colz"); 
   g2D->GetHistogram()->SetMinimum(zMin_h); 
   g2D->GetHistogram()->SetMaximum(zMax_h); 
   g2D->Draw("colz");

   tpPad->cd(); 
   gTrolleyPos->Draw("p same"); 
 
   c1->Update();

   TCanvas *c2 = new TCanvas("c2","Trolley Probes",1200,800); 
   c2->Divide(4,5); 

   for(int i=0;i<NUM_TRLY;i++){
      c2->cd(i+1);
      gTRLY[i]->Draw("alp");
      gm2fieldUtil::Graph::SetGraphLabels(gTRLY[i],Form("Probe %02d",i+1),"","Frequency (Hz)"); 
      gm2fieldUtil::Graph::UseTimeDisplay(gTRLY[i]); 
      gTRLY[i]->Draw("alp");
      for(int j=0;j<NToff;j++) tOff[j]->Draw("same"); 
      for(int j=0;j<NTon;j++)  tOn[j]->Draw("same"); 
      c2->Update(); 
   }

   TCanvas *c3 = new TCanvas("c3","SCC Data",1200,600);
   c3->cd();

   mgSCC->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mgSCC,"SCC Data","","Current (A)"); 
   gm2fieldUtil::Graph::UseTimeDisplay(mgSCC); 
   mgSCC->Draw("a");
   c3->Update(); 

   return 0;
}
//______________________________________________________________________________
int PrintToFile_avg(const char *outpath,int probe,
                    double diff,double diff_err,
                    double diff_aba,double diff_aba_err){

   char outStr[500];
   sprintf(outStr,"%02d,%.3lf,%.3lf,%.3lf,%.3lf",probe,diff,diff_err,diff_aba,diff_aba_err); 

   std::ofstream outfile; 
   outfile.open(outpath,std::ios::app);
   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      outfile << outStr << std::endl;
      outfile.close();
   }
   return 0;
}
//______________________________________________________________________________
int PrintToFile_avg(const char *outpath,int probe,
                    double scc,double scc_err,
                    double bare,double bare_err,
                    double diff,double diff_err){

   char outStr[500];
   sprintf(outStr,"%02d,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf",probe,scc,scc_err,bare,bare_err,diff,diff_err); 

   std::ofstream outfile; 
   outfile.open(outpath,std::ios::app);
   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      outfile << outStr << std::endl;
      outfile.close();
   }
   return 0;
}
//______________________________________________________________________________
int PrintToFile_scc(const char *outpath,int probe,int trial,double scc_on,double scc_on_err,double bare,double bare_err){

   char outStr[500];
   sprintf(outStr,"%02d,%.3lf,%.3lf,%.3lf,%.3lf",probe,scc_on,scc_on_err,bare,bare_err); 

   std::ofstream outfile; 
   outfile.open(outpath,std::ios::app);
   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      outfile << outStr << std::endl;
      outfile.close();
   }
   return 0;
}
//______________________________________________________________________________
int CalculateDiff(std::vector<double> mu1,std::vector<double> sig1,
                  std::vector<double> mu2,std::vector<double> sig2,
                  std::vector<double> &DIFF,std::vector<double> &ERR){
   // compute shift in field
   // mu1,sig1 = mean and stdev for SCC on 
   // mu2,sig2 = mean and stdev for SCC off  
   const int N1 = mu1.size();
   const int N2 = mu2.size();
   if(N1!=N2){
      std::cout << "[CalculateDiff]: Error!  Vector sizes don't match! " << std::endl;
      std::cout << "                 mu1 size = " << N1 << std::endl;
      std::cout << "                 mu2 size = " << N2 << std::endl;
      return 1;
   }
   double diff=0,err=0;
   for(int i=0;i<N1;i++){
      diff = mu1[i] - mu2[i];
      std::cout << Form("%.3lf - %.3lf = %.3lf",mu1[i],mu2[i],diff) << std::endl; 
      err  = TMath::Sqrt( sig1[i]*sig1[i] + sig2[i]*sig2[i] );
      DIFF.push_back(diff); 
      ERR.push_back(err); 
   }
   return 0;
} 
//______________________________________________________________________________
int GetStats(int probe,int nev,std::vector<double> time,std::vector<trolleyAnaEvent_t> Data,
             std::vector<double> &MEAN,std::vector<double> &STDEV){

   const int N = time.size();
   int M=0,rc=0;
   double mean=0,stdev=0;
   std::vector<double> freq;
   std::vector<trolleyAnaEvent_t> Event; 
   for(int i=0;i<N;i++){
      // std::cout << "Looking for time " << gm2fieldUtil::GetStringTimeStampFromUTC(time[i]) << std::endl;
      // find events 
      rc = FilterSingle(probe,nev,time[i],Data,freq); 
      // now get mean of events 
      mean  = gm2fieldUtil::Math::GetMean<double>(freq); 
      stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(freq);  
      // store result
      MEAN.push_back(mean); 
      STDEV.push_back(stdev); 
      // set up for next time 
      freq.clear(); 
   }
  
   // now restrict to less than 20 results if necessary 
   // int cntr = MEAN.size();
   // while(cntr>20){
   //    MEAN.pop_back();
   //    STDEV.pop_back();
   //    cntr = MEAN.size();
   // }

   return 0;
}
//______________________________________________________________________________
TGraph *GetSCCPlot(int type,std::vector<gm2field::surfaceCoils_t> data){
   // construct the sum of the top (type=1), bottom (type=0) or azi (type=-1) coil currents 
   // to see what the SCC config is
   double sum=0;
   std::vector<double> x,y;  
   int M=4; // for azi coils
   if(type==0||type==1) M = 100;
   if( TMath::Abs(type)>1 ) return NULL;  
   int NEV = data.size();
   for(int i=0;i<NEV;i++){
      for(int j=0;j<M;j++){
	 if(type==0)  sum += data[i].BotCurrents[j];
	 if(type==1)  sum += data[i].TopCurrents[j];
	 if(type==-1) sum += data[i].AzCurrents[j];
      }
      if(type==0)  x.push_back( data[i].BotTime[0]/1E+9 ); 
      if(type==1)  x.push_back( data[i].TopTime[0]/1E+9 ); 
      if(type==-1) x.push_back( data[i].TopTime[0]/1E+9 ); // don't have an associated Az time...  
      y.push_back(sum); 
      // set up for next event 
      sum = 0;
   }
   TGraph *g = gm2fieldUtil::Graph::GetTGraph(x,y); 
   return g;
}

