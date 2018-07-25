// Test out the FXPR plotting functions   

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

int CombineResults(std::vector<double> lo,std::vector<double> hi,std::vector<double> &u);
int CorrectForDrift(double f0,std::vector<fixedProbeEvent_t> in,std::vector<fixedProbeEvent_t> &out); 
int GetStats(int nev,std::vector<double> Time,std::vector<fixedProbeEvent_t> fxprDataAvg,
             std::vector<double> &MEAN,std::vector<double> &STDEV); 
int GetStats(int probe,int nev,std::vector<double> time,std::vector<trolleyAnaEvent_t> Data,
             std::vector<double> &MEAN,std::vector<double> &STDEV);

int CalculateDiff(std::vector<double> mu1,std::vector<double> sig1,
                  std::vector<double> mu2,std::vector<double> sig2,
                  std::vector<double> &DIFF,std::vector<double> &ERR);

int SCCTest(){
   
   int M=0;  
   int nev=20;  // number of events to use for each plateau region 

   int rc=0;
   int method = gm2fieldUtil::Constants::kPhaseDerivative;

   std::vector<int> run;
   gm2fieldUtil::Import::GetRunList(run);
   const int NRUNS = run.size();

   // fixed probe data
   std::string fxprPath = "./input/probe-lists/fxpr-list.csv"; 
   std::vector<int> fxprList; 
   gm2fieldUtil::Import::ImportData1<int>(fxprPath,"csv",fxprList);

   // get fixed probe data 
   std::vector<gm2field::fixedProbeFrequency_t> fxprData;
   for(int i=0;i<NRUNS;i++) rc = gm2fieldUtil::RootHelper::GetFPFrequencies(run[i],fxprData);
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

   // double f0 = fxprDataAvg[0].freq; 
   // std::vector<fixedProbeEvent_t> fxprDataAvg_cor; 
   // rc = CorrectForDrift(f0,fxprDataAvg,fxprDataAvg_cor);

   double thr = 30.; // in Hz  

   // now for the trolley
   // std::vector<double> tLo_tr,tHi_tr; 
   // rc = FindTransitionTimes(0,thr,trlyData,tLo_tr,tHi_tr);  

   std::vector< std::vector<double> > tLo_tr,tHi_tr; // first index = probe, second = trial number 
   rc = FindTransitionTimes(thr,trlyData,tLo_tr,tHi_tr); 

   double yMin_tr = 50E+3; 
   double yMax_tr = 60E+3;

   // store all the transitions as TLines 
   std::vector<TLine *> lLine,hLine; 
   std::vector< std::vector<TLine *> > LO_tr,HI_tr;
   for(int i=0;i<NTRLY;i++){
      M = tLo_tr[i].size();
      for(int j=0;j<M;j++){
         // low line 
	 TLine *myLoLine = new TLine(tLo_tr[i][j],yMin_tr,tLo_tr[i][j],yMax_tr); 
         myLoLine->SetLineColor(kGreen+1); 
         myLoLine->SetLineWidth(2); 
         myLoLine->SetLineStyle(2); 
         lLine.push_back(myLoLine);
         // high line  
	 TLine *myHiLine = new TLine(tHi_tr[i][j],yMin_tr,tHi_tr[i][j],yMax_tr); 
         myHiLine->SetLineColor(kRed); 
         myHiLine->SetLineWidth(2); 
         myHiLine->SetLineStyle(2); 
         hLine.push_back(myHiLine); 
      }
      LO_tr.push_back(lLine); 
      HI_tr.push_back(hLine); 
   } 
   
   double yMin = 55E+3; 
   double yMax = 65E+3;

   double mean=0,stdev=0;

   // // transition times for the FXPR 
   // std::vector<double> tLo,tHi; 
   // rc = FindTransitionTimes(thr,fxprDataAvg,tLo,tHi);

   // const int NL = tLo.size();
   // TLine **LO = new TLine*[NL];
   // for(int i=0;i<NL;i++){
   //    LO[i] = new TLine(tLo[i],yMin,tLo[i],yMax);
   //    LO[i]->SetLineColor(kGreen+1); 
   //    LO[i]->SetLineWidth(2); 
   //    LO[i]->SetLineStyle(2); 
   // }   

   // const int NH = tHi.size();
   // TLine **HI= new TLine*[NH];
   // for(int i=0;i<NH;i++){
   //    HI[i] = new TLine(tHi[i],yMin,tHi[i],yMax);
   //    HI[i]->SetLineColor(kRed); 
   //    HI[i]->SetLineWidth(2); 
   //    HI[i]->SetLineStyle(2); 
   // }

   // std::cout << "FIXED PROBES" << std::endl;
   // double mean=0,stdev=0;
   // std::vector<double> loMean,loStdev,hiMean,hiStdev; 
   // rc = GetStats(nev,tLo,fxprDataAvg,loMean,loStdev);
   // rc = GetStats(nev,tHi,fxprDataAvg,hiMean,hiStdev);

   // std::vector<double> trial,diff,err; 
   // rc = CalculateDiff(hiMean,hiStdev,loMean,loStdev,diff,err);
  
   // const int NT = diff.size(); 
   // for(int i=0;i<NT;i++){ 
   //    trial.push_back(i+1); 
   //    std::cout << Form("trial %d: scc on = %.3lf +/- %.3lf, scc off = %.3lf +/- %.3lf, diff = %.3lf +/- %.3lf",
   //                      i+1,hiMean[i],hiStdev[i],loMean[i],loStdev[i],diff[i],err[i]) << std::endl;
   // }

   // double mean_diff  = gm2fieldUtil::Math::GetMean<double>(diff); 
   // double stdev_diff = gm2fieldUtil::Math::GetStandardDeviation<double>(diff); 
   // std::cout << Form("Mean Difference: %.3lf +/- %.3lf",mean_diff,stdev_diff) << std::endl;

   std::cout << "TROLLEY" << std::endl;
   std::vector<double> loMean_tr,loStdev_tr,hiMean_tr,hiStdev_tr; 
   std::vector<double> diff_tr,err_tr; 
   std::vector<double> tL,tH; 

   sccTrlyEvent_t scc_trial;                             // single trial  
   std::vector<sccTrlyEvent_t> scc_probe;                // single probe, a collection of trials 
   std::vector< std::vector<sccTrlyEvent_t> > sccEvent;  // all probes 
  
   for(int i=0;i<NTRLY;i++){
      // find times    
      rc = FindTransitionTimes(i,thr,trlyData,tL,tH); 
      // get statistics  
      rc = GetStats(i,nev,tLo_tr[i],trlyData,loMean_tr,loStdev_tr);
      rc = GetStats(i,nev,tHi_tr[i],trlyData,hiMean_tr,hiStdev_tr);
      rc = CalculateDiff(hiMean_tr,hiStdev_tr,loMean_tr,loStdev_tr,diff_tr,err_tr);
      // fill vectors 
      M = loMean_tr.size(); 
      for(int j=0;j<M;j++){
	 scc_trial.probeID        = i;
         scc_trial.freq_bare      = loMean_tr[j];  
         scc_trial.freq_bare_err  = loStdev_tr[j];  
         scc_trial.freq_scc       = hiMean_tr[j];  
         scc_trial.freq_scc_err   = hiStdev_tr[j];  
         scc_trial.freq_diff      = diff_tr[j];  
         scc_trial.freq_diff_err  = err_tr[j]; 
	 scc_probe.push_back(scc_trial); 
      }
      sccEvent.push_back(scc_probe); 
      // clear vectors 
      scc_probe.clear();
      // tL.clear(); 
      // tH.clear(); 
      loMean_tr.clear(); 
      loStdev_tr.clear(); 
      hiMean_tr.clear(); 
      hiStdev_tr.clear(); 
      diff_tr.clear(); 
      err_tr.clear(); 
   }   

   std::vector<double> x; 
 
   // print results 
   for(int i=0;i<NTRLY;i++){
      std::cout << Form("PROBE %02d",i) << std::endl;
      M = sccEvent[i].size();
      for(int j=0;j<M;j++){
	 x.push_back(sccEvent[i][j].freq_diff); 
	 std::cout << Form("trial %02d: scc on = %.3lf +/- %.3lf, scc off = %.3lf +/- %.3lf, diff = %.3lf +/- %.3lf",
                           j+1,sccEvent[i][j].freq_scc,sccEvent[i][j].freq_scc_err,sccEvent[i][j].freq_bare,sccEvent[i][j].freq_bare_err,
			   sccEvent[i][j].freq_diff,sccEvent[i][j].freq_diff_err) << std::endl;
      }
      mean  = gm2fieldUtil::Math::GetMean<double>(x); 
      stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(x); 
      std::cout << Form("MEAN = %.3lf +/- %.3lf Hz (%.3lf ppb)",mean,stdev,stdev/0.06179) << std::endl; 
      std::cout << "--------------------------------------------" << std::endl;
      x.clear();
   }

   // std::vector<double> trial_tr,diff_tr,err_tr; 
  
   // const int NT_tr = diff_tr.size(); 
   // for(int i=0;i<NT_tr;i++){ 
   //    trial_tr.push_back(i+1); 
   //    std::cout << Form("trial %d: scc on = %.3lf +/- %.3lf, scc off = %.3lf +/- %.3lf, diff = %.3lf +/- %.3lf",
   //                      i+1,hiMean_tr[i],hiStdev_tr[i],loMean_tr[i],loStdev_tr[i],diff_tr[i],err_tr[i]) << std::endl;
   // }

   // double mean_diff_tr  = gm2fieldUtil::Math::GetMean<double>(diff_tr); 
   // double stdev_diff_tr = gm2fieldUtil::Math::GetStandardDeviation<double>(diff_tr); 
   // std::cout << Form("Mean Difference: %.3lf +/- %.3lf",mean_diff_tr,stdev_diff_tr) << std::endl;

   // TGraphErrors *gDiff = gm2fieldUtil::Graph::GetTGraphErrors(trial,diff,err);
   // gm2fieldUtil::Graph::SetGraphParameters(gDiff,20,kBlack); 

   // TGraphErrors *gDiff_tr = gm2fieldUtil::Graph::GetTGraphErrors(trial_tr,diff_tr,err_tr);
   // gm2fieldUtil::Graph::SetGraphParameters(gDiff_tr,20,kBlack); 

   std::vector<fixedProbeEvent_t> lo1,hi1; 
   // std::vector<fixedProbeEvent_t> lo2,hi2; 
   // std::vector<fixedProbeEvent_t> lo3,hi3; 
   // rc = FilterSingle(nev,tLo[13],fxprDataAvg,lo1);
   // rc = FilterSingle(nev,tHi[13],fxprDataAvg,hi1);
   // rc = FilterSingle(nev,tLo[1],fxprDataAvg,lo2);
   // rc = FilterSingle(nev,tHi[1],fxprDataAvg,hi2);
   // rc = FilterSingle(nev,tLo[2],fxprDataAvg,lo3);
   // rc = FilterSingle(nev,tHi[2],fxprDataAvg,hi3);
  
   // Fixed probe average plot 
   TGraph *gFPAVG = GetTGraphNew(fxprDataAvg);
   gm2fieldUtil::Graph::SetGraphParameters(gFPAVG,21,kBlack); 

   // TGraph *gFPAVG_cor = GetTGraphNew(fxprDataAvg_cor);
   // gm2fieldUtil::Graph::SetGraphParameters(gFPAVG_cor,20,kBlue); 

   // TGraph *gLo1 = GetTGraphNew(lo1);
   // gm2fieldUtil::Graph::SetGraphParameters(gLo1,20,kGreen+2); 

   // TGraph *gHi1 = GetTGraphNew(hi1);
   // gm2fieldUtil::Graph::SetGraphParameters(gHi1,20,kRed+2); 

   // TGraph *gLo2 = GetTGraphNew(lo2);
   // gm2fieldUtil::Graph::SetGraphParameters(gLo2,20,kGreen+2); 

   // TGraph *gHi2 = GetTGraphNew(hi2);
   // gm2fieldUtil::Graph::SetGraphParameters(gHi2,20,kRed+2);
 
   // TGraph *gLo3 = GetTGraphNew(lo3);
   // gm2fieldUtil::Graph::SetGraphParameters(gLo3,20,kGreen+2); 

   // TGraph *gHi3 = GetTGraphNew(hi3);
   // gm2fieldUtil::Graph::SetGraphParameters(gHi3,20,kRed+2); 

   TMultiGraph *mg = new TMultiGraph();
   mg->Add(gFPAVG,"lp"); 
   // mg->Add(gFPAVG_cor,"lp"); 
   // mg->Add(gLo1  ,"lp"); 
   // mg->Add(gHi1  ,"lp"); 
   // mg->Add(gLo2  ,"lp"); 
   // mg->Add(gHi2  ,"lp"); 
   // mg->Add(gLo3  ,"lp"); 
   // mg->Add(gHi3  ,"lp"); 

   TCanvas *c1 = new TCanvas("c1","FXPR Data",1200,600);
   c1->Divide(1,2); 

   c1->cd(1);
   mg->Draw("a"); 
   gm2fieldUtil::Graph::SetGraphLabels(mg,"FXPR Avg","","Frequency (Hz)"); 
   gm2fieldUtil::Graph::SetGraphLabelSizes(mg,0.05,0.06); 
   gm2fieldUtil::Graph::UseTimeDisplay(mg); 
   mg->GetYaxis()->SetRangeUser(yMin,yMax); 
   mg->Draw("a");
   // for(int i=0;i<NL;i++) LO[i]->Draw("same"); 
   // for(int i=0;i<NH;i++) HI[i]->Draw("same"); 
   c1->Update();

   // c1->cd(2);
   // gDiff->Draw("alp");
   // gm2fieldUtil::Graph::SetGraphLabels(gDiff,"SCC on - SCC off","Trial","Frequency Difference (Hz)");
   // gm2fieldUtil::Graph::SetGraphLabelSizes(gDiff,0.05,0.06); 
   // gDiff->Draw("alp"); 
   // c1->Update();  

   const int NL_tr = tLo_tr[0].size(); 
   const int NH_tr = tHi_tr[0].size(); 

   TString Title;
   TCanvas *c2 = new TCanvas("c2","TRLY Data (0--3)",1200,600);
   c2->Divide(2,2); 

   TCanvas *c3 = new TCanvas("c3","TRLY Data (4--7)",1200,600);
   c3->Divide(2,2); 

   TCanvas *c4 = new TCanvas("c4","TRLY Data (8--11)",1200,600);
   c4->Divide(2,2);

   TCanvas *c5 = new TCanvas("c5","TRLY Data (12--15)",1200,600);
   c5->Divide(2,2);

   TCanvas *c6 = new TCanvas("c6","TRLY Data (16)",1200,600);

   for(int i=0;i<NTRLY;i++){
      Title = Form("Probe %d",i);
      if(i<=3)          c2->cd(i+1); 
      if(i>3   && i<8)  c3->cd(i-3); 
      if(i>=8  && i<12) c4->cd(i-7); 
      if(i>=12 && i<16) c5->cd(i-11); 
      if(i==16)         c6->cd(); 
      gTRLY[i]->Draw("alp");
      gm2fieldUtil::Graph::SetGraphLabels(gTRLY[i],Title,"","Frequency (Hz)"); 
      gm2fieldUtil::Graph::UseTimeDisplay(gTRLY[i]); 
      gTRLY[i]->GetYaxis()->SetRangeUser(yMin_tr,yMax_tr); 
      gTRLY[i]->Draw("alp");
      for(int j=0;j<NL_tr;j++) LO_tr[i][j]->Draw("same"); 
      for(int j=0;j<NH_tr;j++) HI_tr[i][j]->Draw("same"); 
      if(i<=3)          c2->Update(); 
      if(i>3   && i<8)  c3->Update(); 
      if(i>=8  && i<12) c4->Update(); 
      if(i>=12 && i<16) c5->Update(); 
      if(i==16)         c6->Update(); 
   }

   // TCanvas *c3 = new TCanvas("c3","TRLY Diff Data",1200,600);
   // c3->cd(); 

   // Title = Form("SCC on - SCC off (probe %d)",0);
   // gDiff_tr->Draw("alp"); 
   // gm2fieldUtil::Graph::SetGraphLabels(gDiff_tr,Title,"","Frequency Difference (Hz)"); 
   // gDiff_tr->Draw("alp"); 
   // c3->Update();

   return 0;
}
//______________________________________________________________________________
int CorrectForDrift(double f0,std::vector<fixedProbeEvent_t> in,std::vector<fixedProbeEvent_t> &out){
   fixedProbeEvent_t data; 
   double drift=0;
   const int N = in.size();
   for(int i=0;i<N;i++){
      drift        = in[i].freq - f0; 
      data.time    = in[i].time; 
      data.freq    = in[i].freq - drift; 
      data.freqErr = in[i].freqErr;
      out.push_back(data);  
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
      return 1;
   }
   double diff=0,err=0;
   for(int i=0;i<N1;i++){
      diff = mu1[i] - mu2[i];
      err  = TMath::Sqrt( sig1[i]*sig1[i] + sig2[i]*sig2[i] );
      DIFF.push_back(diff); 
      ERR.push_back(err); 
   }
   return 0;
} 
//______________________________________________________________________________
int GetStats(int nev,std::vector<double> time,std::vector<fixedProbeEvent_t> fxprDataAvg,
             std::vector<double> &MEAN,std::vector<double> &STDEV){

   const int N = time.size();
   int M=0,rc=0;
   double mean=0,stdev=0;
   std::vector<double> freq;
   std::vector<fixedProbeEvent_t> Event,Event_cor; 
   for(int i=0;i<N;i++){
      // find events 
      rc = FilterSingle(nev,time[i],fxprDataAvg,Event); 
      // now get mean of events 
      M = Event.size();
      for(int j=0;j<M;j++) freq.push_back( Event[j].freq ); 
      mean  = gm2fieldUtil::Math::GetMean<double>(freq); 
      stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(freq);  
      // store result
      MEAN.push_back(mean); 
      STDEV.push_back(stdev); 
      // set up for next time 
      Event.clear(); 
      freq.clear(); 
   }
  
   // now restrict to less than 20 results if necessary 
   int cntr = MEAN.size();
   while(cntr>20){
      MEAN.pop_back();
      STDEV.pop_back();
      cntr = MEAN.size();
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
   int cntr = MEAN.size();
   while(cntr>20){
      MEAN.pop_back();
      STDEV.pop_back();
      cntr = MEAN.size();
   }

   return 0;
}

