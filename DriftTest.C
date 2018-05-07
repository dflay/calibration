// Compute DeltaB for TRLY data sets  

#include <cstdlib> 
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>  
#include <cmath> 

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
#include "gm2fieldImport.h"
#include "TemperatureSensor.h"

#include "./include/Constants.h"
#include "./include/drift.h"
#include "./include/fixedProbeEvent.h"
#include "./include/plungingProbeAnaEvent.h"
#include "./include/trolleyAnaEvent.h"

#include "./src/FXPRFuncs.C"
#include "./src/Consolidate.C"
#include "./src/CustomMath.C"
#include "./src/CustomGraph.C"
#include "./src/CustomImport.C"
#include "./src/CustomExport.C"
#include "./src/CustomAlgorithms.C"
#include "./src/DeltaBFuncs.C"

double LagrangePolynomialInterpolation(double x_pr,std::vector<double> x,std::vector<double> f); 
TGraph *GetLagrangeInterpolatedTGraph(int method,std::vector<int> probeList,
                                      unsigned long long tStart,unsigned long long tStop,unsigned long long tStep,
                                      std::vector<gm2field::fixedProbeFrequency_t> fxprData); 

TSpline3 *GetSpline(int method,std::vector<int> probe,std::vector<gm2field::fixedProbeFrequency_t> fxprData); 
// TGraph *GetInterpolatedTGraph(int method,std::vector<int> probe,
//                               unsigned long long tStart,unsigned long long tStop,unsigned long long tStep,
//                               std::vector<gm2field::fixedProbeFrequency_t> fxprData); 

int DriftTest(){

   int rc=0;
   int method = gm2fieldUtil::Constants::kPhaseDerivative;

   std::cout << "----------------------------------" << std::endl;
   std::cout << "Drift Test " << std::endl;

   std::vector<int> run;
   rc = gm2fieldUtil::Import::GetRunList(run);
   const int NRUNS = run.size(); 
   std::cout << "Processing runs ";
   for(int i=0;i<NRUNS;i++) std::cout << run[i] << " ";
   std::cout << std::endl;

   // Fixed Probe data
   std::string fxprPath = "./input/probe-lists/fxpr-list.csv";
   std::vector<int> fxprList;
   gm2fieldUtil::Import::ImportData1<int>(fxprPath,"csv",fxprList);

   // get the fxpr data       
   std::vector<gm2field::fixedProbeFrequency_t> fxpr1,fxpr2; 
   rc = gm2fieldUtil::RootHelper::GetFPFrequencies(run[0],fxpr1);
   rc = gm2fieldUtil::RootHelper::GetFPFrequencies(run[1],fxpr2);

   std::vector<fixedProbeEvent_t> fxprAvg1;
   rc = GetAverageFXPRVectorsNew(method,0,0,0,0,fxprList,fxpr1,fxprAvg1);

   std::vector<fixedProbeEvent_t> fxprAvg2;
   rc = GetAverageFXPRVectorsNew(method,0,0,0,0,fxprList,fxpr2,fxprAvg2);

   const int N1 = fxprAvg1.size(); 
   const int N2 = fxprAvg2.size(); 

   if (rc!=0) {
      std::cout << "No data.  Exiting..." << std::endl;
      return 1;
   }
 
   std::vector<double> runStart,runStop;
   runStart.push_back( fxprAvg1[0].time ); 
   runStart.push_back( fxprAvg2[0].time ); 
   runStop.push_back( fxprAvg1[N1-1].time ); 
   runStop.push_back( fxprAvg2[N2-1].time ); 

   std::cout << "FXPR Time Stamps" << std::endl;
   for(int i=0;i<NRUNS;i++){
      std::cout << Form("Run %d, start = %s, stop = %s",run[i],
                        gm2fieldUtil::GetStringTimeStampFromUTC((int)runStart[i]).c_str(),
                        gm2fieldUtil::GetStringTimeStampFromUTC((int)runStop[i]).c_str() ) << std::endl;
   }

   std::cout << "Getting R2R graph..." << std::endl;
   std::vector<double> stats; 
   TGraph *gR2R = GetDriftTGraphR2R(method,run,fxprList,stats);
   gm2fieldUtil::Graph::SetGraphParameters(gR2R,20,kRed); 
   std::cout << "--> Done." << std::endl;

   std::cout << Form("SCC effect    = %.3lf +/- %.3lf Hz",stats[0],stats[1]) << std::endl;
   std::cout << Form("drift applied = %.3lf +/- %.3lf Hz",stats[2],stats[3]) << std::endl;

   std::cout << "Getting graphs..." << std::endl;

   // Fixed probe average plot
   TGraph *gFPAVG1 = GetTGraphNew(fxprAvg1);
   TGraph *gFPAVG2 = GetTGraphNew(fxprAvg2);
   std::cout << "--> Done." << std::endl;

   gm2fieldUtil::Graph::SetGraphParameters(gFPAVG1,21,kBlack);
   gm2fieldUtil::Graph::SetGraphParameters(gFPAVG2,21,kBlack);

   // gFPAVG1->SetMarkerSize(1.0); 
   // gFPAVG2->SetMarkerSize(1.0); 
   // gR2R->SetMarkerSize(1.0); 

   // Fixed probe interpolated plots 
   // std::vector<double> stats2; 
   // TGraph *gFPINT = GetInterpolatedTGraph(method,fxprList,runStop[0],runStart[1],tStep,fxprData,stats2);
   // gm2fieldUtil::Graph::SetGraphParameters(gFPINT,20,kBlack);

   // TGraph *gFPINT2 = GetLagrangeInterpolatedTGraph(method,fxprList,tStart,tStop,1E+9,fxprData); 
   // gm2fieldUtil::Graph::SetGraphParameters(gFPINT2,21,kGreen+1);

   // TLine *xStart = new TLine(runStop[0]/1E+9 ,MIN,runStop[0]/1E+9 ,MAX);
   // xStart->SetLineColor(kRed); 
 
   // TLine *xStop  = new TLine(runStart[1]/1E+9,MIN,runStart[1]/1E+9,MAX); 
   // xStop->SetLineColor(kRed);  

   // TSpline3 *gSpline = GetSpline(method,fxprList,fxprData); 
   // gSpline->SetMarkerColor(kBlue);  
   // gSpline->SetLineWidth(2); 
   // gSpline->SetLineColor(kBlue);  

   TMultiGraph *mg = new TMultiGraph();
   // mg->Add(gFPINT2,"lp");
   mg->Add(gFPAVG1,"p");
   mg->Add(gFPAVG2,"p");
   mg->Add(gR2R   ,"p");  
   // mg->Add(gFPINT ,"p");

   TLegend *L = new TLegend(0.6,0.6,0.8,0.8); 
   L->AddEntry(gFPAVG1,"Average"      ,"p");  
   L->AddEntry(gR2R   ,"R2R corrected","p"); 
   // L->AddEntry(gFPINT ,"Linear Interpolated","p"); 
   // L->AddEntry(gFPINT2,"Lagrange Polynomial Interpolated","p"); 
   // L->AddEntry(gSpline,"Spline Fit","l");  

   TCanvas *c1 = new TCanvas("c1","Fixed Probe Data",1200,600); 
   // c1->Divide(1,3);

   c1->cd();
   mg->Draw("ap"); 
   gm2fieldUtil::Graph::SetGraphLabels(mg,"Fixed Probe Avg","","Frequency (Hz)");
   gm2fieldUtil::Graph::UseTimeDisplay(mg); 
   gm2fieldUtil::Graph::SetGraphLabelSizes(mg,0.04,0.04);
   // mg->GetYaxis()->SetRangeUser(MIN,MAX);  
   mg->Draw("ap");
   // xStart->Draw("same"); 
   // xStop->Draw("same"); 
   // L->Draw("same"); 
   c1->Update();  

   return 0;
}
//______________________________________________________________________________
TGraph *GetLagrangeInterpolatedTGraph(int method,std::vector<int> probeList,
                                      unsigned long long tStart,unsigned long long tStop,unsigned long long tStep,
                                      std::vector<gm2field::fixedProbeFrequency_t> fxprData){
   // use the Lagrange polynomial interpolation method to fit the data
   // grab vectors of the data 
   unsigned long long t0=tStart;
   std::vector<unsigned long long> time;
   std::vector<double> TIME,FREQ;  
   GetAverageFXPRVectors(method,t0,tStart,tStop,tStep,probeList,fxprData,time,FREQ); 
   const int NPTS = time.size(); 
   for(int i=0;i<NPTS;i++) TIME.push_back(time[i]/1E+9); 

   std::cout << "Using " << NPTS << " points" << std::endl;

   // now build the function
   // use a step size that's half the original for a test... 
   double t_pr=0,freq=0,t_arg=0;
   unsigned long long fineStep = tStep/3;
   int N2 = (tStop-tStart)/fineStep; 
   std::cout << "Interpolating for " << N2 << " points" << std::endl;
   std::vector<double> X,Y;
   for(int i=0;i<N2;i++){
      t_pr = (tStart + ( (double)i )*fineStep - t0)/1E+9; 
      freq = LagrangePolynomialInterpolation(t_pr,TIME,FREQ); // find freq at t = t'
      t_arg = t_pr + t0/1E+9;  
      if(freq<0) freq = 0;
      X.push_back(t_arg); 
      Y.push_back(freq);
      // std::cout << Form("%s: %.3lf",gm2fieldUtil::GetStringTimeStampFromUTC(t_arg).c_str(),freq) << std::endl; 
   } 

   TGraph *g = gm2fieldUtil::Graph::GetTGraph(X,Y);
   return g;
}
//______________________________________________________________________________
double LagrangePolynomialInterpolation(double x_pr,std::vector<double> x,std::vector<double> f){
   // construct the lagrange polynomial of order n to interpolate between the data 
   const int n = x.size();  
   std::vector<double> L; 
   double num=1,den=1,arg=0;
   for(int i=0;i<n;i++){
      for(int j=0;j<n;j++){
	 if(i!=j){
	    num *= (x_pr-x[j]);
	    den *= (x[i]-x[j]);
            // std::cout << Form("num *= (%.3lf-%.3lf), den *= (%.3lf-%.3lf)",x_pr,x[j],x[i],x[j]) << std::endl;
         }
         // std::cout << Form("i = %d, j = %d, num = %.3lf, den = %.3lf",i,j,num,den) << std::endl;
      }
      arg = num/den;
      L.push_back(arg);
      // reset 
      num=1;
      den=1;
      // std::cout << Form("L[%d](%.0lf) = %.31f",i,x_pr,arg) << std::endl;  
   }

   // compute the function P(x) = sum_i[ fi(x)*Li(x) ]  
   double P=0; 
   for(int i=0;i<n;i++) P += f[i]*L[i]; 
 
   return P; 
}
//______________________________________________________________________________
TSpline3 *GetSpline(int method,std::vector<int> probe,std::vector<gm2field::fixedProbeFrequency_t> fxprData){

   // determine average FXPR frequency across a list of probes,  
   // using a range defined by tStart < t < tStop with a step size tStep 
   // the times tStart, tStop are necessarily end/start points of runs that 
   // are far in time (more than a second or two)   

   // populate all the data we have in fxprData into vectors 
   std::vector<double> TIME,FREQ;
   unsigned long long aTime;
   double mean=0,stdev=0;
   const int NEvents = fxprData.size();
   for(int i=0;i<NEvents;i++){
      aTime = fxprData[i].GpsTimeStamp[probe[0]]; // use the zeroth probe time...  
      GetAverageFXPR(method,aTime,probe,fxprData,mean,stdev);
      TIME.push_back(aTime/1E+9); 
      FREQ.push_back(mean); 
   } 

   // now get a graph 
   TGraph *g = gm2fieldUtil::Graph::GetTGraph(TIME,FREQ);

   TSpline3 *gSpl = new TSpline3("gSpl",g);
   return gSpl;

}

