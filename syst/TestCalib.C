// test out our calibration ABA method

#include <cstdlib>

#include "TCanvas.h"
#include "TRandom3.h"

#include "gm2fieldMath.h"
#include "gm2fieldGraph.h"

#include "./src/FXPRFuncs.C"
#include "./src/Consolidate.C"
#include "./src/CustomAlgorithms.C"

double gXsize = 0.05; 
double gYsize = 0.06; 
 
double getRandomNumber(double range);
double driftFunc(double *par,double x); 
double fieldOscillation(double *par,double x); 

int TestCalib(){

   bool enableOsc=false;
   bool enableDrift=false;
   bool smearTime=false;

   char ans='t'; 
   std::cout << "Smear the timestamps? [y or n]: ";
   std::cin  >> ans;

   double smearPCT=0;
   if(ans=='y'){
      smearTime = true; 
      std::cout << "By how much? (enter percent): ";
      std::cin  >> smearPCT; 
   } 

   smearPCT /= 100; 

   ans = 't';
   std::cout << "Enable field drift? [y or n]: "; 
   std::cin  >> ans; 
   if(ans=='y') enableDrift = true; 

   ans = 't';
   std::cout << "Enable field oscillation? [y or n]: "; 
   std::cin  >> ans; 
   if(ans=='y') enableOsc = true; 
 
   // drift function parameters
   const int NPAR = 2;
   double offset = 5.; // 5 Hz 
   double slope  = 60./(3600.); // this is 1 ppm/hr = (60 Hz/(60 sec*60 min) 2.; // 2 Hz/sec  
   double par[NPAR] = {offset,slope}; 
  
   // field oscillation parameters 
   const int oscNPAR = 3; 
   double ampl  = 10.*0.06179;  // 10 ppb
   double freq  = 1./(2.*60.);  // 2 mins 
   double phase = getRandomNumber( 2.*TMath::Pi() );  // random phase on the range of 2pi rad
   double oscPar[oscNPAR] = {ampl,freq,phase};  
 
   // set up the test data
   double baseline = 6.; 
   double signal = 60.;
   double df=0;
   double dt = 4.*60.; // 4 mins 5.;  
   double arg_tA=0,arg_a=0,arg_eA=0,arg_tAs=0,arg_tA_last=0;
   double arg_tB=0,arg_b=0,arg_eB=0,arg_tBs=0,arg_tB_last=0;
   std::vector<double> tA,tB,fA,fB,eA,eB; 

   const int NPTS = 10; 
   for(int i=0;i<NPTS;i++){
      arg_tA = arg_tB_last + dt; 
      arg_tB = arg_tA + dt; 
      if(smearTime){
	 arg_tAs = arg_tA*(1. + getRandomNumber(smearPCT) );
	 arg_tBs = arg_tB*(1. + getRandomNumber(smearPCT) );
         while( arg_tAs<arg_tB_last ) {
	    arg_tAs = arg_tA*(1. + getRandomNumber(smearPCT) );
         }
         while( arg_tBs<arg_tAs ) {
	    arg_tBs = arg_tB*(1. + getRandomNumber(smearPCT) );
         }
	 // now we have valid timing that's smeared 
         arg_tA = arg_tAs;  
         arg_tB = arg_tBs;  
      }
      if(enableDrift){
	 arg_a = baseline + signal + driftFunc(par,arg_tA); 
	 arg_b = baseline + driftFunc(par,arg_tB);
      }else{
	 arg_a = baseline + signal;
	 arg_b = baseline; 
      }
      if(enableOsc){
	 arg_a += fieldOscillation(oscPar,arg_tA); 
	 arg_b += fieldOscillation(oscPar,arg_tB); 
      }
      arg_eA  = 1. + getRandomNumber(0.10); 
      arg_eB  = 1. + getRandomNumber(0.10); 
      tA.push_back(arg_tA);
      tB.push_back(arg_tB);
      fA.push_back(arg_a); 
      fB.push_back(arg_b);
      eA.push_back(arg_eA); 
      eB.push_back(arg_eB); 
      std::cout << Form("tA = %.3lf, fA = %.3lf ± %.3lf, tB = %.3lf, fB = %.3lf ± %.3lf",tA[i],fA[i],eA[i],tB[i],fB[i],eB[i]) << std::endl; 
      arg_tA_last = arg_tA;
      arg_tB_last = arg_tB;
   }

   std::cout << "=================================================" << std::endl; 

   bool useTimeWeight = true; 
   std::vector<double> trial,diff,diffErr; 
   int rc = GetDifference_ABA_final(useTimeWeight,tA,fA,eA,tB,fB,eB,diff,diffErr); 

   double mean=0,err=0,stdev=0; 
   rc = GetWeightedAverageStats(diff,diffErr,mean,err,stdev); 

   std::cout << Form("MEAN = %.3lf +/- %.3lf",mean,stdev) << std::endl;

   double min =-1E+3,max=1E+3; 
   const int N = diff.size();
   for(int i=0;i<N;i++){
      trial.push_back(i+1); 
      std::cout << Form("trial %02d: %.3lf",i+1,diff[i]) << std::endl;
      if(diff[i]<min) min = diff[i]; 
      if(diff[i]>max) max = diff[i]; 
   } 
 
   min *= 1. - 0.10; 
   max *= 1. + 0.10; 
  
   std::vector<double> TRIAL;
   for(int i=0;i<NPTS;i++) TRIAL.push_back(i+1); 

   // now make a plot of drift and oscillation data 
   double tMIN=1E+6,tMAX=-1E+6; 
   for(int i=0;i<NPTS;i++){
      if(tMIN>tA[i]) tMIN = tA[i]; 
      if(tMAX<tA[i]) tMAX = tA[i]; 
      if(tMIN>tB[i]) tMIN = tB[i]; 
      if(tMAX<tB[i]) tMAX = tB[i]; 
   }

   double arg_t=0,arg_d=0,arg_o=0; 
   std::vector<double> tt,dd,oo; 
   const int NN = 500; 
   double tStep = (tMAX-tMIN)/( (double)NN ); 
   for(int i=0;i<NN;i++){
      arg_t = tMIN + ( (double)i )*tStep;
      arg_d = driftFunc(par,arg_t); 
      arg_o = fieldOscillation(oscPar,arg_t); 
      tt.push_back(arg_t); 
      dd.push_back(arg_d); 
      oo.push_back(arg_o); 
   }

   TGraphErrors *gA = gm2fieldUtil::Graph::GetTGraphErrors(tA,fA,eA); 
   TGraphErrors *gB = gm2fieldUtil::Graph::GetTGraphErrors(tB,fB,eB);

   TGraph *gOsc     = gm2fieldUtil::Graph::GetTGraph(tt,oo); 
   TGraph *gDrift   = gm2fieldUtil::Graph::GetTGraph(tt,dd);  

   gm2fieldUtil::Graph::SetGraphParameters(gA    ,21,kBlue);  
   gm2fieldUtil::Graph::SetGraphParameters(gB    ,20,kRed);  
   gm2fieldUtil::Graph::SetGraphParameters(gOsc  ,21,kMagenta);  
   gm2fieldUtil::Graph::SetGraphParameters(gDrift,20,kGreen+2);  

   TMultiGraph *mg = new TMultiGraph();
   mg->Add(gA,"lp"); 
   mg->Add(gB,"lp"); 

   TMultiGraph *mgs = new TMultiGraph();
   mgs->Add(gDrift,"lp"); 
   mgs->Add(gOsc,"lp"); 

   TLegend *L = new TLegend(0.6,0.6,0.8,0.8);
   L->AddEntry(gA,"A Measurement","p"); 
   L->AddEntry(gB,"B Measurement","p"); 

   TString driftLabel = Form("Field Drift (%.3lf Hz/sec)",par[1]);
   TString oscLabel   = Form("Field Oscillation (A = %.3lf ppb, T = %.3lf sec, #phi = %.3lf rad)",oscPar[0]/0.06179,1./oscPar[1],oscPar[2]);
 
   TLegend *LS = new TLegend(0.6,0.6,0.8,0.8);
   LS->AddEntry(gDrift,driftLabel,"p"); 
   LS->AddEntry(gOsc  ,oscLabel,"p"); 

 
   TGraphErrors *gDiff = gm2fieldUtil::Graph::GetTGraphErrors(trial,diff,diffErr);
   gm2fieldUtil::Graph::SetGraphParameters(gDiff,20,kBlack);  

   TString TitleM     = Form("Measurements (Known Test Value = %.3lf Hz)",signal);
   TString TitleD     = Form("ABA Difference");
   TString xAxisTitle = Form("Trial");
   TString yAxisTitle = Form("Frequency Difference (Hz)");

   if(smearTime) TitleM += Form(" [Timestamps smeared by up to %.0lf%]",smearPCT*100);

   TCanvas *c1 = new TCanvas("c1","ABA Difference Trial",1200,600);
   c1->Divide(1,2); 

   c1->cd(1);
   mg->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mg,TitleM,"Time","Frequency (Hz)");
   gm2fieldUtil::Graph::SetGraphLabelSizes(mg,gXsize,gYsize); 
   mg->GetYaxis()->SetRangeUser(min,max);
   mg->Draw("a");
   L->Draw("same"); 
   c1->Update();

   c1->cd(2);
   gDiff->Draw("alp");
   gm2fieldUtil::Graph::SetGraphLabels(gDiff,TitleD,xAxisTitle,yAxisTitle);
   gm2fieldUtil::Graph::SetGraphLabelSizes(gDiff,gXsize,gYsize); 
   gDiff->GetYaxis()->SetRangeUser(min,max);
   gDiff->Draw("alp");
   c1->Update(); 

   TCanvas *c2 = new TCanvas("c2","Environmental Effects",1200,600);

   c2->cd();
   mgs->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mgs,"Environment","Time","Frequency (Hz)");
   // gm2fieldUtil::Graph::SetGraphLabelSizes(mgs,gXsize,gYsize); 
   // mgs->GetYaxis()->SetRangeUser(min,max);
   mgs->Draw("a");
   LS->Draw("same"); 
   c2->Update();

   return 0;
}
//______________________________________________________________________________
double getRandomNumber(double range){
   TRandom3 *R = new TRandom3(0); 
   double r = R->Rndm();
   double val = (-1.)*range + 2.*r*range; 
   delete R;
   return val; 
}
//______________________________________________________________________________
double driftFunc(double *par,double x){
   double f = par[0] + par[1]*x; 
   return f;
}
//______________________________________________________________________________
double fieldOscillation(double *par,double x){
   double f = par[0]*TMath::Sin( 2.*TMath::Pi()*par[1]*x + par[2] ); 
   return f;
}
