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
#include "gm2fieldFunc.h"
#include "TemperatureSensor.h"

#include "./include/Constants.h"
#include "./include/fixedProbeEvent.h"
#include "./include/trolleyAnaEvent.h"
#include "./include/sccEvent.h" 

#include "./src/BlindFuncs.C"
#include "./src/InputManager.C"
#include "./src/FitFuncs.C"
#include "./src/CustomUtilities.C"
#include "./src/CustomImport.C"
#include "./src/CustomExport.C"
#include "./src/CustomMath.C"
#include "./src/CustomGraph.C"
#include "./src/CustomAlgorithms.C"

int PrintToFile_csv(const char *outpath,std::vector<double> off,std::vector<double> on); 

int PrintToFile_avg(const char *outpath,int probe,double scc_on,double scc_on_err,double bare,double bare_err);
int PrintToFile_avg(const char *outpath,int probe,double diff,double diff_err,double diff_aba,double diff_aba_err);
int PrintToFile_scc(const char *outpath,int probe,int trial,double scc_on,double scc_on_err,double bare,double bare_err); 

int CombineResults(std::vector<double> lo,std::vector<double> hi,std::vector<double> &u);

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

   int probe=-1,coilSet=-9,axis=-1; 
   std::cout << "Use bottom (0), top (1), or azi coils (-1)? ";
   std::cin  >> coilSet;
   std::cout << "Enter axis (0 = x, 1 = y, 2 = z): ";
   std::cin  >> axis;
   std::cout << "Enter probe: ";
   std::cin  >> probe; 

   char AXIS = 't'; 
   if(axis==0) AXIS = 'x';
   if(axis==1) AXIS = 'y';
   if(axis==2) AXIS = 'z';

   std::string prodVersion = "v9_21_01"; 

   // get trolley data
   std::vector<trolleyAnaEvent_t> trlyData;
   for(int i=0;i<NRUNS;i++) rc = GetTrolleyData(run[i],method,trlyData,prodVersion);
   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   TGraph *gTR = GetTRLYTGraph(probe-1,"GpsTimeStamp","freq",trlyData);
   gm2fieldUtil::Graph::SetGraphParameters(gTR,20,kBlack);

   // determine the correct ordering of the SCC on/off cycles 
   std::vector<surfaceCoilEvent_t> sccData;
   for(int i=0;i<NRUNS;i++) rc = GetSurfaceCoilData(run[i],sccData,prodVersion);
   if(rc!=0) return 1;

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

   double thr   = 10E-3;         // in A 
   double delta = 2.*60. + 40.;  // 2 min 40 sec   
   std::vector<double> sccOff,sccOn;
   // rc = FindTransitionTimes(coilSet,axis,thr,delta,sccData,sccOff,sccOn);

   int runPeriod = 1;
   char trl[200]; 
   sprintf(trl,"tr%c",AXIS);
   std::string trLabel = trl;
   rc = LoadSCCTimes(probe,runPeriod,prodVersion,trLabel,sccOff,sccOn);

   bool sccStartOn = false;
   // if(rc==1) sccStartOn = true;

   if(sccOn[0]<sccOff[0]) sccStartOn = true;

   if(sccStartOn){
      std::cout << "SCC was ON at start of analysis" << std::endl;
   }else{
      std::cout << "SCC was OFF at start analysis" << std::endl;
   }

   double yMin_scc = -30; 
   double yMax_scc =  30;

   double yMin_tr = 50E+3; 
   double yMax_tr = 60E+3;

   const int NToff = sccOff.size(); 
   const int NTon  = sccOn.size(); 
   TLine **tOff = new TLine*[NToff]; 
   TLine **tOn  = new TLine*[NTon]; 
   TLine **tOff_tr = new TLine*[NToff]; 
   TLine **tOn_tr  = new TLine*[NTon]; 
   for(int i=0;i<NToff;i++){
      tOff[i] = new TLine(sccOff[i],yMin_scc,sccOff[i],yMax_scc); 
      tOff[i]->SetLineColor(kRed);
      tOff[i]->SetLineWidth(2);  
      tOff[i]->SetLineStyle(2);  
      tOff_tr[i] = new TLine(sccOff[i],yMin_tr,sccOff[i],yMax_tr); 
      tOff_tr[i]->SetLineColor(kRed);
      tOff_tr[i]->SetLineWidth(2);  
      tOff_tr[i]->SetLineStyle(2);  
   }
   for(int i=0;i<NTon;i++){
      tOn[i]  = new TLine(sccOn[i] ,yMin_scc,sccOn[i] ,yMax_scc); 
      tOn[i]->SetLineColor(kGreen+1); 
      tOn[i]->SetLineWidth(2);  
      tOn[i]->SetLineStyle(2);  
      tOn_tr[i]  = new TLine(sccOn[i] ,yMin_tr,sccOn[i] ,yMax_tr); 
      tOn_tr[i]->SetLineColor(kGreen+1); 
      tOn_tr[i]->SetLineWidth(2);  
      tOn_tr[i]->SetLineStyle(2);  
   }

   TCanvas *c1 = new TCanvas("c1","SCC Data",1200,600);
   c1->Divide(1,2); 

   c1->cd(1);
   mgSCC->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mgSCC,"SCC Data","","Total Current (A)"); 
   gm2fieldUtil::Graph::UseTimeDisplay(mgSCC); 
   mgSCC->Draw("a");
   for(int i=0;i<NToff;i++) tOff[i]->Draw("same"); 
   for(int i=0;i<NTon ;i++) tOn[i]->Draw("same"); 
   c1->Update();

   c1->cd(2);
   gTR->Draw("alp");
   gm2fieldUtil::Graph::SetGraphLabels(gTR,"TRLY Data","","Frequency (Hz)"); 
   gm2fieldUtil::Graph::UseTimeDisplay(gTR); 
   gTR->Draw("alp");
   for(int i=0;i<NToff;i++) tOff_tr[i]->Draw("same"); 
   for(int i=0;i<NTon ;i++) tOn_tr[i]->Draw("same"); 
   c1->Update();

   rc = PrintToFile_csv("test.csv",sccOff,sccOn);  

   return 0;
}
//______________________________________________________________________________
int PrintToFile_csv(const char *outpath,std::vector<double> off,std::vector<double> on){

   const int NOFF = off.size(); 
   const int NON  = on.size();
   if(NOFF!=NON){
      std::cout << "Unequal number of transitions! " << std::endl;
      std::cout << "off = " << NOFF << " on = " << NON << std::endl;
      return 1;
   } 

   unsigned long theTime=0; 
   char outStr[200]; 

   std::ofstream outfile;
   outfile.open(outpath);

   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      for(int i=0;i<NOFF;i++){
         theTime = (unsigned long)off[i];
	 sprintf(outStr,"%s",gm2fieldUtil::GetStringTimeStampFromUTC(theTime).c_str());
         outfile << outStr << std::endl; 
         theTime = (unsigned long)on[i]; 
	 sprintf(outStr,"%s",gm2fieldUtil::GetStringTimeStampFromUTC(theTime).c_str());
         outfile << outStr << std::endl; 
      }
      outfile.close();
   } 
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
