// Compare results from E989 to E821      

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
#include "./include/grad_meas.h"
#include "./include/deltab.h"
#include "./include/results.h"
#include "./include/perturbation.h" 

#include "./src/InputManager.C"
#include "./src/FitFuncs.C"
// #include "./src/FXPRFuncs.C"
#include "./src/TRLYFuncs.C"
#include "./src/CalibFuncs.C"
// #include "./src/Consolidate.C"
#include "./src/CustomUtilities.C"
#include "./src/CustomMath.C"
#include "./src/CustomImport.C"
#include "./src/CustomExport.C"
#include "./src/CustomAlgorithms.C"
#include "./src/CustomGraph.C"

enum myUnits{
   Hz  = 0,
   ppm = 1,
   ppb = 2
};

enum expType{
   kE821 = 0,
   kE989 = 1
};

double gMarkerSize = 0.8;
double gXSize      = 0.05;  
double gYSize      = 0.06;  

int Convert(int units,int exp,bool relToTR1,std::vector<result_prod_t> &data,bool flipSign=false); 
int LoadTableData(const char *inpath,std::vector<result_prod_t> &data);
int PrintToScreen(TString yAxis,std::vector<result_prod_t> data); 

TGraphErrors *GetGraph(TString yAxis,std::vector<result_prod_t> data); 

int Compare_prod(bool isBlind,std::string date,std::string blindLabel){

   int rc=0;

   char outDir[200],inpath_e989[200];
   sprintf(outDir,"./output"); 
   if(isBlind)  sprintf(outDir,"%s/blinded/%s",outDir,blindLabel.c_str());
   if(!isBlind) sprintf(outDir,"%s/unblinded",outDir);
   sprintf(inpath_e989,"%s/%s/results_all-probes.csv",outDir,date.c_str()); 

   // load E989 data 
   std::vector<result_prod_t> e989; 
   rc = LoadTableData(inpath_e989,e989);
   if(rc!=0) return 1;

   // load E821 data 
   char inpath_e821[200]; 
   sprintf(inpath_e821,"./input/E821/results.csv");
   std::vector<result_prod_t> e821; 
   rc = LoadTableData(inpath_e821,e821);

   // convert the E821 data to match our convention 
   bool relativeToTR1 = false;
   bool flipSign      = true;  
   rc = Convert(ppm,kE821,relativeToTR1,e821,flipSign); 

   // convert E989 to match E821 units 
   relativeToTR1 = false;
   flipSign      = false;  
   rc = Convert(ppm,kE989,relativeToTR1,e989,flipSign); 

   std::cout << "---------------- Calibration Coefficients ----------------  " << std::endl;
   std::cout << "E821" << std::endl;
   rc = PrintToScreen("raw",e821); 
   std::cout << "E989" << std::endl; 
   rc = PrintToScreen("raw",e989);  

   // build plots 
   TGraphErrors *gE821 = GetGraph("raw",e821); 
   TGraphErrors *gE989 = GetGraph("raw",e989); 

   gm2fieldUtil::Graph::SetGraphParameters(gE821,21,kBlack); 
   gm2fieldUtil::Graph::SetGraphParameters(gE989,20,kRed); 

   TMultiGraph *mg = new TMultiGraph();
   // mg->Add(gE821,"lp"); 
   mg->Add(gE989,"lp"); 

   TLegend *L = new TLegend(0.6,0.6,0.8,0.8);
   L->AddEntry(gE821,"E821","p"); 
   L->AddEntry(gE989,"E989","p"); 

   TCanvas *c1 = new TCanvas("c1","Comparing E821 to E989",1200,600); 

   c1->cd(); 
   mg->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mg,"Calibration Comparison (Relative to Probe 1)","Probe Number","Calibration (ppm)"); 
   mg->Draw("a");
   // L->Draw("same");
   c1->Update();

   return 0;
}

//______________________________________________________________________________
int Convert(int units,int exp,bool relToTR1,std::vector<result_prod_t> &data,bool flipSign){
   // convert to the E821 convention
   // units: ppm, subtract off center probe result
   double ref     = data[0].diff; 
   double refFree = data[0].diffFree; 
   double arg=0;
   double sign=1; 
   double sf = 1.;

   if(flipSign) sign = -1; 

   const int n = data.size();
   for(int i=0;i<n;i++){
      if(exp==kE989){
	 if(relToTR1){
	    data[i].diff     -= ref;
	    data[i].diffFree -= refFree;
	 }
	 if(units==ppm){
	    sf = 61.79; 
	 }else if(units==ppb){
	    sf = 0.06179; 
	 }
	 data[i].diff     /= sf; 
	 data[i].diffFree /= sf;
	 data[i].diffErr  /= sf;  
	 data[i].mErr     /= sf;  
	 data[i].pErr     /= sf;
      }  
      data[i].diff     *= sign; 
      data[i].diffFree *= sign;
   } 
   return 0;
}
//______________________________________________________________________________
int LoadTableData(const char *inpath,std::vector<result_prod_t> &data){

   result_prod_t aPoint; 

   std::string sp,sd,sdf,sse,sme,sfe;

   std::ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      std::cout << "Cannot open the file: " << inpath << std::endl;
      return 1;
   }else{
      while( !infile.eof() ){
	 std::getline(infile,sp ,',');
	 std::getline(infile,sd ,',');
	 std::getline(infile,sdf,',');
	 std::getline(infile,sse,',');
	 std::getline(infile,sme,',');
	 std::getline(infile,sfe);
	 aPoint.diff     = std::atof( sd.c_str()  ); 
	 aPoint.diffFree = std::atof( sdf.c_str() ); 
	 aPoint.diffErr  = std::atof( sse.c_str() ); 
	 aPoint.mErr     = std::atof( sme.c_str() ); 
	 aPoint.pErr     = std::atof( sfe.c_str() ); 
	 data.push_back(aPoint);
      }
      data.pop_back();
      infile.close();
   }
   return 0;
}
//______________________________________________________________________________
TGraphErrors *GetGraph(TString yAxis,std::vector<result_prod_t> data){
   double arg=0,argErr=0;
   std::vector<double> x,y,ey;
   const int N = data.size();
   for(int i=0;i<N;i++){
      x.push_back( (double)i+1 );  
      if(yAxis=="raw"){
	 arg    = data[i].diff;
	 argErr = TMath::Sqrt(data[i].diffErr*data[i].diffErr + data[i].mErr*data[i].mErr);
      }else if(yAxis=="free"){
	 arg    = data[i].diffFree;
	 argErr = TMath::Sqrt(data[i].diffErr*data[i].diffErr + data[i].mErr*data[i].mErr + data[i].pErr*data[i].pErr);
      }
      y.push_back(arg);
      ey.push_back(argErr);  
   }
   TGraphErrors *g = gm2fieldUtil::Graph::GetTGraphErrors(x,y,ey);
   return g;
}
//______________________________________________________________________________
int PrintToScreen(TString yAxis,std::vector<result_prod_t> data){
   double arg=0,argErr=0;
   const int N = data.size();
   for(int i=0;i<N;i++){
      if(yAxis=="raw"){
	 arg    = data[i].diff;
	 argErr = TMath::Sqrt(data[i].diffErr*data[i].diffErr + data[i].mErr*data[i].mErr);
      }else if(yAxis=="free"){
	 arg    = data[i].diffFree;
	 argErr = TMath::Sqrt(data[i].diffErr*data[i].diffErr + data[i].mErr*data[i].mErr + data[i].pErr*data[i].pErr);
      }
      std::cout << Form("%02d: %.3lf +/- %.3lf",i+1,arg,argErr) << std::endl;
   }
   return 0;
}
