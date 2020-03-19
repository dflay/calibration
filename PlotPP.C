// Plot the Plunging Probe data   

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

int AddPPMultiGraph(std::string xAxis,std::string yAxis,std::vector<plungingProbeAnaEvent_t> data,
                    TMultiGraph *mg,TLegend *L); 

TGraph *getTGraph(int run,std::string xAxis,std::string yAxis,std::vector<plungingProbeAnaEvent_t> data); 

int FillPPVector3(int ppDAQ_run,std::string axis,std::vector<plungingProbeAnaEvent_t> data,
                  std::vector<double> &x); 

int PlotPP(){

   int rc=0;

   std::vector<int> run;
   run.push_back(6734);
   run.push_back(6735);

   int runPeriod   = 2;
   int probeNumber = 7;
   int axis        = 2; 
   int prMethod    = gm2fieldUtil::Constants::kPhaseDerivative;
   int ppMethod    = plungingProbeAnalysis::kLeastSquaresPhase;

   std::string cutFile = "cutData.json";

   char cutPath[200];
   sprintf(cutPath,"./input/json/run-%d/%s",runPeriod,cutFile.c_str());
   std::string cutpath = cutPath;

   std::string prodVersion   = "v9_30_00_dev";
   std::string nmrAnaVersion = "v02_02";

   double tempCorValue = -0.340;

   bool useNMRANA = true;
   int probe = 1;

   const int N = run.size();
   std::vector<plungingProbeAnaEvent_t> ppData;
   for(int i=0;i<N;i++){
      std::cout << "Getting PP data for run " << run[i] << "..." << std::endl;
      rc = GetPlungingProbeData(run[i],prMethod,ppMethod,ppData,prodVersion,nmrAnaVersion,cutpath,useNMRANA,tempCorValue);
      if(rc!=0){
         std::cout << "No data!" << std::endl;
         return 1;
      }
   }

   std::vector<int> fxprList;  
   rc = gm2fieldUtil::Import::ImportData1<int>("./input/probe-lists/fxpr-list_set-1.csv","csv",fxprList); 
   
   // time cut for the FXPR data
   bool isMax=false;
   cutpath = "./input/json/run-2/extra-cuts.json"; 
   unsigned long long t0 = GetFXPRCutTime(cutpath,probeNumber,axis,isMax); 
   
   bool subtractDrift = true;
   int period         = 10;
   std::vector<averageFixedProbeEvent_t> fxprData,fxprFltr;
   for(int i=0;i<N;i++){
      rc = GetFixedProbeData_avg(run[i],prMethod,fxprList,fxprData,prodVersion,subtractDrift,period);
      if(isMax){
	 rc = GetFixedProbeData_avg(run[i],prMethod,fxprList,fxprFltr,prodVersion,subtractDrift,period,0,t0);
      }else{
	 rc = GetFixedProbeData_avg(run[i],prMethod,fxprList,fxprFltr,prodVersion,subtractDrift,period,t0);
      }
      if(rc!=0){
         std::cout << "No data!" << std::endl;
         return 1;
      }
   }

   // Fixed probe plot 
   TGraph *gFXPR = GetFXPRTGraph_avg("GpsTimeStamp","freq","NONE",fxprData);
   gm2fieldUtil::Graph::SetGraphParameters(gFXPR,21,kBlack);

   TGraph *gFXPRf = GetFXPRTGraph_avg("GpsTimeStamp","freq","NONE",fxprFltr);
   gm2fieldUtil::Graph::SetGraphParameters(gFXPRf,20,kRed);

   // Plunging probe plots 
   TGraph *g1 = GetPPTGraph1("TimeStamp","freq",ppData);
   gm2fieldUtil::Graph::SetGraphParameters(g1,25,kBlack);

   TLegend *L = new TLegend(0.6,0.6,0.8,0.8); 
   L->AddEntry(g1,"All Data","p"); 

   TMultiGraph *mg = new TMultiGraph();
   mg->Add(g1,"lp"); 

   TMultiGraph *mgfp = new TMultiGraph();
   mgfp->Add(gFXPR ,"lp"); 
   mgfp->Add(gFXPRf,"lp"); 
  
   AddPPMultiGraph("TimeStamp","freq",ppData,mg,L);   

   const int NFP = fxprData.size(); 
   double xMin = fxprData[0].time/1E+9; 
   double xMax = fxprData[NFP-1].time/1E+9; 

   TCanvas *c1 = new TCanvas("c1","PP & FXPR Data",1200,600);
   c1->Divide(1,2); 
  
   c1->cd(1); 
   mg->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mg,"PP Data","","Frequency (Hz)"); 
   gm2fieldUtil::Graph::UseTimeDisplay(mg); 
   mg->GetXaxis()->SetLimits(xMin,xMax);  
   gm2fieldUtil::Graph::SetGraphLabelSizes(mg,0.05,0.06); 
   mg->Draw("a");
   L->Draw("same"); 
   c1->Update();

   c1->cd(2);
   mgfp->Draw("alp"); 
   gm2fieldUtil::Graph::SetGraphLabels(mgfp,"FXPR Data","","Frequency (Hz)"); 
   gm2fieldUtil::Graph::UseTimeDisplay(mgfp); 
   mgfp->GetXaxis()->SetLimits(xMin,xMax);  
   gm2fieldUtil::Graph::SetGraphLabelSizes(mgfp,0.05,0.06); 
   mgfp->Draw("alp"); 
   c1->Update(); 

   return 0;
}
//______________________________________________________________________________
int AddPPMultiGraph(std::string xAxis,std::string yAxis,std::vector<plungingProbeAnaEvent_t> data,
                    TMultiGraph *mg,TLegend *L){
   // add to a multigraph that is sorted by PP-DAQ run number
   int color=2;
   TString label; 
   const int N = data.size(); 
   for(int i=0;i<N;i++){
      do{
	 color++; 
      }while(color==10||color==19); 
      TGraph *g = getTGraph(data[i].run,xAxis,yAxis,data);
      gm2fieldUtil::Graph::SetGraphParameters(g,20,2+i);
      label = Form("PP-DAQ run %d",data[i].run);  
      mg->Add(g,"p"); 
      L->AddEntry(g,label,"p"); 
   }
   return 0;
}
//______________________________________________________________________________
TGraph *getTGraph(int run,std::string xAxis,std::string yAxis,std::vector<plungingProbeAnaEvent_t> data){
   std::vector<double> x,y;
   FillPPVector3(run,xAxis,data,x); 
   FillPPVector3(run,yAxis,data,y);
   TGraph *g = gm2fieldUtil::Graph::GetTGraph(x,y);
   return g;  
}
//______________________________________________________________________________
int FillPPVector3(int ppDAQ_run,std::string axis,std::vector<plungingProbeAnaEvent_t> data,
                  std::vector<double> &x){

   // find the index for ppDAQ_run 
   int k=0; 
   const int N = data.size();
   for(int i=0;i<N;i++) if(data[i].run==ppDAQ_run) k = i;

   const int M = data[k].numTraces; 
   for(int j=0;j<M;j++){
      if(axis=="event")      x.push_back( (double)(j+1)       );
      if(axis=="TimeStamp")  x.push_back(data[k].time[j]/1E+9 );
      if(axis=="x")          x.push_back(data[k].r[j]         );
      if(axis=="y")          x.push_back(data[k].y[j]         );
      if(axis=="z")          x.push_back(data[k].phi[j]       );
      if(axis=="temp")       x.push_back(data[k].temp[j]      );
      if(axis=="temp_err")   x.push_back(data[k].temp_err[j]  );
      if(axis=="freq")       x.push_back(data[k].freq[j]      );
      if(axis=="freq_err")   x.push_back(data[k].freq_err[j]  );
      if(axis=="freq_LO")    x.push_back(data[k].freq_LO[j]   );
      if(axis=="freq_RF")    x.push_back(data[k].freq_RF[j]   );
   } 
   return 0;
}
