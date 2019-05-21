// Plot fixed probe data      

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
#include "./include/plungingProbeAnaEvent.h"
#include "./include/trolleyAnaEvent.h"
#include "./include/fixedProbeEvent.h"

#include "./src/FitFuncs.C"
#include "./src/TRLYFuncs.C"
#include "./src/CustomMath.C"
#include "./src/CustomImport.C"
#include "./src/CustomExport.C"
#include "./src/CustomGraph.C"
#include "./src/CustomAlgorithms.C"
#include "./src/CustomUtilities.C"

int PrintProbeInfo(int run,int method,std::string version,std::vector<int> probe); 

int PlotFXPR(){

   int rc=0;
   int method = gm2fieldUtil::Constants::kPhaseDerivative;

   int runPeriod = 1;  
   std::string prodVersion = "v9_21_01"; 

   int run=0;
   std::cout << "Enter run: ";
   std::cin  >> run; 
 
   std::vector<int> probe;
   for(int i=0;i<378;i++) probe.push_back(i); 
   // probe.push_back(201); 
   // probe.push_back(202);
   // probe.push_back(211); 
   // probe.push_back(212); 
   // probe.push_back(216); 
   // probe.push_back(217);
   // probe.push_back(226);
   // probe.push_back(227);

   rc = PrintProbeInfo(run,method,prodVersion,probe);
   return rc;

   std::vector<averageFixedProbeEvent_t> fxprData;   

   rc = GetFixedProbeData_avg(run,method,probe,fxprData,prodVersion);
   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   std::vector<fixedProbeEvent_t> fxpr1,fxpr2; 
   rc = GetFixedProbeData(run,method,probe[0],fxpr1,prodVersion); 
   rc = GetFixedProbeData(run,method,probe[1],fxpr2,prodVersion); 

   TGraphErrors *gt = GetFXPRTGraph_avg("GpsTimeStamp","freq","NONE",fxprData);
   gm2fieldUtil::Graph::SetGraphParameters(gt,20,kBlack);

   TGraph *g1 = GetFXPRTGraph("GpsTimeStamp","freq",fxpr1);
   gm2fieldUtil::Graph::SetGraphParameters(g1,20,kBlue);

   TGraph *g2 = GetFXPRTGraph("GpsTimeStamp","freq",fxpr2);
   gm2fieldUtil::Graph::SetGraphParameters(g2,20,kRed);

   TMultiGraph *mg = new TMultiGraph();
   mg->Add(g1,"lp"); 
   mg->Add(g2,"lp"); 
   mg->Add(gt,"lp"); 

   TLegend *L = new TLegend(0.6,0.6,0.8,0.8); 
   L->AddEntry(g1,Form("Probe %03d",probe[0]),"p"); 
   L->AddEntry(g2,Form("Probe %03d",probe[1]),"p"); 
   L->AddEntry(gt,"Average","p"); 

   TString Title = Form("Run %d, FXPR AVG",run);
   
   // TString plotPath = Form("tr-%02d_adev_run-%05d.png",probe,run);

   TCanvas *c1 = new TCanvas("c1","TRLY Data",1200,600);
   
   c1->cd();
   mg->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mg,"FXPR Data","","Frequency (Hz)");
   gm2fieldUtil::Graph::SetGraphLabelSizes(mg,0.05,0.06); 
   gm2fieldUtil::Graph::UseTimeDisplay(mg); 
   mg->Draw("a");
   L->Draw("same"); 
   c1->Update();

   return 0;
}
//______________________________________________________________________________
int PrintProbeInfo(int run,int method,std::string version,std::vector<int> probe){
   char outpath[200],outStr[200];
   sprintf(outpath,"probe-info.csv");

   int rc=0;
   std::vector<fixedProbeEvent_t> data;
   const int N = probe.size();
   
   ofstream outfile;
   outfile.open(outpath);

   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
   }else{
      for(int i=0;i<N;i++){
	 rc = GetFixedProbeData(run,method,probe[i],data,version);
         sprintf(outStr,"%03d,%c,%c,%c,%d,%02d",
	       probe[i],data[0].yokeID,data[0].layerID,data[0].radID,data[0].aziID,data[0].stationID);
	 std::cout << outStr << std::endl;
         outfile << outStr << std::endl; 
	 data.clear();
      }
      outfile.close();
   }

   std::cout << "--> Done" << std::endl;

   return 0;
}

