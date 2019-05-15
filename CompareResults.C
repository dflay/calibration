// compare calibration analysis results based on production tags 

#include <cstdlib> 
#include <iostream>
#include <fstream> 

#include "TPad.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "gm2fieldMath.h"
#include "gm2fieldImport.h"
#include "gm2fieldExport.h"
#include "gm2fieldGraph.h"

TGraphErrors *GetTGraphErrors(std::string xAxis,std::string yAxis,std::string yAxisErr,json data); 

int AddToMultiGraph(int color,std::string label,std::string xAxis,std::string yAxis,std::string yAxisErr,
                    json data,TMultiGraph *mg,TLegend *L); 

int CompareResults(){

   json input; 

   std::string inpath = "compare.json"; 
   int rc = gm2fieldUtil::Import::ImportJSON(inpath,input);

   const int NT = input["ana-set"].size(); 

   bool isBlind = (bool)input["blinding"]["enable"]; 
   std::string blindLabel = input["blinding"]["label"]; 
   std::string date       = input["date"];  

   std::string prefix; 
   if(isBlind){
      prefix = "./output/blinded/" + blindLabel + "/" + date + "/"; 
   }else{
      prefix = "./output/"; 
   }

   std::string xAxis    = "probe"; 
   std::string yAxis    = "calibCoeff_aba"; 
   std::string yAxisErr = "calibCoeffErr_aba"; 

   TString Title      = Form("Calibration Results"); 
   TString xAxisTitle = Form("%s",xAxis.c_str());
   TString yAxisTitle = Form("%s (Hz)",yAxis.c_str());

   TMultiGraph *mg = new TMultiGraph(); 
   TLegend *L = new TLegend(0.6,0.6,0.8,0.8); 

   std::cout << "Comparing " << NT << " data sets..." << std::endl; 

   std::string label;  
   int color=0;
   json result; 
   for(int i=0;i<NT;i++){
      color++;
      label  = input["ana-set"][i]; 
      inpath = prefix + label + "/calibData_" + date + ".json"; 
      std::cout << "Processing " << label << std::endl; 
      rc = gm2fieldUtil::Import::ImportJSON(inpath,result);
      AddToMultiGraph(color,label,xAxis,yAxis,yAxisErr,result,mg,L); 
      result.clear();
   }   

   TCanvas *c1 = new TCanvas("c1","Calibration Comparison",1200,600);
   
   c1->cd();
   mg->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mg,Title,xAxisTitle,yAxisTitle); 
   mg->Draw("a");
   L->Draw("same"); 

   return 0;
}
//______________________________________________________________________________
int AddToMultiGraph(int color,std::string label,std::string xAxis,std::string yAxis,std::string yAxisErr,
                    json data,TMultiGraph *mg,TLegend *L){

   TGraphErrors *g = GetTGraphErrors(xAxis,yAxis,yAxisErr,data);
   gm2fieldUtil::Graph::SetGraphParameters(g,20,color); 

   mg->Add(g,"lp"); 
   L->AddEntry(g,label.c_str(),"p"); 

   return 0;
}
//______________________________________________________________________________
TGraphErrors *GetTGraphErrors(std::string xAxis,std::string yAxis,std::string yAxisErr,json data){

   std::vector<double> x,y,ey; 
   const int N = data[yAxis].size();
   for(int i=0;i<N;i++){
      if(xAxis.compare("probe")==0){
	 x.push_back(i+1); 
      }else{
	 x.push_back( data[xAxis][i] ); 
      } 
      y.push_back( data[yAxis][i] ); 
      ey.push_back( data[yAxisErr][i] ); 
   }

   TGraphErrors *g = gm2fieldUtil::Graph::GetTGraphErrors(x,y,ey);
   return g; 
}
