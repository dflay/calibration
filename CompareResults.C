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

TMultiGraph *GetTMultiGraph_diff(std::vector< std::vector<double> > X,
                                 std::vector< std::vector<double> > Y,std::vector< std::vector<double> > EY,double sf=1.); 

int PrintToFile_csv(const char *outpath,std::string xAxis,std::string yAxis,std::string yAxisErr,json data); 

int FillVector(std::string axis,json data,std::vector<double> &v); 
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
   std::string yAxis    = input["y-axis"]; 
   std::string yAxisErr = input["y-axis-err"]; 

   TString Title      = Form("Calibration Results"); 
   TString xAxisTitle = Form("%s",xAxis.c_str());
   TString yAxisTitle = Form("%s (Hz)",yAxis.c_str());

   TMultiGraph *mg = new TMultiGraph(); 
   TLegend *L = new TLegend(0.6,0.6,0.8,0.8); 

   std::cout << "Comparing " << NT << " data sets..." << std::endl; 

   std::vector<double> x,y,ey;  
   std::vector< std::vector<double> > X,Y,EY; 
 
   std::string label;  
   int color=0;
   json result; 
   for(int i=0;i<NT;i++){
      color++;
      label  = input["ana-set"][i]; 
      inpath = prefix + label + "/calibData_" + date + ".json"; 
      std::cout << "Processing data set: " << label << std::endl; 
      rc = gm2fieldUtil::Import::ImportJSON(inpath,result);
      AddToMultiGraph(color,label,xAxis,yAxis,yAxisErr,result,mg,L);
      // also gather all data we care about into three doubly-indexed vectors.  
      FillVector(xAxis,result,x); 
      FillVector(yAxis,result,y); 
      FillVector(yAxisErr,result,ey); 
      X.push_back(x);  
      Y.push_back(y);  
      EY.push_back(ey); 
      // set up for next data set 
      result.clear();
      x.clear();
      y.clear();
      ey.clear();
   }   

   double sf = 1; 
   TMultiGraph *mgDiff = GetTMultiGraph_diff(X,Y,EY,sf); 

   TCanvas *c1 = new TCanvas("c1","Calibration Comparison",1200,600);
   c1->Divide(1,2); 
   
   c1->cd(1);
   mg->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mg,Title,xAxisTitle,yAxisTitle); 
   gm2fieldUtil::Graph::SetGraphLabelSizes(mg,0.05,0.06); 
   mg->Draw("a");
   L->Draw("same");
   c1->Update(); 

   c1->cd(2);
   mgDiff->Draw("a");
   gm2fieldUtil::Graph::SetGraphLabels(mgDiff,"Differences",xAxisTitle,"Difference (Hz)");
   gm2fieldUtil::Graph::SetGraphLabelSizes(mgDiff,0.05,0.06); 
   mgDiff->GetYaxis()->SetRangeUser(-10,10);  
   mgDiff->Draw("a");
   c1->Update(); 

   return 0;
}
//______________________________________________________________________________
TMultiGraph *GetTMultiGraph_diff(std::vector< std::vector<double> > X,
                                 std::vector< std::vector<double> > Y,std::vector< std::vector<double> > EY,
                                 double sf){

    TMultiGraph *mg = new TMultiGraph(); 

   // now get differences
   const int N = X.size(); 
   int M=0;
   double arg=0,arg_err=0; 
   std::vector<double> x,diff,diff_err; 
   std::vector< std::vector<double> > DIFF,DIFF_ERR;
   for(int i=1;i<N;i++){
      M = X[0].size();
      for(int j=0;j<M;j++){
	 arg     = (Y[i][j] - Y[0][j])/sf; 
         // arg_err = TMath::Sqrt( EY[i][j]*EY[i][j] + EY[0][j]*EY[0][j] )/sf; // estimate
         arg_err = 0.5*( EY[i][j] + EY[0][j] )/sf; // estimate.  This seems more realistic
	 // std::cout << X[0][j] << " " << arg << " " << arg_err << std::endl;
	 x.push_back(X[0][j]);  
         diff.push_back(arg); 
         diff_err.push_back(arg_err); 
      }
      TGraphErrors *g = gm2fieldUtil::Graph::GetTGraphErrors(x,diff,diff_err);
      gm2fieldUtil::Graph::SetGraphParameters(g,20,i); 
      mg->Add(g,"lp"); 
      // set up for next data 
      x.clear(); 
      diff.clear();
      diff_err.clear(); 
   }
   return mg; 
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
   FillVector(xAxis,data,x);  
   FillVector(yAxis,data,y);  
   FillVector(yAxisErr,data,ey);  
   TGraphErrors *g = gm2fieldUtil::Graph::GetTGraphErrors(x,y,ey);
   return g; 
}
//______________________________________________________________________________
int FillVector(std::string axis,json data,std::vector<double> &v){
   const int N = data["calibCoeff_aba"].size();
   if(N==0){
      std::cout << "[FillVector]: No data! " << std::endl;
      return 1;
   }
   for(int i=0;i<N;i++){
      if(axis.compare("probe")==0){
	 v.push_back(i+1); 
      }else{
	 v.push_back( data[axis][i] ); 
      } 
   }
   return 0;
}
//______________________________________________________________________________
int PrintToFile_csv(const char *outpath,std::string xAxis,std::string yAxis,std::string yAxisErr,json data){

   char outStr[200];
   double x=0,y=0,ey=0; 
   const int N = data[yAxis].size();

   std::ofstream outfile;
   outfile.open(outpath);

   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      for(int i=0;i<N;i++){
         if(xAxis.compare("probe")==0){
	    x = i+1;
         }else{
	    x = data[xAxis][i]; 
         } 
         y  = data[yAxis][i]; 
         ey = data[yAxisErr][i]; 
	 sprintf(outStr,"%02d,%.3lf,%.3lf",(int)x,y,ey);
	 outfile << outStr << std::endl;      
      }
      outfile.close();
   }
   return 0;
}
