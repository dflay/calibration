// compare calibration analysis results based on various settings  

#include <cstdlib> 
#include <iostream>
#include <fstream> 

#include "TPad.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1F.h"

#include "gm2fieldMath.h"
#include "gm2fieldImport.h"
#include "gm2fieldExport.h"
#include "gm2fieldGraph.h"

#include "./include/Color.h"

TGraphErrors *GetTGraphErrors(std::string xAxis,std::string yAxis,std::string yAxisErr,json data);

TMultiGraph *GetTMultiGraph_diff(std::vector< std::vector<double> > X,
                                 std::vector< std::vector<double> > Y,std::vector< std::vector<double> > EY,
                                 std::vector<std::string> label,double sf,double rho,TH1F *h); 

int PrintToScreen(std::string label,std::string axis,std::string axisErr,json data); 
int PrintToFile_csv(const char *outpath,std::string xAxis,std::string yAxis,std::string yAxisErr,json data); 

int FillVector(std::string axis,json data,std::vector<double> &v); 
int AddToMultiGraph(int color,int marker,std::string label,std::string xAxis,std::string yAxis,std::string yAxisErr,
                    json data,TMultiGraph *mg,TLegend *L); 

int CompareResults(){

   json input; 

   std::string inpath = "./input/json/compare.json"; 
   int rc = gm2fieldUtil::Import::ImportJSON(inpath,input);

   const int NT = input["ana-set"].size();

   bool isBlind           = (bool)input["blinding"]["enable"];
   bool plotBand          = (bool)input["blinding"]["plot-band"];  
   std::string blindLabel = input["blinding"]["label"]; 
   
   std::string date="NONE";
   const int ND = input["date"].size(); 
   if(ND==0) date = input["date"][0];  
   
   bool plotSim = (bool)( (int)input["plot-sim"] ); 

   std::string bl_prefix; 
   if(isBlind){
      bl_prefix = "./output/blinded/" + blindLabel; 
   }else{
      bl_prefix = "./output/unblinded"; 
   }

   std::string xAxis    = "probe";
   std::string yAxis    = input["y-axis"]; 
   std::string yAxisErr = input["y-axis-err"];

   // load simulation data 
   std::vector<double> sx,sy,sye;
   std::string inpath_sim = "./output/simulation/osc-change.csv"; 
   rc = gm2fieldUtil::Import::ImportData2<double,double>(inpath_sim,"csv",sy,sye);
   const int NS = sy.size();
   for(int i=0;i<NS;i++) sx.push_back(i+1);

   TGraphErrors *gSim = gm2fieldUtil::Graph::GetTGraphErrors(sx,sy,sye);
   gm2fieldUtil::Graph::SetGraphParameters(gSim,21,kGreen+2); 

   TString Title      = Form("Calibration Results"); 
   TString xAxisTitle = Form("%s",xAxis.c_str());
   TString yAxisTitle = Form("%s (Hz)",yAxis.c_str());

   TMultiGraph *mg = new TMultiGraph(); 
   TLegend *L = new TLegend(0.6,0.6,0.8,0.8); 

   std::cout << "Comparing " << NT << " data sets..." << std::endl; 

   std::vector<double> x,y,ey;  
   std::vector< std::vector<double> > X,Y,EY; 
 
   char filename[200],alabel[200];
   std::string dateStr,anaStr; 
   std::vector<std::string> label; 
   int mkColor=0,marker=19;
   json result; 
   for(int i=0;i<NT;i++){
      mkColor = color_df::color[i];
      marker++;
      if(ND>1) date = input["date"][i];
      dateStr = input["date"][i]; 
      anaStr  = input["ana-set"][i]; 
      sprintf(alabel,"%s, %s",anaStr.c_str(),dateStr.c_str() );  
      label.push_back(alabel); 
      sprintf(filename,"%s_%s.csv",anaStr.c_str(),dateStr.c_str() ); 
      inpath = bl_prefix + "/" + date + "/" + anaStr + "/calibData_" + date + ".json"; 
      std::cout << "Processing data set: " << label[i] << std::endl; 
      rc = gm2fieldUtil::Import::ImportJSON(inpath,result);
      AddToMultiGraph(mkColor,marker,label[i],xAxis,yAxis,yAxisErr,result,mg,L);
      rc = PrintToFile_csv(filename,xAxis,yAxis,yAxisErr,result);
      rc = PrintToScreen(label[i],yAxis,yAxisErr,result); 
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

   // load Ran's data
   json rData; 
   std::string ran_fn     = input["compare-ran"]["filename"];  
   std::string inpath_ran = "./output/blinded/bingzhi-li/" + ran_fn;
   rc = gm2fieldUtil::Import::ImportJSON(inpath_ran,rData);
   if(rc!=0) return 1; 

   // load Bingzhi's data
   json bData;  
   std::string bing_fn     = input["compare-bingzhi"]["filename"];  
   std::string inpath_bing = "./output/blinded/bingzhi-li/" + bing_fn; 
   rc = gm2fieldUtil::Import::ImportJSON(inpath_bing,bData);
   if(rc!=0) return 1; 

   int NL = label.size();

   bool plotRan  = (bool)( (int)input["compare-ran"]["enable"] );
   bool plotBing = (bool)( (int)input["compare-bingzhi"]["enable"] );

   if(plotRan){
      label.push_back("Ran"); 
      NL = label.size(); 
      std::cout << "Adding Ran's data to the plots..." << std::endl;
      AddToMultiGraph(kOrange+1,20,label[NL-1],xAxis,yAxis,yAxisErr,rData,mg,L);
      if(yAxis.compare("calibCoeff_opt")==0){
	 yAxis    = "calibCoeff_cor"; 
	 yAxisErr = "calibCoeffErr_cor"; 
      }
      // also gather all data we care about into three doubly-indexed vectors.  
      FillVector(xAxis,rData,x); 
      FillVector(yAxis,rData,y); 
      FillVector(yAxisErr,rData,ey); 
      X.push_back(x);  
      Y.push_back(y);  
      EY.push_back(ey); 
      std::cout << "--> Done!" << std::endl;
   } 

   x.clear();
   y.clear();
   ey.clear();

   yAxis    = "calibCoeff_opt"; 
   yAxisErr = "calibCoeffErr_opt"; 

   if(plotBing){
      label.push_back("Bingzhi"); 
      NL = label.size(); 
      std::cout << "Adding Bingzhi's data to the plots..." << std::endl;
      AddToMultiGraph(kCyan+1,20,label[NL-1],xAxis,yAxis,yAxisErr,bData,mg,L);
      if(yAxis.compare("calibCoeff_opt")==0){
	 yAxis    = "calibCoeff"; 
	 yAxisErr = "calibCoeffErr"; 
      }
      // also gather all data we care about into three doubly-indexed vectors.  
      FillVector(xAxis,bData,x); 
      FillVector(yAxis,bData,y); 
      FillVector(yAxisErr,bData,ey); 
      X.push_back(x);  
      Y.push_back(y);  
      EY.push_back(ey); 
      std::cout << "--> Done!" << std::endl;
   } 

   double sf  = 1.; 
   double rho = input["correlation"]; 
   TH1F *h = new TH1F("hDiff","hDiff",100,-10,10);  
   TMultiGraph *mgDiff = GetTMultiGraph_diff(X,Y,EY,label,sf,rho,h);
   if(plotSim) mgDiff->Add(gSim,"lp");  
   
   TLegend *LD = new TLegend(0.6,0.6,0.8,0.8);
   if(plotSim) LD->AddEntry(gSim,"Simulation"); 

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
   if(plotSim) LD->Draw("same"); 
   c1->Update();
 
   TCanvas *c2 = new TCanvas("c2","Calibration Comparison Histo",1200,600);
   c2->cd();

   h->Draw("");
   c2->Update();  

   return 0;
}
//______________________________________________________________________________
TMultiGraph *GetTMultiGraph_diff(std::vector< std::vector<double> > X,
                                 std::vector< std::vector<double> > Y,std::vector< std::vector<double> > EY,
				 std::vector<std::string> label,double sf,double rho,
                                 TH1F *h){

   TMultiGraph *mg = new TMultiGraph(); 

   // now get differences
   const int N = X.size(); 
   int M=0,mkColor=0;
   double arg=0,arg_err=0; 
   std::vector<double> x,diff,diff_err; 
   std::vector< std::vector<double> > DIFF,DIFF_ERR;
   for(int i=1;i<N;i++){
      M = X[0].size();
      for(int j=0;j<M;j++){
	 arg     = (Y[i][j] - Y[0][j])/sf; 
         arg_err = TMath::Sqrt( EY[i][j]*EY[i][j] + EY[0][j]*EY[0][j] - 2.*rho*EY[i][j]*EY[0][j])/sf; // estimate
	 if(label[i].compare("Ran")==0 || label[i].compare("Bingzhi")==0 ){
	    arg     = (Y[i][j] - Y[0][j])/sf; 
	    arg_err = TMath::Sqrt( EY[i][j]*EY[i][j] + EY[0][j]*EY[0][j] - 2.*rho*EY[i][j]*EY[0][j])/sf; // estimate
	 }
         // arg_err = TMath::Sqrt( EY[i][j]*EY[i][j] + EY[0][j]*EY[0][j] )/sf; // estimate
         // arg_err = 0.5*( EY[i][j] + EY[0][j] )/sf; // estimate.  This seems more realistic
	 // std::cout << X[0][j] << " " << arg << " " << arg_err << std::endl;
	 x.push_back(X[0][j]);  
         diff.push_back(arg); 
         diff_err.push_back(arg_err);
	 // fill histogram 
	 h->Fill(arg); 
      }
      mkColor = color_df::color[i]; 
      if(label[i].compare("Ran")==0 )     mkColor = kOrange+1; 
      if(label[i].compare("Bingzhi")==0 ) mkColor = kCyan+1; 
      TGraphErrors *g = gm2fieldUtil::Graph::GetTGraphErrors(x,diff,diff_err);
      gm2fieldUtil::Graph::SetGraphParameters(g,20,mkColor); 
      mg->Add(g,"lp");
      // set up for next data 
      x.clear(); 
      diff.clear();
      diff_err.clear(); 
   }
   return mg; 
}
//______________________________________________________________________________
int AddToMultiGraph(int color,int marker,std::string label,std::string xAxis,std::string yAxis,std::string yAxisErr,
                    json data,TMultiGraph *mg,TLegend *L){

   std::cout << label << std::endl;
   if( (label.compare("Ran")==0||label.compare("Bingzhi")==0) && yAxis.compare("calibCoeff_opt")==0 ){
      // Ran's numbers aren't necessarily ABA, but do have drift
      yAxis    = "calibCoeff_cor"; 
      yAxisErr = "calibCoeffErr_cor";
      std::cout << Form("[AddToMultiGraph]: Adding %s's data, using axes = %s, %s",
                        label.c_str(),yAxis.c_str(),yAxisErr.c_str()) << std::endl; 
   }

   TGraphErrors *g = GetTGraphErrors(xAxis,yAxis,yAxisErr,data);
   gm2fieldUtil::Graph::SetGraphParameters(g,marker,color); 

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
   const int N = data["calibCoeff"].size();
   if(N==0){
      std::cout << "[FillVector]: No data for axis = " << axis << "! " << std::endl;
      return 1;
   }
   // std::cout << "[FillVector]: Looping over " << N << " entries for axis " << axis << std::endl;
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
int PrintToScreen(std::string label,std::string axis,std::string axisErr,json data){
   std::cout << Form("%s for data set %s",axis.c_str(),label.c_str()) << std::endl;
   double x=0,dx1=0,dx2=0,dx3=0,dx4=0,dx5=0;
   const int N = data[axis].size();
   for(int i=0;i<N;i++){
      x   = data[axis][i]; 
      dx1 = data[axisErr][i];
      dx2 = data["misCor_err"][i];    
      dx3 = data["freeErr"][i];    
      dx4 = data["systErr"][i];   
      dx5 = TMath::Sqrt( dx1*dx1 + dx2*dx2 + dx3*dx3 + dx4*dx4);  
      std::cout << Form("%02d: %.3lf ± %.3lf ± %.3lf ± %.3lf ± %.3lf (tot = %.3lf)",
                        i+1,x,dx1,dx2,dx3,dx4,dx5) << std::endl; 
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
         // std::cout << outStr << std::endl;
	 outfile << outStr << std::endl;      
      }
      outfile.close();
   }
   return 0;
}
