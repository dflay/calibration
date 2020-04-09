// test free proton code 

#include <cstdlib>
#include <iostream>
#include <vector>
#include <string> 

#include "TCanvas.h"

#include "gm2fieldGraph.h"
#include "gm2fieldImport.h"

#include "./src/FreeProton.C"

int FillVector(int probe,std::string type,json data,std::vector<double> &x); 
TGraphErrors *GetTGraphErrors(int probe,std::string xAxis,std::string yAxis,std::string yAxisErr,json data); 

int TestFree(){

   std::string probeID = "PP-145-01";
   int runPeriod       = 1;  

   FreeProton *fp = new FreeProton(probeID,runPeriod);
   fp->Print(); 



   const int NPTS = 10; 
   double min = 19;
   double max = 30; 
   double step = (max-min)/( (double)NPTS ); 

   double it=0,arg=0,arg_err=0;
   std::vector<double> T,chi,delta_b,delta_t,sigma; 
   std::vector<double> chi_err,delta_b_err,delta_t_err,sigma_err; 

   double freq = 61792272.689; 
   double temp = 25.0; 
   double freqFree    = fp->GetOmegaP_free(freq,temp);
   double DELTA_B     = fp->GetDelta_b(temp)/1E-9; 
   double DELTA_B_ERR = fp->GetDelta_b_err(temp)/1E-9; 

   std::cout << Form("T = %.1lf deg C, freq = %.3lf Hz, freqFree = %.3lf Hz, free - raw = %.3lf Hz, delta_b = %.1lf ± %.1lf",
                temp,freq,freqFree,freqFree-freq,DELTA_B,DELTA_B_ERR) << std::endl; 

   for(int i=0;i<NPTS;i++){
      it      = min + ( (double)i )*step;
      T.push_back(it);
      arg     = fp->GetDelta_b(it)/1E-9; 
      arg_err = fp->GetDelta_b_err(it)/1E-9; 
      delta_b.push_back(arg); 
      delta_b_err.push_back(arg_err); 
      arg     = fp->GetDelta_t(it)/1E-9; 
      arg_err = fp->GetDelta_t_err(it)/1E-9;   
      delta_t.push_back(arg); 
      delta_t_err.push_back(arg_err);
      arg     = fp->GetChi(it)/1E-9; 
      arg_err = fp->GetChi_err(it)/1E-9;
      chi.push_back(arg); 
      chi_err.push_back(arg_err);
      arg     = fp->GetSigma(it)/1E-9; 
      arg_err = fp->GetSigma_err(it)/1E-9;
      std::cout << Form("T = %.1lf, sigma = %.1lf ± %.1lf ppb, delta_b = %.1lf ± %.1lf ppb, delta_t = %.1lf ± %.1lf ppb",
                        it,arg,arg_err,delta_b[i],delta_b_err[i],delta_t[i],delta_t_err[i]) << std::endl;
      sigma.push_back(arg); 
      sigma_err.push_back(arg_err); 
   }
  
   TGraphErrors *gDB  = gm2fieldUtil::Graph::GetTGraphErrors(T,delta_b,delta_b_err);  
   TGraphErrors *gDT  = gm2fieldUtil::Graph::GetTGraphErrors(T,delta_t,delta_t_err);  
   TGraphErrors *gCHI = gm2fieldUtil::Graph::GetTGraphErrors(T,chi,chi_err);  
   TGraphErrors *gSIG = gm2fieldUtil::Graph::GetTGraphErrors(T,sigma,sigma_err);  

   gm2fieldUtil::Graph::SetGraphParameters(gDB ,20,kBlack); 
   gm2fieldUtil::Graph::SetGraphParameters(gDT ,20,kBlack); 
   gm2fieldUtil::Graph::SetGraphParameters(gCHI,20,kBlack); 
   gm2fieldUtil::Graph::SetGraphParameters(gSIG,20,kBlack); 
  
   // load free proton, TRLY data from calibration analysis 
   json ccData;
   std::string inpath_data = "./output/blinded/flay/04-06-20/run-1/calibData_04-06-20.json"; 
   int rc = gm2fieldUtil::Import::ImportJSON(inpath_data,ccData);

   TLegend *L = new TLegend(0.6,0.6,0.8,0.8);
   TMultiGraph *mg = new TMultiGraph(); 
   mg->Add(gDT,"lp"); 

   // build graphs
   int color=0,color_last=2; 
   const int NP = 17;  
   TGraphErrors **gpp = new TGraphErrors*[NP];
   for(int i=0;i<NP;i++){
      color = color_last + i;
      if(color==10 || color==19) color++;
      gpp[i] = GetTGraphErrors(i+1,"ppTemp","fpCor","fpCorErr",ccData); 
      gm2fieldUtil::Graph::SetGraphParameters(gpp[i],21,color);
      color = color_last;
      mg->Add(gpp[i],"p"); 
      L->AddEntry(gpp[i],Form("Probe %02d",i+1),"p"); 
   }
    
   TGraphErrors *gTR  = GetTGraphErrors(-1,"probe","trTemp","trTempErr",ccData);
   gm2fieldUtil::Graph::SetGraphParameters(gTR,20,kBlack);  
 
   double xSize = 0.05; 
   double ySize = 0.06; 

   TCanvas *c1 = new TCanvas("c1","Free Proton Terms",1200,600);
   c1->Divide(2,2); 

   c1->cd(1);
   gDB->Draw("alp"); 
   gm2fieldUtil::Graph::SetGraphLabels(gDB,"Bulk Magnetic Susceptibility","Temperature (#circC)","#delta_{b} (ppb)"); 
   gm2fieldUtil::Graph::SetGraphLabelSizes(gDB,xSize,ySize);  
   gDB->Draw("alp"); 
   c1->Update();

   c1->cd(2);
   gSIG->Draw("alp"); 
   gm2fieldUtil::Graph::SetGraphLabels(gSIG,"Water Diamagnetic Shielding","Temperature (#circC)","#sigma (ppb)"); 
   gm2fieldUtil::Graph::SetGraphLabelSizes(gSIG,xSize,ySize);  
   gSIG->Draw("alp"); 
   c1->Update();

   c1->cd(3);
   gCHI->Draw("alp"); 
   gm2fieldUtil::Graph::SetGraphLabels(gCHI,"Water Magnetic Susceptibility","Temperature (#circC)","#chi (ppb)"); 
   gm2fieldUtil::Graph::SetGraphLabelSizes(gCHI,xSize,ySize);  
   gCHI->Draw("alp"); 
   c1->Update();

   c1->cd(4);
   mg->Draw("a"); 
   gm2fieldUtil::Graph::SetGraphLabels(mg,"Total Correction","Temperature (#circC)","#delta_{t} (ppb)"); 
   gm2fieldUtil::Graph::SetGraphLabelSizes(mg,xSize,ySize);  
   mg->Draw("a");
   L->Draw("same");
   c1->Update();

   TCanvas *c2 = new TCanvas("c2","TRLY Temperatures",1200,600);

   c2->cd(1);
   gTR->Draw("alp"); 
   gm2fieldUtil::Graph::SetGraphLabels(gTR,"TRLY Temperatures","Probe","Temperature (#circC)"); 
   // gm2fieldUtil::Graph::SetGraphLabelSizes(gTR,xSize,ySize);  
   gTR->Draw("alp"); 
   c2->Update();

   delete fp; 

   return 0;
}
//______________________________________________________________________________
TGraphErrors *GetTGraphErrors(int probe,std::string xAxis,std::string yAxis,std::string yAxisErr,json data){
   // get a TGraph of a given data type for a given probe number
   std::vector<double> x,y,ey; 
   FillVector(probe,xAxis   ,data,x);  
   FillVector(probe,yAxis   ,data,y);  
   FillVector(probe,yAxisErr,data,ey);  

   TGraphErrors *g = gm2fieldUtil::Graph::GetTGraphErrors(x,y,ey);
   return g; 
}
//______________________________________________________________________________
int FillVector(int probe,std::string type,json data,std::vector<double> &x){
   double arg=0;
   int N=0;
   if(type.compare("probe")==0){
      N = 17; 
      for(int i=0;i<N;i++) x.push_back(i+1);
   }else{
      N = data[type].size();
      for(int i=0;i<N;i++){
	 arg = (double)data[type][i];
	 if(probe>0){ 
	    if( probe==(i+1) ) x.push_back(arg);
	 }else{
	    // don't care about probe number
	    x.push_back(arg);
	 }
      }
   }
   return 0;
}
