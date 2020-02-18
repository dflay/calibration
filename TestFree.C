// test free proton code 

#include <cstdlib>
#include <iostream>
#include <vector>
#include <string> 

#include "TCanvas.h"

#include "gm2fieldGraph.h"

#include "./src/FreeProton.C"

int TestFree(){

   FreeProton *fp = new FreeProton("PP-145-01",1);
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
   gDT->Draw("alp"); 
   gm2fieldUtil::Graph::SetGraphLabels(gDT,"Total Correction","Temperature (#circC)","#delta_{t} (ppb)"); 
   gm2fieldUtil::Graph::SetGraphLabelSizes(gDT,xSize,ySize);  
   gDT->Draw("alp");
   c1->Update();

   delete fp; 

   return 0;
}
