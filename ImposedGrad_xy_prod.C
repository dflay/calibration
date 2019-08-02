// Calculate the imposed gradient   

#include <cstdlib> 
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>  
#include <cmath> 

#include "TStyle.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSpline.h"

#include "gm2fieldMath.h"
#include "gm2fieldGraph.h"
#include "gm2fieldRootHelper.h"
#include "gm2fieldUnits.h"
#include "TemperatureSensor.h"
#include "MovingAverage.h"

#include "./include/Constants.h"
#include "./include/fixedProbeEvent.h"
#include "./include/trolleyAnaEvent.h"
#include "./include/sccEvent.h" 
#include "./include/date.h"

#include "./src/CustomUtilities.C"
#include "./src/CustomMath.C"
#include "./src/CustomImport.C"
#include "./src/CustomExport.C"
#include "./src/CustomAlgorithms.C"
#include "./src/CustomGraph.C"
#include "./src/InputManager.C"
#include "./src/OscFuncs.C"
#include "./src/MyFits.C"
#include "./src/FitFuncs.C"
#include "./src/FitErr.C"
#include "./src/TRLYFuncs.C"

int PrintToFile(const char *outpath,std::vector<double> x,std::vector<double> dx); 
int GetData(char axis,int i,std::vector< std::vector<double> > dB,std::vector< std::vector<double> > dB_err, 
            std::vector<int> &PR,std::vector<double> &v,std::vector<double> &w,std::vector<double> &ew); 

int ImposedGrad_xy_prod(std::string configFile){
  
   std::cout << "------------------------------------" << std::endl;
   std::cout << "GET IMPOSED GRADIENT (USING TRLY DELTA-B VALUES)" << std::endl;

   InputManager *inputMgr = new InputManager();
   inputMgr->UseAxis();         // need to grab the axis data in the JSON file 
   inputMgr->Load(configFile);
   inputMgr->Print();

   std::string blindLabel    = inputMgr->GetBlindLabel();

   bool isBlind              = inputMgr->IsBlind();
   int runPeriod             = inputMgr->GetRunPeriod();

   date_t theDate; 
   GetDate(theDate); 

   std::string plotDir = GetPath("plots" ,isBlind,blindLabel,theDate.getDateString());
   std::string outDir  = GetPath("output",isBlind,blindLabel,theDate.getDateString());

   std::vector<double> db,db_err;
   std::vector< std::vector<double> > dB,dB_err;  

   std::vector<std::string> gradName;
   gradName.push_back("rad"); 
   gradName.push_back("vert"); 

   const int NG = gradName.size(); 

   char inpath[200]; 
   std::cout << "TRLY DELTA-B VALUES: " << std::endl; 
   std::vector<deltab_t> trly;
   for(int i=0;i<17;i++){
      for(int j=0;j<NG;j++){
	 sprintf(inpath,"%s/dB-trly_final-location_%s-grad_pr-%02d.csv",outDir.c_str(),gradName[j].c_str(),i+1);
	 LoadDeltaBData(inpath,trly);
	 db.push_back( trly[j].dB_fxpr ); 
	 db_err.push_back( trly[j].dB_fxpr_err ); 
      }
      dB.push_back(db); 
      dB_err.push_back(db_err);
      std::cout << Form("probe %02d: dBx = %.3lf +/- %.3lf, dBy = %.3lf +/- %.3lf",
                        i+1,dB[i][0],dB_err[i][0],dB[i][1],dB_err[i][1]) << std::endl;
      // clean up
      db.clear(); 
      db_err.clear(); 
      trly.clear();
   }

   // use all probes in a specific height or radius
   int npar=0,MM=0,rc=0;
   const int NP = 5;
   std::vector<double> X,Y,EY; 
   double intercept=0,r=0,sum_sq=0,SLOPE=0,err=0;
   double pos[NP] = {30.3,17.5,0,-17.5,-30.3};

   TCanvas *c1 = new TCanvas("c1","Radial Gradients",600,1200); 
   c1->Divide(1,5); 
  
   TGraphErrors **gx = new TGraphErrors*[NP];  

   std::string fitFunc;
   std::vector<int> PR;
   std::vector<double> xpar,xparErr; 
   std::vector<double> VER,RG,ERG; 
   std::cout << "RADIAL GRADIENT vs HEIGHT" << std::endl;
   for(int i=0;i<NP;i++){
      // grab data from dB vector   
      GetData('x',i+1,dB,dB_err,PR,X,Y,EY);
      MM = X.size(); 
      std::cout << "----------------------------------------------------------" << std::endl; 
      std::cout << Form("DATA for h = %.3lf mm: ",pos[i]) << std::endl;
      for(int j=0;j<MM;j++) std::cout << Form("probe %02d, x = %.3lf, dBx = %.3lf +/- %.3lf",PR[j],X[j],Y[j],EY[j]) << std::endl; 
      // make a TGraph 
      gx[i] = gm2fieldUtil::Graph::GetTGraphErrors(X,Y,EY);
      gm2fieldUtil::Graph::SetGraphParameters(gx[i],20,i+1);
      c1->cd(i+1);
      gStyle->SetOptFit(111); 
      gx[i]->Draw("alp"); 
      gm2fieldUtil::Graph::SetGraphLabels(gx[i],Form("h = %.3lf mm",pos[i]),"Radius (mm)","dB_{x} (Hz)"); 
      gm2fieldUtil::Graph::SetGraphLabelSizes(gx[i],0.05,0.06); 
      gx[i]->Draw("alp");
      if(MM<=2){
	 fitFunc = "pol1"; 
      }else{
	 fitFunc = "pol2"; 
      }
      TFitResultPtr fitResult = gx[i]->Fit(fitFunc.c_str(),"QS"); 
      TF1 *myFit = gx[i]->GetFunction(fitFunc.c_str());
      npar = myFit->GetNpar();
      for(int j=0;j<npar;j++){
	 xpar.push_back( myFit->GetParameter(j) ); 
	 xparErr.push_back( myFit->GetParError(j) ); 
      }
      c1->Update();
      // linear fit  
      rc = gm2fieldUtil::Math::LeastSquaresFitting(X,Y,intercept,SLOPE,r);
      // compute errors 
      for(int j=0;j<MM;j++) sum_sq += EY[j]*EY[j]; 
      err = sqrt(sum_sq)/fabs(X[0]-X[MM-1]); // error estimate 
      // store results and print to screen
      VER.push_back(pos[i]);  
      RG.push_back(SLOPE);  
      ERG.push_back(err); 
      std::cout << Form("least squares:") << std::endl;
      std::cout << Form("pos = %.3lf mm, grad = %.3lf Hz/mm, err = %.3lf Hz/mm",VER[i],SLOPE,err) << std::endl;
      std::cout << Form("ROOT fit") << std::endl;
      std::cout << Form("pos = %.3lf mm, grad = %.3lf Hz/mm, err = %.3lf Hz/mm",VER[i],xpar[1],xparErr[1]) << std::endl;
      // clean up
      sum_sq = 0;
      PR.clear();
      X.clear();
      Y.clear();
      EY.clear();
      xpar.clear();
      xparErr.clear();
   }

   c1->cd(); 
   TString plotPath = Form("%s/rad-grad-plots.png",plotDir.c_str());
   c1->Print(plotPath); 

   TCanvas *c2 = new TCanvas("c2","Vertical Gradients",600,1200); 
   c2->Divide(1,5); 
  
   TGraphErrors **gy = new TGraphErrors*[NP];  

   std::vector<double> ypar,yparErr;
   std::vector<double> RAD,VG,EVG; 
   std::cout << "VERTICAL GRADIENT vs RADIUS" << std::endl;
   for(int i=0;i<NP;i++){
      // grab data from dB vector 
      GetData('y',i+1,dB,dB_err,PR,X,Y,EY);
      MM = X.size();
      std::cout << "----------------------------------------------------------" << std::endl; 
      std::cout << Form("DATA for r = %.3lf mm: ",pos[i]) << std::endl;
      for(int j=0;j<MM;j++) std::cout << Form("probe %02d, y = %.3lf, dBy = %.3lf +/- %.3lf",PR[j],X[j],Y[j],EY[j]) << std::endl; 
      // make a TGraph and fit  
      gy[i] = gm2fieldUtil::Graph::GetTGraphErrors(X,Y,EY);
      gm2fieldUtil::Graph::SetGraphParameters(gy[i],20,i+1);
      c2->cd(i+1);
      gStyle->SetOptFit(111); 
      gy[i]->Draw("alp"); 
      gm2fieldUtil::Graph::SetGraphLabels(gy[i],Form("r = %.3lf mm",pos[i]),"Height (mm)","dB_{y} (Hz)"); 
      gm2fieldUtil::Graph::SetGraphLabelSizes(gy[i],0.05,0.06); 
      gy[i]->Draw("alp");
      if(MM<=2){
	 fitFunc = "pol1"; 
      }else{
	 fitFunc = "pol2"; 
      }
      TFitResultPtr fitResult = gy[i]->Fit(fitFunc.c_str(),"QS"); 
      TF1 *myFit = gy[i]->GetFunction(fitFunc.c_str());
      npar = myFit->GetNpar();
      for(int j=0;j<npar;j++){
	 ypar.push_back( myFit->GetParameter(j) ); 
	 yparErr.push_back( myFit->GetParError(j) ); 
      }
      c2->Update();
      // linear fit  
      rc = gm2fieldUtil::Math::LeastSquaresFitting(X,Y,intercept,SLOPE,r);
      for(int j=0;j<MM;j++) sum_sq += EY[j]*EY[j]; 
      err = sqrt(sum_sq)/fabs(X[0]-X[MM-1]); // error estimate 
      RAD.push_back(pos[i]);  
      VG.push_back(SLOPE);  
      EVG.push_back(err); 
      std::cout << Form("least squares:") << std::endl; 
      std::cout << Form("pos = %.3lf mm, grad = %.3lf Hz/mm, err = %.3lf Hz/mm",RAD[i],SLOPE,err) << std::endl;
      std::cout << Form("ROOT fit") << std::endl;
      std::cout << Form("pos = %.3lf mm, grad = %.3lf Hz/mm, err = %.3lf Hz/mm",RAD[i],ypar[1],yparErr[1]) << std::endl;
      // clean up
      sum_sq = 0;
      PR.clear();
      X.clear();
      Y.clear();
      EY.clear();
      ypar.clear();
      yparErr.clear();
   }

   c2->cd(); 
   plotPath = Form("%s/vert-grad-plots.png",plotDir.c_str());
   c2->Print(plotPath); 

   TGraphErrors *gRadGrad_vs_Height  = gm2fieldUtil::Graph::GetTGraphErrors(VER,RG,ERG);
   gm2fieldUtil::Graph::SetGraphParameters(gRadGrad_vs_Height,20,kBlack); 

   TGraphErrors *gVertGrad_vs_Radius = gm2fieldUtil::Graph::GetTGraphErrors(RAD,VG,EVG);
   gm2fieldUtil::Graph::SetGraphParameters(gVertGrad_vs_Radius,20,kBlack); 

   TString TitleRvH = Form("Radial Gradient vs Height");
   TString TitleVvR = Form("Vertical vs Radius");

   TCanvas *c3 = new TCanvas("c3","Transverse Gradients",1200,600);
   c3->Divide(1,2);
 
   c3->cd(1); 
   gRadGrad_vs_Height->Draw("ap");
   gm2fieldUtil::Graph::SetGraphLabels(gRadGrad_vs_Height,TitleRvH,"Height Above Midplane (mm)","Gradient (Hz/mm)"); 
   gm2fieldUtil::Graph::SetGraphLabelSizes(gRadGrad_vs_Height,0.05,0.06); 
   gRadGrad_vs_Height->Draw("ap");
   TFitResultPtr fitResult_rg = gRadGrad_vs_Height->Fit("pol2","QS"); 
   TF1 *myFit_rg = gRadGrad_vs_Height->GetFunction("pol2");
   c3->Update(); 

   c3->cd(2); 
   gVertGrad_vs_Radius->Draw("ap");
   gm2fieldUtil::Graph::SetGraphLabels(gVertGrad_vs_Radius,TitleVvR,"Radius (mm)","Gradient (Hz/mm)"); 
   gm2fieldUtil::Graph::SetGraphLabelSizes(gVertGrad_vs_Radius,0.05,0.06); 
   gVertGrad_vs_Radius->Draw("ap");
   TFitResultPtr fitResult_vg = gVertGrad_vs_Radius->Fit("pol2","QS"); 
   TF1 *myFit_vg = gVertGrad_vs_Radius->GetFunction("pol2");
   c3->Update(); 

   c3->cd(); 
   plotPath = Form("%s/imposed-grad_xy-plots.png",plotDir.c_str());
   c3->Print(plotPath); 

   // get fit parameters 
   npar = myFit_rg->GetNpar();
   std::vector<double> par_rg,parErr_rg; 
   for(int i=0;i<npar;i++){
      par_rg.push_back( myFit_rg->GetParameter(i) ); 
      parErr_rg.push_back( myFit_rg->GetParError(i) ); 
   }

   npar = myFit_vg->GetNpar();
   std::vector<double> par_vg,parErr_vg; 
   for(int i=0;i<npar;i++){
      par_vg.push_back( myFit_vg->GetParameter(i) ); 
      parErr_vg.push_back( myFit_vg->GetParError(i) ); 
   }

   // determine fit error at each location
   // radial gradient as a function of height 
   double arg=0; 
   double XX[3] = {0,0,0}; 
   const int NVER = VER.size();
   std::vector<double> radGrad_err; 
   for(int i=0;i<NVER;i++){
      XX[0] = VER[i]; 
      arg   = GetFitError(myFit_rg,fitResult_rg,MyPolyFitFuncDerivative,XX); 
      radGrad_err.push_back(arg); 
   }
   // vertical gradient as a function of radius 
   const int NRAD = RAD.size(); 
   std::vector<double> vertGrad_err; 
   for(int i=0;i<NRAD;i++){
      XX[0] = RAD[i]; 
      arg   = GetFitError(myFit_vg,fitResult_vg,MyPolyFitFuncDerivative,XX); 
      vertGrad_err.push_back(arg); 
   }
 
   // print to file
   char outpath_rad[200],outpath_vert[200]; 
   sprintf(outpath_rad ,"%s/imposed-grad-x_fit-pars.csv",outDir.c_str()); 
   sprintf(outpath_vert,"%s/imposed-grad-y_fit-pars.csv",outDir.c_str()); 

   rc = PrintToFile(outpath_rad ,par_rg,parErr_rg); 
   rc = PrintToFile(outpath_vert,par_vg,parErr_vg); 

   // fit errors 
   sprintf(outpath_rad ,"%s/imposed-grad-x_fit-err.csv",outDir.c_str()); 
   sprintf(outpath_vert,"%s/imposed-grad-y_fit-err.csv",outDir.c_str()); 

   rc = PrintToFile(outpath_rad ,VER,radGrad_err); 
   rc = PrintToFile(outpath_vert,RAD,vertGrad_err); 

   return rc;
}
//______________________________________________________________________________
int PrintToFile(const char *outpath,std::vector<double> x,std::vector<double> dx){
   char header[200],outStr[200]; 
   sprintf(header,"#par,par-err"); 
   const int N = x.size();
   std::ofstream outfile; 
   outfile.open(outpath); 

   if( outfile.fail() ){
      std::cout << "[ImposedGrad_xy_prod]: Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      outfile << header << std::endl;
      for(int i=0;i<N;i++){
	 sprintf(outStr,"%.3lf,%.3lf",x[i],dx[i]);
	 outfile << outStr << std::endl; 
      }
      outfile.close();
   }
   return 0; 
}
//______________________________________________________________________________
int GetData(char axis,int i,std::vector< std::vector<double> > dB,std::vector< std::vector<double> > dB_err, 
            std::vector<int> &PR,std::vector<double> &v,std::vector<double> &w,std::vector<double> &ew){

   // axis = x, y 
   // i    = coordinate index

   if(axis=='x'){
      if(i==1){
	 // highest vertical spot
	 PR.push_back(13);
	 PR.push_back(11);
      }else if(i==2){ 
	 // second highest spot 
	 PR.push_back(14);
	 PR.push_back(4);
	 PR.push_back(10);
      }else if(i==3){
	 // center line 
	 PR.push_back(15);
	 PR.push_back(5);
	 PR.push_back(1);
	 PR.push_back(3);
	 PR.push_back(9);
      }else if(i==4){
	 // second lowest 
	 PR.push_back(16);
	 PR.push_back(2);
	 PR.push_back(8);
      }else if(i==5){
	 // lowest spot 
	 PR.push_back(17);
	 PR.push_back(7);
      }
   }else if(axis=='y'){
      if(i==1){
	 // largest radius
	 PR.push_back(8); 
	 PR.push_back(10);
      }else if(i==2){ 
	 // second highest  
	 PR.push_back(7);
	 PR.push_back(3);
	 PR.push_back(11);
      }else if(i==3){
	 // center line 
	 PR.push_back(6);
	 PR.push_back(2);
	 PR.push_back(1);
	 PR.push_back(4);
	 PR.push_back(12);
      }else if(i==4){
	 // second lowest 
	 PR.push_back(17);
	 PR.push_back(5);
	 PR.push_back(13);
      }else if(i==5){
	 // lowest radius 
	 PR.push_back(16);
	 PR.push_back(14);
      }
   }

   trolleyProbePosition_t trlyPos;
   int rc = GetTrolleyProbePositions(trlyPos);
 
   int pr=0; 
   const int NN = PR.size();
   for(int i=0;i<NN;i++){
      pr = PR[i]-1; 
      if(axis=='x'){
	 v.push_back( trlyPos.r[pr] ); 
	 w.push_back( dB[pr][0] ); 
	 ew.push_back( dB_err[pr][0] ); 
      }else if(axis=='y'){
	 v.push_back( trlyPos.y[pr] ); 
	 w.push_back( dB[pr][1] ); 
	 ew.push_back( dB_err[pr][1] ); 
      }
   } 

   return 0; 

}

