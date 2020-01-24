// Calculate the imposed gradient  
// Use a 2D fit to extract imposed gradients for each probe   

#include <cstdlib> 
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>  
#include <cmath> 

#include "TF1.h"
#include "TF2.h"
#include "TH2F.h"
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
// #include "./src/TRLYCoordinates.C"
#include "./src/SystFuncs.C"

int PrintToFile(const char *outpath,std::vector<double> x,std::vector<double> dx);
int CalculateResidual(TGraph2D *g,TF2 *myFit,TGraph2D *gr); 

double myFitFunc2D(double *x,double *par); 
double myFitFunc2D_deriv_x(double *x,double *par);
double myFitFunc2D_deriv_y(double *x,double *par);

double myFitFunc2D_3rdOrder(double *x,double *par); 
double myFitFunc2D_3rdOrder_deriv_x(double *x,double *par);
double myFitFunc2D_3rdOrder_deriv_y(double *x,double *par);
double myFitFunc2D_3rdOrder_deriv_x_err(double *x,double *parErr); 
double myFitFunc2D_3rdOrder_deriv_y_err(double *x,double *parErr); 

int ImposedGrad_xy_2D_prod(std::string configFile){
 
   int rc=0;
 
   std::cout << "------------------------------------" << std::endl;
   std::cout << "GET IMPOSED GRADIENT (USING TRLY DELTA-B VALUES, 2D FIT)" << std::endl;

   InputManager *inputMgr = new InputManager();
   inputMgr->UseAxis();         // need to grab the axis data in the JSON file 
   inputMgr->Load(configFile);
   inputMgr->Print();

   std::string blindLabel    = inputMgr->GetBlindLabel();

   bool isBlind              = inputMgr->IsBlind();
   int runPeriod             = inputMgr->GetRunPeriod();
   int fitOrder              = inputMgr->GetImpGradFitOrder(); 
   // systematics 
   bool isSyst               = inputMgr->GetSystStatus();
   bool varyFit              = inputMgr->GetSystFitStatus("imp-grad");
   int systDirNum            = inputMgr->GetSystDirNum();

   // date_t theDate; 
   // GetDate(theDate);
   std::string theDate = inputMgr->GetAnaDate();

   std::string plotDir = GetPath("plots" ,isBlind,blindLabel,theDate,isSyst,systDirNum);
   std::string outDir  = GetPath("output",isBlind,blindLabel,theDate,isSyst,systDirNum);

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

   int NP = 6;   // for fitOrder = 2 
   if(fitOrder==3) NP = 10;

   const int npar = NP; 
   double xMin    = -40;
   double xMax    =  40;
   double yMin    = -40;
   double yMax    =  40;
 
   double par[npar]; 
   for(int i=0;i<npar;i++) par[i] = 1.; 

   TF2 *myFit_x = NULL;
   TF2 *myFit_y = NULL;

   if(fitOrder==2){
      myFit_x = new TF2("myFit_x",myFitFunc2D         ,xMin,xMax,yMin,yMax,npar); 
      myFit_y = new TF2("myFit_y",myFitFunc2D         ,xMin,xMax,yMin,yMax,npar); 
   }else if(fitOrder==3){
      myFit_x = new TF2("myFit_x",myFitFunc2D_3rdOrder,xMin,xMax,yMin,yMax,npar); 
      myFit_y = new TF2("myFit_y",myFitFunc2D_3rdOrder,xMin,xMax,yMin,yMax,npar); 
   }else{
      std::cout << "[ImposedGrad_xy_2D_prod]: Invalid fit order " << fitOrder << std::endl;
      return 1;
   }

   for(int i=0;i<npar;i++) myFit_x->SetParameter(i,par[i]); 
   for(int i=0;i<npar;i++) myFit_y->SetParameter(i,par[i]); 

   TGraph2D *gDBx = new TGraph2D();
   TGraph2D *gDBy = new TGraph2D();

   double xc=0,yc=0;
   const int NBin = 20; 
   TH2F *hDBx = new TH2F("hDBx","#DeltaB_{x}",NBin,xMin,xMax,NBin,yMin,yMax);
   for(int i=0;i<17;i++){
      xc = GetTrolleyProbeTransverseCoordinate(i+1,"r"); 
      yc = GetTrolleyProbeTransverseCoordinate(i+1,"y"); 
      xc+=0; 
      yc+=0; 
      // std::cout << Form("probe %02d: x = %.3lf, y = %.3lf, Delta-Bx = %.3lf",i+1,xc,yc,dB[i][0]) << std::endl; 
      gDBx->SetPoint(i,xc,yc,dB[i][0]);
   }

   TH2F *hDBy = new TH2F("hDBy","#DeltaB_{y}",NBin,xMin,xMax,NBin,yMin,yMax);
   for(int i=0;i<17;i++){
      xc = GetTrolleyProbeTransverseCoordinate(i+1,"r"); 
      yc = GetTrolleyProbeTransverseCoordinate(i+1,"y");
      xc+=0; 
      yc+=0; 
      // std::cout << Form("probe %02d: x = %.3lf, y = %.3lf, Delta-By = %.3lf",i+1,xc,yc,dB[i][1]) << std::endl; 
      gDBy->SetPoint(i,xc,yc,dB[i][1]);
   }

   gDBx->SetMarkerStyle(20); 
   gDBy->SetMarkerStyle(20); 

   TGraph2D *gDBx_res = new TGraph2D();
   TGraph2D *gDBy_res = new TGraph2D();

   TString drawOption = "colz"; // "cont4"; // "pcol"; 
 
   TCanvas *c1 = new TCanvas("c1","DeltaBx",1200,600);
   c1->Divide(2,2); 

   c1->cd(1);
   gDBx->Draw(drawOption); 
   gDBx->SetTitle("#DeltaB_{x}(x,y)"); 
   gDBx->GetXaxis()->SetTitle("x (mm)"); 
   gDBx->GetXaxis()->CenterTitle(); 
   gDBx->GetYaxis()->SetTitle("y (mm)"); 
   gDBx->GetYaxis()->CenterTitle(); 
   gDBx->Draw(drawOption);
   gDBx->Fit(myFit_x,"","R");
   // myFit_x->Draw("same surf1");  
   c1->Update(); 

   c1->cd(2);
   gDBy->Draw(drawOption); 
   gDBy->SetTitle("#DeltaB_{y}(x,y)"); 
   gDBy->GetXaxis()->SetTitle("x (mm)"); 
   gDBy->GetXaxis()->CenterTitle(); 
   gDBy->GetYaxis()->SetTitle("y (mm)"); 
   gDBy->GetYaxis()->CenterTitle(); 
   gDBy->Draw(drawOption); 
   gDBy->Fit(myFit_y,"","R"); 
   // myFit_y->Draw("same surf1");  
   c1->Update();

   CalculateResidual(gDBx,myFit_x,gDBx_res);
   CalculateResidual(gDBy,myFit_y,gDBy_res);

   c1->cd(3);
   gDBx_res->Draw(drawOption); 
   gDBx_res->SetTitle("#DeltaB_{x}(x,y) - Fit"); 
   gDBx_res->GetXaxis()->SetTitle("x (mm)"); 
   gDBx_res->GetXaxis()->CenterTitle(); 
   gDBx_res->GetYaxis()->SetTitle("y (mm)"); 
   gDBx_res->GetYaxis()->CenterTitle(); 
   gDBx_res->Draw(drawOption);
   c1->Update(); 

   c1->cd(4);
   gDBy_res->Draw(drawOption); 
   gDBy_res->SetTitle("#DeltaB_{y}(x,y) - Fit"); 
   gDBy_res->GetXaxis()->SetTitle("x (mm)"); 
   gDBy_res->GetXaxis()->CenterTitle(); 
   gDBy_res->GetYaxis()->SetTitle("y (mm)"); 
   gDBy_res->GetYaxis()->CenterTitle(); 
   gDBy_res->Draw(drawOption); 
   c1->Update();

   TString plotPath = Form("%s/imposed-grad_2D-fits.png",plotDir.c_str());
   c1->cd();
   c1->Print(plotPath);
 
   std::vector<double> xPar,xParErr; 
   std::vector<double> yPar,yParErr; 
 
   for(int i=0;i<npar;i++){
      xPar.push_back( myFit_x->GetParameter(i) ); 
      xParErr.push_back( myFit_x->GetParError(i) ); 
      yPar.push_back( myFit_y->GetParameter(i) ); 
      yParErr.push_back( myFit_y->GetParError(i) ); 
   }

   double DBX=0,DBY=0,DBX_err=0,DBY_err=0;
   double X[2],parx[npar],pary[npar],parxErr[npar],paryErr[npar]; 
   for(int i=0;i<npar;i++) parx[i] = xPar[i]; 
   for(int i=0;i<npar;i++) pary[i] = yPar[i];
   for(int i=0;i<npar;i++) parxErr[i] = xParErr[i]; 
   for(int i=0;i<npar;i++) paryErr[i] = yParErr[i];

   std::vector<double> igX,igY,igXe,igYe;  
   // std::cout << "FIT RESULTS" << std::endl;
   for(int i=0;i<17;i++){
      X[0] = GetTrolleyProbeTransverseCoordinate(i+1,"r"); 
      X[1] = GetTrolleyProbeTransverseCoordinate(i+1,"y");
      if(fitOrder==2){ 
	 DBX  = myFitFunc2D_deriv_x(X,parx);  
	 DBY  = myFitFunc2D_deriv_y(X,pary);
	 DBX_err = 0;  
	 DBY_err = 0;
      }else if(fitOrder==3){
	 DBX     = myFitFunc2D_3rdOrder_deriv_x(X,parx);  
	 DBY     = myFitFunc2D_3rdOrder_deriv_y(X,pary);
	 DBX_err = myFitFunc2D_3rdOrder_deriv_x_err(X,parxErr);  
	 DBY_err = myFitFunc2D_3rdOrder_deriv_y_err(X,paryErr);
      } 
      // std::cout << Form("probe %02d: dB_I/dx = %.3lf ± %.3lf Hz/mm, dB_I/dy = %.3lf ± %.3lf Hz/mm",i+1,DBX,DBX_err,DBY,DBY_err) << std::endl;
      igX.push_back(DBX);  
      igY.push_back(DBY); 
      igXe.push_back(DBX_err);  
      igYe.push_back(DBY_err);  
   }
 
   if(isSyst && varyFit){
      std::cout << "[ImposedGrad_xy_2D_prod]: SYSTEMATIC VARIATION! Will vary fit result within (Gaussian) uncertainties" << std::endl;
      rc = systFunc::RandomizeFitValues(igX,igXe);
      rc = systFunc::RandomizeFitValues(igY,igYe);
   }

   // print to file
   char outpath_rad[200],outpath_vert[200]; 
   sprintf(outpath_rad ,"%s/imposed-grad-x_fit-pars_2D.csv",outDir.c_str()); 
   sprintf(outpath_vert,"%s/imposed-grad-y_fit-pars_2D.csv",outDir.c_str()); 

   rc = PrintToFile(outpath_rad ,xPar,xParErr); 
   rc = PrintToFile(outpath_vert,yPar,yParErr); 

   sprintf(outpath_rad ,"%s/imposed-grad-x_2D.csv",outDir.c_str()); 
   sprintf(outpath_vert,"%s/imposed-grad-y_2D.csv",outDir.c_str()); 

   rc = PrintToFile(outpath_rad ,igX,igXe); 
   rc = PrintToFile(outpath_vert,igY,igYe); 

   // fit errors 
   // sprintf(outpath_rad ,"%s/imposed-grad-x_fit-err.csv",outDir.c_str()); 
   // sprintf(outpath_vert,"%s/imposed-grad-y_fit-err.csv",outDir.c_str()); 

   // rc = PrintToFile(outpath_rad ,VER,radGrad_err); 
   // rc = PrintToFile(outpath_vert,RAD,vertGrad_err); 

   return rc;
}
//______________________________________________________________________________
int CalculateResidual(TGraph2D *g,TF2 *myFit,TGraph2D *gr){

   double *X = g->GetX(); 
   double *Y = g->GetY(); 
   double *Z = g->GetZ(); 

   int N = g->GetN(); 
   // std::cout << "RESIDUALS: " << std::endl;
   // std::cout << "NPTS = " << N << std::endl;

   double res=0,func=0,x,y,val;
   for(int i=0;i<N;i++){
      func = myFit->Eval(X[i],Y[i]); 
      res  = Z[i] - func;
      // std::cout << Form("point %02d: x = %9.3lf, y = %9.3lf, func = %9.3lf, val = %9.3lf, res = %9.3lf",i+1,X[i],Y[i],func,Z[i],res) << std::endl;
      gr->SetPoint(i,X[i],Y[i],res); 
   }
   return 0;
}
//______________________________________________________________________________
int PrintToFile(const char *outpath,std::vector<double> x,std::vector<double> dx){
   char outStr[200]; 
   const int N = x.size();
   std::ofstream outfile; 
   outfile.open(outpath); 

   if( outfile.fail() ){
      std::cout << "[ImposedGrad_xy_2D_prod]: Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      for(int i=0;i<N;i++){
	 sprintf(outStr,"%.3lf,%.3lf",x[i],dx[i]);
	 outfile << outStr << std::endl; 
      }
      outfile.close();
   }
   return 0; 
}
//______________________________________________________________________________
double myFitFunc2D(double *x,double *par){
   double X  = x[0]; 
   double Y  = x[1];
   double X2 = X*X; 
   double Y2 = Y*Y;  
   double f  = par[0] + par[1]*X + par[2]*Y + par[3]*X2 + par[4]*Y2 + par[5]*X*Y; 
   return f;
}
//______________________________________________________________________________
double myFitFunc2D_deriv_x(double *x,double *par){
   double X  = x[0]; 
   double Y  = x[1];
   double X2 = X*X;
   double Y2 = Y*Y;   
   double f  = par[1] + 2.*par[3]*X + par[5]*Y; 
   return f;
}
//______________________________________________________________________________
double myFitFunc2D_deriv_y(double *x,double *par){
   double X  = x[0]; 
   double Y  = x[1]; 
   double X2 = X*X;
   double Y2 = Y*Y;
   double f  = par[2] + 2.*par[4]*Y + par[5]*X; 
   return f;
}
//______________________________________________________________________________
double myFitFunc2D_3rdOrder(double *x,double *par){
   double X = x[0]; 
   double Y = x[1]; 
   double X2 = X*X;
   double X3 = X*X*X; 
   double Y2 = Y*Y; 
   double Y3 = Y*Y*Y;
   double f  = par[0] + par[1]*X + par[2]*Y + par[3]*X2 + par[4]*Y2 + par[5]*X*Y 
             + par[6]*X3 + par[7]*Y3 + par[8]*X2*Y + par[9]*X*Y2; 
   return f;
}
//______________________________________________________________________________
double myFitFunc2D_3rdOrder_deriv_x(double *x,double *par){
   double X = x[0]; 
   double Y = x[1]; 
   double X2 = X*X;
   double X3 = X*X*X; 
   double Y2 = Y*Y; 
   double Y3 = Y*Y*Y; 
   double f  = par[1] + 2.*par[3]*X + par[5]*Y + 3.*par[6]*X2 + 2.*par[8]*X*Y + par[9]*Y2;
   return f;
}
//______________________________________________________________________________
double myFitFunc2D_3rdOrder_deriv_y(double *x,double *par){
   double X = x[0]; 
   double Y = x[1]; 
   double X2 = X*X;
   double X3 = X*X*X; 
   double Y2 = Y*Y; 
   double Y3 = Y*Y*Y; 
   double f  = par[2] + 2.*par[4]*Y + par[5]*X + 3.*par[7]*Y2 + par[8]*X2 + 2.*par[9]*X*Y; 
   return f;
}
//______________________________________________________________________________
double myFitFunc2D_3rdOrder_deriv_x_err(double *x,double *parErr){
   double X = x[0]; 
   double Y = x[1]; 
   const int npar = 6;
   int index[npar]    = {1,3,5,6,8,9};
   double coeff[npar] = {1.,2.,1.,3.,2.,1.}; 
   double err=0,err_sq=0;
   if(X==0 && Y!=0){
      err = TMath::Sqrt( parErr[1]*parErr[1] + parErr[5]*parErr[5] + parErr[9]*parErr[9] ); 
   }else if(X!=0 && Y==0){
      err = TMath::Sqrt( parErr[1]*parErr[1] + 4.*parErr[3]*parErr[3] + 9.*parErr[6]*parErr[6] ); 
   }else if(X==0 && Y==0){
      err = parErr[1]; 
   }else{
      for(int i=0;i<npar;i++){
	 err_sq += coeff[i]*coeff[i]*parErr[index[i]]*parErr[index[i]]; 
      }
      err = TMath::Sqrt(err_sq); 
   } 
   return err;
}
//______________________________________________________________________________
double myFitFunc2D_3rdOrder_deriv_y_err(double *x,double *parErr){
   double X = x[0]; 
   double Y = x[1]; 
   const int npar = 6;
   int index[npar]    = {2,4,5,7,8,9};
   double coeff[npar] = {1.,2.,1.,3.,1.,2.}; 
   double err=0,err_sq=0;
   if(X==0 && Y!=0){
      err = TMath::Sqrt( parErr[2]*parErr[2] + 4.*parErr[4]*parErr[4] + 9.*parErr[7]*parErr[7] ); 
   }else if(X!=0 && Y==0){
      err = TMath::Sqrt( parErr[2]*parErr[2] + parErr[5]*parErr[5] + parErr[8]*parErr[8] ); 
   }else if(X==0 && Y==0){
      err = parErr[2]; 
   }else{
      for(int i=0;i<npar;i++){
	 err_sq += coeff[i]*coeff[i]*parErr[index[i]]*parErr[index[i]]; 
      }
      err = TMath::Sqrt(err_sq); 
   }
   return err;
}
