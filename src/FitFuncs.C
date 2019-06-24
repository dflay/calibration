#include "../include/FitFuncs.h"
// //______________________________________________________________________________
// TF1 *GetFitToFXPR(TString fitName,double (*fitFunc)(double *,double *),
//                   int npar,int method,std::vector<int> fxprList,std::vector<gm2field::fixedProbeFrequency_t> fxprData){
//    // fit the fixed probe data to an arbitrary function 
//    // fit assumes the zeroth parameter is the time offset 
//    // to subtract off from the time axis (ROOT has trouble with large x values) 
// 
//    const int N     = fxprData.size();
//    const int NFXPR = fxprList.size();
//    int startIndex  = fxprList[0];
//    int stopIndex   = fxprList[NFXPR-1];
// 
//    unsigned long long tMin = fxprData[0].GpsTimeStamp[startIndex]; 
//    unsigned long long tMax = fxprData[N-1].GpsTimeStamp[stopIndex]; 
// 
//    // get the average fixed probe plot 
//    TGraph *gFPAVG = GetTGraph(method,tMin,tMax,1E+9,fxprList,fxprData);
//    gm2fieldUtil::Graph::SetGraphParameters(gFPAVG,20,kBlack);
// 
//    // time offset (to get fit to work) 
//    double t0 = tMin/1E+9;
// 
//    // custom fit -- first parameter must be the t0 offset  
//    const int NPAR = npar;
//    double par[NPAR];
//    for(int i=0;i<NPAR;i++){
//       if(i==0) par[i] = t0;
//       if(i!=0) par[i] = 1.;
//    }
// 
//    double TMIN = tMin/1E+9;
//    double TMAX = tMax/1E+9;
// 
//    std::cout << Form("start time: %.0lf (%s)",TMIN,gm2fieldUtil::GetStringTimeStampFromUTC( (int)TMIN ).c_str() ) << std::endl;  
//    std::cout << Form("end time:   %.0lf (%s)",TMAX,gm2fieldUtil::GetStringTimeStampFromUTC( (int)TMAX ).c_str() ) << std::endl;  
// 
//    TF1 *myFit = new TF1(fitName,fitFunc,TMIN,TMAX,NPAR);
//    for(int i=0;i<NPAR;i++) myFit->SetParameter(i,par[i]);
//    myFit->SetLineColor(kRed);
// 
//    myFit->FixParameter(0,par[0]);
// 
//    // draw and fit
//    gROOT->SetBatch(kTRUE);  // turn off this canvas, we don't need it  
// 
//    TCanvas *c = new TCanvas("fitCanvas","Fit Canvas",500,500);
//    c->cd();
//    gFPAVG->Draw("ap");
//    gFPAVG->Fit(fitName,"Q");
//    c->Update();
// 
//    delete c;
// 
//    gROOT->SetBatch(kFALSE);  // turn canvas making back on  
// 
//    return myFit;
// }
// //______________________________________________________________________________
// TF1 *GetPolyFitToFXPR(TString fitName,int nOrder,TGraph *gFPAVG,std::vector< std::vector<double> > &eps){
//    // fit the fixed probe data to an nth-order polynomial function 
//    // fit assumes the zeroth parameter is the time offset, first parameter is order  
//    // ROOT has trouble with large x values, so we subtract off a time constant defined as 
//    // the first event in the time series   
// 
//    const int NEvents = gFPAVG->GetN(); 
//    double *x = gFPAVG->GetX();
//    double *y = gFPAVG->GetY();
//    // time offset (to get fit to work) 
//    double t0 = x[0];
// 
//    // custom fit -- first parameter must be the t0 offset  
//    const int NPAR = nOrder+2;
//    double par[NPAR];
//    for(int i=0;i<NPAR;i++){
//       if(i==0) par[i] = t0;
//       if(i==1) par[i] = nOrder;
//       if(i==2) par[i] = y[0]; 
//       if(i>1)  par[i] = 1.;
//    }
// 
//    double TMIN = t0;
//    double TMAX = x[NEvents-1];
// 
//    // std::cout << Form("start time: %.0lf (%s)",TMIN,gm2fieldUtil::GetStringTimeStampFromUTC( (int)TMIN ).c_str() ) << std::endl;  
//    // std::cout << Form("end time:   %.0lf (%s)",TMAX,gm2fieldUtil::GetStringTimeStampFromUTC( (int)TMAX ).c_str() ) << std::endl;  
// 
//    double parLo[NPAR],parHi[NPAR]; 
//    for(int i=0;i<NPAR;i++){
//       if(i<2){
// 	 parLo[i] = 0;
// 	 parHi[i] = 0;
//       }else if(i==2){ 
// 	 parLo[i] = par[i] - 1E+3;
// 	 parHi[i] = par[i] + 1E+3;
//       }else{
// 	 parLo[i] = -10;
// 	 parHi[i] =  10;
//       }
//    }
// 
//    TF1 *myFit = new TF1(fitName,MyFitFunc_poly,TMIN,TMAX,NPAR);
//    myFit->SetLineColor(kRed+1);
//    myFit->FixParameter(0,par[0]);
//    myFit->FixParameter(1,par[1]);
//    for(int i=2;i<NPAR;i++){
//       myFit->SetParameter(i,par[i]);
//       myFit->SetParLimits(i,parLo[i],parHi[i]); 
//    }
// 
//    // draw and fit
//    gROOT->SetBatch(kTRUE);  // turn off this canvas, we don't need it  
// 
//    TCanvas *c = new TCanvas("fitCanvas","Fit Canvas",500,500);
//    c->cd();
//    gFPAVG->Draw("ap");
//    TFitResultPtr fitResult = gFPAVG->Fit(fitName,"QS");
//    c->Update();
//    delete c;
//    gROOT->SetBatch(kFALSE);  // turn canvas making back on  
// 
//    std::cout << "Saving the covariance matrix..." << std::endl; 
//    TMatrixDSym myCov = fitResult->GetCovarianceMatrix();
//    // TMatrixD cor = fitResult->GetCorrelationMatrix();
//    // std::cout << "-- Covariance matrix: " << std::endl;
//    // myCov.Print(); 
// 
//    std::vector<double> COL;
//    for(int i=0;i<=nOrder+1;i++){    // row index
//       for(int j=0;j<=nOrder+1;j++){  // column index
// 	 COL.push_back(myCov[i][j]); 
// 	 // std::cout << myCov[i][j] << std::endl;
//       }
//       eps.push_back(COL);
//       COL.clear();
//    }
//    std::cout << "--> Done." << std::endl;
// 
//    return myFit;
// }
// //______________________________________________________________________________
// TF1 *GetPolyFitToFXPR(TString fitName,int nOrder,int method,std::vector<int> fxprList,std::vector<gm2field::fixedProbeFrequency_t> fxprData,
//                       std::vector< std::vector<double> > &eps){
//    // fit the fixed probe data to an nth-order polynomial function 
//    // fit assumes the zeroth parameter is the time offset, first parameter is order  
//    // ROOT has trouble with large x values, so we subtract off a time constant defined as 
//    // the first event in the time series   
// 
//    const int N     = fxprData.size();
//    const int NFXPR = fxprList.size();
//    int startIndex  = fxprList[0];
//    int stopIndex   = fxprList[NFXPR-1];
// 
//    unsigned long long tStep = 1E+9; 
//    unsigned long long tMin  = fxprData[0].GpsTimeStamp[startIndex]; 
//    unsigned long long tMax  = fxprData[N-1].GpsTimeStamp[stopIndex]; 
// 
//    // get the average fixed probe plot 
//    TGraphErrors *gFPAVG = GetTGraphErrors(method,tMin,tMax,tStep,fxprList,fxprData);
//    gm2fieldUtil::Graph::SetGraphParameters(gFPAVG,20,kBlack);
// 
//    // time offset (to get fit to work) 
//    double t0 = tMin/1E+9;
// 
//    // custom fit -- first parameter must be the t0 offset  
//    const int NPAR = nOrder+2;
//    double par[NPAR];
//    for(int i=0;i<NPAR;i++){
//       if(i==0) par[i] = t0;
//       if(i==1) par[i] = nOrder;
//       if(i>1) par[i] = 1.;
//    }
// 
//    double TMIN = tMin/1E+9;
//    double TMAX = tMax/1E+9;
// 
//    // std::cout << Form("start time: %.0lf (%s)",TMIN,gm2fieldUtil::GetStringTimeStampFromUTC( (int)TMIN ).c_str() ) << std::endl;  
//    // std::cout << Form("end time:   %.0lf (%s)",TMAX,gm2fieldUtil::GetStringTimeStampFromUTC( (int)TMAX ).c_str() ) << std::endl;  
// 
//    double parLo[NPAR],parHi[NPAR]; 
//    for(int i=0;i<NPAR;i++){
//       if(i<2){
// 	 parLo[i] = 0;
// 	 parHi[i] = 0;
//       }else if(i==2){ 
// 	 parLo[i] = 40E+3;
// 	 parHi[i] = 70E+3;
//       }else{
// 	 parLo[i] = -10;
// 	 parHi[i] =  10;
//       }
//    }
// 
//    TF1 *myFit = new TF1(fitName,MyFitFunc_poly,TMIN,TMAX,NPAR);
//    myFit->SetLineColor(kRed+1);
//    myFit->FixParameter(0,par[0]);
//    myFit->FixParameter(1,par[1]);
//    for(int i=2;i<NPAR;i++) myFit->SetParameter(i,par[i]);
//    for(int i=2;i<NPAR;i++) myFit->SetParLimits(i,parLo[i],parHi[i]);
// 
//    // draw and fit
//    gROOT->SetBatch(kTRUE);  // turn off this canvas, we don't need it  
// 
//    TCanvas *c = new TCanvas("fitCanvas","Fit Canvas",500,500);
//    c->cd();
//    gFPAVG->Draw("ap");
//    TFitResultPtr fitResult = gFPAVG->Fit(fitName,"QS");
//    c->Update();
//    delete c;
//    gROOT->SetBatch(kFALSE);  // turn canvas making back on  
// 
//    std::cout << "Saving the covariance matrix..." << std::endl; 
//    TMatrixDSym myCov = fitResult->GetCovarianceMatrix();
//    // TMatrixD cor = fitResult->GetCorrelationMatrix();
//    // std::cout << "-- Covariance matrix: " << std::endl;
//    // myCov.Print(); 
// 
//    std::vector<double> COL;
//    for(int i=0;i<=nOrder+1;i++){    // row index
//       for(int j=0;j<=nOrder+1;j++){  // column index
// 	 COL.push_back(myCov[i][j]); 
// 	 // std::cout << myCov[i][j] << std::endl;
//       }
//       eps.push_back(COL);
//       COL.clear();
//    }
//    std::cout << "--> Done." << std::endl;
// 
//    return myFit;
// }
//______________________________________________________________________________
TGraphErrors *GetFitErrorBand(std::vector<double> x,TF1 *fit,std::vector< std::vector<double> > eps,double sf){
   // get fit parameters 
   const int NPAR = fit->GetNpar(); 
   double par[NPAR],parErr[NPAR];
   for(int i=0;i<NPAR;i++){
      par[i]    = fit->GetParameter(i); 
      parErr[i] = fit->GetParError(i); 
      // std::cout << par[i] << std::endl; 
   }
 
   // now get the fit error 
   double df=0;
   const int N = x.size();
   std::vector<double> ex,f,ef;
   for(int i=0;i<N;i++){
      df = sf*GetFitError(1,x[i],NPAR,par,parErr,eps);
      ef.push_back(df);
      f.push_back( fit->Eval(x[i]) );
      ex.push_back(0); 
   }

   TGraphErrors *g = gm2fieldUtil::Graph::GetTGraphErrors(x,f,ef);
   g->SetFillColor(kRed);
   g->SetFillStyle(3002);

   return g;
}
//______________________________________________________________________________
double MyFitFunc_pol2_simple(double *x,double *par){
   double X = x[0];
   double f = par[0] + par[1]*X + par[2]*X*X;
   return f;
}
//______________________________________________________________________________
double MyFitFunc_pol2(double *x,double *par){
   double X = x[0]-par[0];
   double f = par[1] + par[2]*X + par[3]*X*X;
   return f;
}
//______________________________________________________________________________
double MyFitFunc_pol3(double *x,double *par){
   double X = x[0]-par[0];
   double f = par[1] + par[2]*X + par[3]*X*X + par[4]*X*X*X;
   return f;
}
//______________________________________________________________________________
double MyFitFunc_pol4(double *x,double *par){
   double X = x[0]-par[0];
   const int nOrder = 4; 
   double f=0;
   for(int i=1;i<=nOrder+1;i++) f += par[i]*TMath::Power(X, (double)i-1 ); 
   return f;
}
//______________________________________________________________________________
double MyFitFunc_poly(double *x,double *par){
   // arbitrary polynomial function of nth order
   // par[0] = t0 
   // par[1] = polynomial order   
   double X = x[0]-par[0];
   const int nOrder = par[1]; 
   double f=0;
   for(int i=2;i<=nOrder+2;i++) f += par[i]*TMath::Power(X, (double)i-2 );  
   // std::cout << Form("x = %.3lf, f = %.3lf",X,f) << std::endl;  
   return f;
}
//______________________________________________________________________________
double GetFitError(int ifunc,double x,int npar,double *p,double *perr,vector< vector<double> > eps){
   const int NPAR = npar;
   double df_dp[NPAR];
   for(int i=0;i<npar;i++){
      df_dp[i] = FitFuncDerivative(ifunc,i,x,p);
   }  
   
   // diagonal elements 
   double T1=0;
   for(int i=0;i<npar;i++){
      T1 += TMath::Power(df_dp[i]*perr[i],2.)*eps[i][i];
   }  
   
   // off-diagonal elements 
   double T2=0;
   for(int i=0;i<npar;i++){
      for(int j=0;j<npar;j++){
         if(i!=j){
            T2 += df_dp[i]*df_dp[j]*eps[i][j]*perr[i]*perr[j];
         }  
      }  
   }  
   
   double df_sq = T1 + T2;
   double df    = TMath::Sqrt(df_sq);
   
   // check if we get a bad result (i.e., arg < 0), and return something huge if we do...   
   int IsNAN = TMath::IsNaN(df);
   if(IsNAN) df = 1E+6;
   
   return df;
}
//______________________________________________________________________________
double FitFuncDerivative(int type,int ipar,double x,double *p){
   double df=0;
   if(type==1) df = PolyFitFuncDerivative(ipar,x,p);
   return df;
}
//______________________________________________________________________________
double PolyFitFuncDerivative(int ipar,double x,double *p){
   // derivative of the polynomial fit function with respect to parameter p[ipar]. 
   double df=0;
   double X = x-p[0];   // remember p[0] is the time offset 
   if(ipar==0) df = 0;  // fixed parameter  
   if(ipar==1) df = 0;  // fixed parameter 
   if(ipar==2) df = 1.; 
   if(ipar==3) df = TMath::Power(X,1.);
   if(ipar==4) df = TMath::Power(X,2.);
   if(ipar==5) df = TMath::Power(X,3.);
   if(ipar==6) df = TMath::Power(X,4.);
   return df;
}
