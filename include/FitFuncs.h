#ifndef FIT_FUNCS_H
#define FIT_FUNCS_H

// fit functions 

#include <cstdlib> 
#include <iostream>
#include <vector>
#include <string> 

#include "TROOT.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TString.h"
#include "TFitResult.h"
#include "TMatrixDSym.h"

#include "RootTreeStructs.h"

#include "CustomGraph.h"

// TF1 *GetFitToFXPR(TString fitName,double (*fitFunc)(double *,double *),
//                   int npar,int method,std::vector<int> fxprList,std::vector<gm2field::fixedProbeFrequency_t> fxprData); 
// 
// TF1 *GetPolyFitToFXPR(TString fitName,int nOrder,int method,
//                       std::vector<int> fxprList,std::vector<gm2field::fixedProbeFrequency_t> fxprData,
//                       std::vector< std::vector<double> > &eps);
// 
// TF1 *GetPolyFitToFXPR(TString fitName,int nOrder,TGraph *gFPAVG,std::vector< std::vector<double> > &eps);

TGraphErrors *GetFitErrorBand(std::vector<double> x,TF1 *fit,std::vector< std::vector<double> > eps,double sf=1);

double MyFitFunc_pol2_simple(double *x,double *par); 
double MyFitFunc_pol2(double *x,double *par); 
double MyFitFunc_pol3(double *x,double *par); 
double MyFitFunc_pol4(double *x,double *par); 
double MyFitFunc_poly(double *x,double *par); 

double GetFitError(int ifunc,double x,int npar,double *p,double *perr,vector< vector<double> > eps);
double FitFuncDerivative(int type,int ipar,double x,double *p);
double PolyFitFuncDerivative(int ipar,double x,double *p); 

#endif 
