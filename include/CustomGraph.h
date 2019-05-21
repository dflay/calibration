#ifndef CUSTOM_GRAPH_H
#define CUSTOM_GRAPH_H

// Custom plotting functions 

#include <cstdlib>
#include <vector>
#include <string>  

#include "TString.h"
#include "TGraph2D.h"

#include "gm2fieldGraph.h"

#include "grad_meas.h"
#include "plungingProbeAnaEvent.h"
#include "trolleyAnaEvent.h"
#include "fixedProbeEvent.h"
#include "sccEvent.h"

#include "TRLYFuncs.h"

// scc plots
TGraph *GetSCCPlot(int type,std::vector<surfaceCoilEvent_t> data); 

// gradient plots
TGraphErrors *GetTGraphErrors(std::vector<imposed_gradient_t> data); 

// plunging probe plots 
TGraph *GetPPTGraph1(TString xAxis,TString yAxis,std::vector<plungingProbeAnaEvent_t> data); 
TGraph *GetPPTGraph2(TString xAxis,TString yAxis,std::vector<plungingProbeAnaEvent_t> data); 
TGraph *GetPPTGraph3(TString xAxis,TString yAxis,plungingProbeAnaEvent_t data); 
TGraphErrors *GetPPScanGraph(TString xAxis,TString yAxis,std::vector<plungingProbeAnaEvent_t> data,double x0=0);
 
// TGraphErrors *GetTGraphErrors(TString xAxis,TString yAxis,std::vector<plungingProbeAnaEvent_t> data);
// TGraph2D *GetAzimuthalProjection(std::vector<plungingProbeAnaEvent_t> data,int units=gm2fieldUtil::Units::Hz); 
int FillPPVector1(TString axis,std::vector<plungingProbeAnaEvent_t> ppData,std::vector<double> &x);
int FillPPVector2(TString axis,std::vector<plungingProbeAnaEvent_t> ppData,std::vector<double> &x);
int FillPPVector3(TString axis,plungingProbeAnaEvent_t ppData,std::vector<double> &x);

// trolley plots
TGraph *GetTRLYPositionsGraph(); 
TGraph *GetTRLYTGraph(int probe,TString xAxis,TString yAxis,std::vector<trolleyAnaEvent_t> data,double sf=1);
TGraphErrors *GetSlicePlot(char axis,std::vector<trolleyAnaEvent_t> trlyData); 
TGraphErrors *GetSlicePlot(char axis,std::vector<trolleyAnaEvent_t> trlyData,
                           std::vector<double> &X,std::vector<double> &Y,std::vector<double> &EY); 
TGraph2D *GetAzimuthalProjection(std::vector<trolleyAnaEvent_t> data,int units=gm2fieldUtil::Constants::Hz); 
int FillTRVector(int probe,TString axis,std::vector<trolleyAnaEvent_t> data,double sf,std::vector<double> &x);

TGraphErrors *GetSCCTestGraphTRLY(int probe,TString xAxis,TString yAxis,std::vector< std::vector<sccTrlyEvent_t> > data);

TGraphErrors *GetTRLYTGraph_aziScan(int probe,double thr,TString yAxis,std::vector<trolleyAnaEvent_t> data);

// fxpr plots 
TGraph *GetFXPRTGraph(std::string xAxis,std::string yAxis,std::vector<fixedProbeEvent_t> data);
TGraphErrors *GetFXPRTGraph_avg(std::string xAxis,std::string yAxis,std::string yAxisErr,std::vector<averageFixedProbeEvent_t> data);

int FillFPVector_avg(std::string axis,std::vector<averageFixedProbeEvent_t> data,std::vector<double> &v);  
int FillFPVector(std::string axis,std::vector<fixedProbeEvent_t> data,std::vector<double> &v);  

// fixed probe plots 
// TGraph *GetTGraphNew(int method,unsigned long long timeStart,unsigned long long timeStop,unsigned long long timeStep,
//                   std::vector<int> fxprList,std::vector<gm2field::fixedProbeFrequency_t> fxprData);
// TGraph *GetTGraph(int method,unsigned long long timeStart,unsigned long long timeStop,unsigned long long timeStep,
//                   std::vector<int> fxprList,std::vector<gm2field::fixedProbeFrequency_t> fxprData);
// TGraph *GetTGraph(int method,unsigned long long t0,unsigned long long timeStart,unsigned long long timeStop,unsigned long long timeStep,
//                   std::vector<int> fxprList,std::vector<gm2field::fixedProbeFrequency_t> fxprData);
// 
// TGraph *GetInterpolatedTGraph(int method,std::vector<int> probe,
//                               unsigned long long tStart,unsigned long long tStop,unsigned long long tStep,
//                               std::vector<gm2field::fixedProbeFrequency_t> fxprData,
//                               std::vector<double> &stats);
// 
// TGraph *GetTGraph2Runs(int method,std::vector<int> probe, 
//                        unsigned long long tStart,unsigned long long tStop,unsigned long long tStep,
//                        std::vector<gm2field::fixedProbeFrequency_t> fxprData);
// 
// TGraphErrors *GetTGraphErrors(int method,unsigned long long timeStart,unsigned long long timeStop,unsigned long long timeStep,
//                               std::vector<int> fxprList,std::vector<gm2field::fixedProbeFrequency_t> fxprData);

#endif 
