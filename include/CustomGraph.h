#ifndef CUSTOM_GRAPH_H
#define CUSTOM_GRAPH_H

// Custom plotting functions 

#include <cstdlib>
#include <vector>
#include <string>  

#include "TString.h"
#include "TGraph2D.h"

#include "gm2fieldGraph.h"

#include "plungingProbeAnaEvent.h"
#include "trolleyAnaEvent.h"

#include "FXPRFuncs.h"

// plunging probe plots 
TGraph *GetPPTGraph1(TString xAxis,TString yAxis,std::vector<plungingProbeAnaEvent_t> data); 
TGraph *GetPPTGraph2(TString xAxis,TString yAxis,std::vector<plungingProbeAnaEvent_t> data); 
TGraph *GetPPTGraph3(TString xAxis,TString yAxis,plungingProbeAnaEvent_t data); 
// TGraphErrors *GetTGraphErrors(TString xAxis,TString yAxis,std::vector<plungingProbeAnaEvent_t> data);
// TGraph2D *GetAzimuthalProjection(std::vector<plungingProbeAnaEvent_t> data,int units=gm2fieldUtil::Units::Hz); 
int FillPPVector1(TString axis,std::vector<plungingProbeAnaEvent_t> ppData,std::vector<double> &x);
int FillPPVector2(TString axis,std::vector<plungingProbeAnaEvent_t> ppData,std::vector<double> &x);
int FillPPVector3(TString axis,plungingProbeAnaEvent_t ppData,std::vector<double> &x);

// trolley plots
TGraph *GetTRLYTGraph(int probe,TString xAxis,TString yAxis,std::vector<trolleyAnaEvent_t> data);
TGraphErrors *GetSlicePlot(char axis,std::vector<trolleyAnaEvent_t> trlyData); 
TGraphErrors *GetSlicePlot(char axis,std::vector<trolleyAnaEvent_t> trlyData,
                           std::vector<double> &X,std::vector<double> &Y,std::vector<double> &EY); 
TGraph2D *GetAzimuthalProjection(std::vector<trolleyAnaEvent_t> data,int units=gm2fieldUtil::Constants::Hz); 
int FillTRVector(int probe,TString axis,std::vector<trolleyAnaEvent_t> data,std::vector<double> &x);

// fixed probe plots 
TGraph *GetTGraphNew(int method,unsigned long long timeStart,unsigned long long timeStop,unsigned long long timeStep,
                  std::vector<int> fxprList,std::vector<gm2field::fixedProbeFrequency_t> fxprData);
TGraph *GetTGraph(int method,unsigned long long timeStart,unsigned long long timeStop,unsigned long long timeStep,
                  std::vector<int> fxprList,std::vector<gm2field::fixedProbeFrequency_t> fxprData);
TGraph *GetTGraph(int method,unsigned long long t0,unsigned long long timeStart,unsigned long long timeStop,unsigned long long timeStep,
                  std::vector<int> fxprList,std::vector<gm2field::fixedProbeFrequency_t> fxprData);

TGraph *GetInterpolatedTGraph(int method,std::vector<int> probe,
                              unsigned long long tStart,unsigned long long tStop,unsigned long long tStep,
                              std::vector<gm2field::fixedProbeFrequency_t> fxprData,
                              std::vector<double> &stats);

TGraph *GetTGraph2Runs(int method,std::vector<int> probe, 
                       unsigned long long tStart,unsigned long long tStop,unsigned long long tStep,
                       std::vector<gm2field::fixedProbeFrequency_t> fxprData);

TGraphErrors *GetTGraphErrors(int method,unsigned long long timeStart,unsigned long long timeStop,unsigned long long timeStep,
                              std::vector<int> fxprList,std::vector<gm2field::fixedProbeFrequency_t> fxprData);

#endif 
