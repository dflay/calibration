#ifndef CUSTOM_ALGORITHMS_H
#define CUSTOM_ALGORITHMS_H

// useful analysis functions 
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>

#include "TF1.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSpline.h"

#include "RootTreeStructs.h"
#include "gm2fieldMath.h"
#include "gm2fieldGraph.h"
#include "gm2fieldRootHelper.h"
#include "gm2fieldUnits.h"
#include "TemperatureSensor.h"

#include "Constants.h"
#include "drift.h"
#include "plungingProbeAnaEvent.h"
#include "trolleyAnaEvent.h"
#include "fixedProbeEvent.h"
#include "FXPRFuncs.h"
#include "CustomGraph.h"

TGraph *RemoveTrend(TGraph *g1,TF1 *func); 
TGraph *GetDiffPlot(TGraph *g1,TGraph *g2);
TGraph *GetDiffPlot(TGraphErrors *g1,TGraphErrors *g2);
TGraph *GetDriftTGraph(int method,std::vector<int> driftRun,std::vector<int> fxprList,std::vector<double> &stats); 
TGraph *GetDriftTGraphR2R(int method,std::vector<int> driftRun,std::vector<int> fxprList,std::vector<double> &stats);

int CalculateTRLYAvg_Stationary(int probeNumber,std::vector<trolleyAnaEvent_t> Event,double &mean,double &stdev);

// PP functions 

int CalculateAveragePP(bool isDriftCor,std::vector<int> trlyList,std::vector<trolleyAnaEvent_t> trlyData,
                       std::vector<plungingProbeAnaEvent_t> ppData,double &B,double &B_err);

int CalculateAveragePP(bool isDriftCor,std::vector<int> fxprList,std::vector<gm2field::fixedProbeFrequency_t> fxprData,
                       std::vector<plungingProbeAnaEvent_t> ppData,double &B,double &B_err);  


// old approaches (depricated) 
int CorrectPPForDriftDuringMeasurement(int method,std::vector<int> fxprList,std::vector<gm2field::fixedProbeFrequency_t> fxprData,
                                       plungingProbeAnaEvent_t ppEvent,plungingProbeAnaEvent_t &ppEventCor);

int CorrectPPForDriftDuringMeasurement(int method,std::vector<int> fxprList,std::vector<gm2field::fixedProbeFrequency_t> fxprData,
                                       std::vector<plungingProbeAnaEvent_t> ppEvent,std::vector<plungingProbeAnaEvent_t> &ppEventCor);

// fit approaches 
int CorrectPPForDriftDuringMeasurement(int method,TF1 *fxprFit,
                                       plungingProbeAnaEvent_t ppEvent,plungingProbeAnaEvent_t &ppEventCor); 

int CorrectPPForDriftDuringMeasurement(int method,TF1 *fxprFit,
                                       std::vector<plungingProbeAnaEvent_t> ppEvent,std::vector<plungingProbeAnaEvent_t> &ppEventCor);
 
// new approaches
int CorrectPPForDriftDuringMeasurement(std::vector<fixedProbeEvent_t> fxprData,
                                       plungingProbeAnaEvent_t ppEvent,
                                       plungingProbeAnaEvent_t &ppEventCor,
                                       bool isScan=false);
int CorrectPPForDriftDuringMeasurement(std::vector<fixedProbeEvent_t> fxprData,
                                       std::vector<plungingProbeAnaEvent_t> ppEvent,
                                       std::vector<plungingProbeAnaEvent_t> &ppEventCor,
                                       bool isScan=false); 

// using the trolley 
int CorrectPPForDriftDuringMeasurementAlt(int method,std::vector<int> fxprList,std::vector<gm2field::fixedProbeFrequency_t> fxprData,
                                          std::vector<plungingProbeAnaEvent_t> ppEvent,std::vector<plungingProbeAnaEvent_t> &ppEventCor,
                                          std::vector<drift_t> &drift);  

// trolley functions
int CorrectTRLYForDriftDuringMeasurement(std::vector<fixedProbeEvent_t> fxprData,
                                         std::vector<trolleyAnaEvent_t> Event,std::vector<trolleyAnaEvent_t> &EventCor); 
 
int CorrectTRLYForDriftDuringMeasurement(int method,std::vector<int> fxprList,std::vector<gm2field::fixedProbeFrequency_t> fxprData,
                                         std::vector<trolleyAnaEvent_t> Event,std::vector<trolleyAnaEvent_t> &EventCor);

int CorrectTRLYForDriftDuringMeasurement(int method,TF1 *fxprFit,
                                         std::vector<trolleyAnaEvent_t> Event,std::vector<trolleyAnaEvent_t> &EventCor); 

int CorrectTRLYForDriftDuringMeasurement(int method,std::vector<int> trlyList,std::vector<trolleyAnaEvent_t> trlyData,
                                         std::vector<trolleyAnaEvent_t> Event,std::vector<trolleyAnaEvent_t> &EventCor); 

int GetAverageTRLY(unsigned long long time,std::vector<int> trlyList,std::vector<trolleyAnaEvent_t> trlyData,double &mean,double &stdev);

int FindTrolleyEvent(unsigned long long time, std::vector<trolleyAnaEvent_t> trly, int TP);

int CorrectTRLYForDriftAcrossRuns(unsigned long long tStart,unsigned long long tStop,
                                  TGraph *gDrift,std::vector<trolleyAnaEvent_t> &trlyData,
                                  double &drift,double &drift_err); 

#endif 
