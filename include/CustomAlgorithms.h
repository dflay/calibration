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
#include "nmrAnaEvent.h"
#include "plungingProbeAnaEvent.h"
#include "trolleyAnaEvent.h"
// #include "fixedProbeEvent.h"
#include "CustomGraph.h"

TGraph *RemoveTrend(TGraph *g1,TF1 *func); 
TGraph *GetDiffPlot(TGraph *g1,TGraph *g2);
TGraph *GetDiffPlot(TGraphErrors *g1,TGraphErrors *g2);
// TGraph *GetDriftTGraph(int method,std::vector<int> driftRun,std::vector<int> fxprList,std::vector<double> &stats); 
// TGraph *GetDriftTGraphR2R(int method,std::vector<int> driftRun,std::vector<int> fxprList,std::vector<double> &stats);

int GetAvgSCCMagnitude(std::string type,std::vector<surfaceCoilEvent_t> data,double &mean,double &err,unsigned long long t0=0); 

int CheckShimGrad(int probe,std::vector<grad_meas_t> &data); 

int GetWeightedAverageStats(std::vector<double> x,std::vector<double> dx,double &mean,double &err,double &stdev); 

int FindGalilEvent(int probe,unsigned long long trlyTime,
                   std::vector<trolleyAnaEvent_t> trly,
                   std::vector<gm2field::galilTrolley_t> galil);
int SwapEntries(int i,std::vector<double> &x,std::vector<double> &y); 
int CheckDifference(std::vector<double> &x,std::vector<double> &dx);

int FindTRLYStopTimes(int probe,double angle,std::vector<trolleyAnaEvent_t> trlyData,
                      std::vector<gm2field::galilTrolley_t> trlyGalil,
                      std::vector<double> &time); 

int FindTransitionTimes(int type,int gradType,double thr,double delta,std::vector<surfaceCoilEvent_t> data,
                        std::vector<double> &timeOff,std::vector<double> &timeOn,int coilNum=-1); 

int CalculateTRLYAvg_Stationary(int probeNumber,std::vector<trolleyAnaEvent_t> Event,double &mean,double &stdev);

// copy functions
int CopyPlungingProbe(plungingProbeAnaEvent_t x,plungingProbeAnaEvent_t &y);
int CopyPlungingProbe(std::vector<plungingProbeAnaEvent_t> x,std::vector<plungingProbeAnaEvent_t> &y);
// int CopyPlungingProbe_nmrAna_to_ppAna_single(int method,nmrAnaEvent_t data,plungingProbeAnaEvent_t &y);
int CopyPlungingProbe_nmrAna_to_ppAna(int method,std::vector<nmrAnaEvent_t> data,plungingProbeAnaEvent_t &y);
int CopyTrolleyProbe(std::vector<trolleyAnaEvent_t> x,std::vector<trolleyAnaEvent_t> &y);
int FilterPlungingProbeData(std::vector<int> subRun,
                            std::vector<plungingProbeAnaEvent_t> x,
                            std::vector<plungingProbeAnaEvent_t> &y);

// Difference calculations 
int GetDifference(std::vector<double> A,std::vector<double> A_err,
                  std::vector<double> B,std::vector<double> B_err,
                  std::vector<double> &diff,std::vector<double> &diff_err);

int GetDifference_ABA(bool useTimeWeight,
                      std::vector<double> sccTime ,std::vector<double> scc ,std::vector<double> scc_err,
                      std::vector<double> bareTime,std::vector<double> bare,std::vector<double> bare_err,
                      std::vector<double> &diff_aba,std::vector<double> &diff_aba_err); 

int GetDifference_ABA_sccFirst(bool useTimeWeight,
                               std::vector<double> sccTime ,std::vector<double> scc ,std::vector<double> scc_err,
                               std::vector<double> bareTime,std::vector<double> bare,std::vector<double> bare_err,
                               std::vector<double> &diff_aba,std::vector<double> &diff_aba_err);

// same as "sccFirst", but better notation 
int GetDifference_ABA_final(bool useTimeWeight,
                            std::vector<double> A_time,std::vector<double> A,std::vector<double> A_err,
                            std::vector<double> B_time,std::vector<double> B,std::vector<double> B_err,
                            std::vector<double> &diff_aba,std::vector<double> &diff_aba_err);

// PP functions 
int CalculateAveragePP(std::vector<plungingProbeAnaEvent_t> ppData,double &B,double &B_err); 
int CalculateAveragePP(bool isDriftCor,std::vector<int> trlyList,std::vector<trolleyAnaEvent_t> trlyData,
                       std::vector<plungingProbeAnaEvent_t> ppData,double &B,double &B_err);

// int CalculateAveragePP(bool isDriftCor,std::vector<int> fxprList,std::vector<gm2field::fixedProbeFrequency_t> fxprData,
//                        std::vector<plungingProbeAnaEvent_t> ppData,double &B,double &B_err);  


// old approaches (depricated) 
// int CorrectPPForDriftDuringMeasurement(int method,std::vector<int> fxprList,std::vector<gm2field::fixedProbeFrequency_t> fxprData,
//                                        plungingProbeAnaEvent_t ppEvent,plungingProbeAnaEvent_t &ppEventCor);
// 
// int CorrectPPForDriftDuringMeasurement(int method,std::vector<int> fxprList,std::vector<gm2field::fixedProbeFrequency_t> fxprData,
//                                        std::vector<plungingProbeAnaEvent_t> ppEvent,std::vector<plungingProbeAnaEvent_t> &ppEventCor);

// fit approaches 
int CorrectPPForDriftDuringMeasurement(int method,TF1 *fxprFit,
                                       plungingProbeAnaEvent_t ppEvent,plungingProbeAnaEvent_t &ppEventCor); 

int CorrectPPForDriftDuringMeasurement(int method,TF1 *fxprFit,
                                       std::vector<plungingProbeAnaEvent_t> ppEvent,std::vector<plungingProbeAnaEvent_t> &ppEventCor);
 
// new approaches
// int CorrectPPForDriftDuringMeasurement(std::vector<fixedProbeEvent_t> fxprData,
//                                        plungingProbeAnaEvent_t ppEvent,
//                                        plungingProbeAnaEvent_t &ppEventCor,
//                                        bool isScan=false);
// int CorrectPPForDriftDuringMeasurement(std::vector<fixedProbeEvent_t> fxprData,
//                                        std::vector<plungingProbeAnaEvent_t> ppEvent,
//                                        std::vector<plungingProbeAnaEvent_t> &ppEventCor,
//                                        bool isScan=false); 

// using the trolley 
// int CorrectPPForDriftDuringMeasurementAlt(int method,std::vector<int> fxprList,std::vector<gm2field::fixedProbeFrequency_t> fxprData,
//                                           std::vector<plungingProbeAnaEvent_t> ppEvent,std::vector<plungingProbeAnaEvent_t> &ppEventCor,
//                                           std::vector<drift_t> &drift);  

// trolley functions
// int CorrectTRLYForDriftDuringMeasurement(std::vector<fixedProbeEvent_t> fxprData,
//                                          std::vector<trolleyAnaEvent_t> Event,std::vector<trolleyAnaEvent_t> &EventCor); 
 
// int CorrectTRLYForDriftDuringMeasurement(int method,std::vector<int> fxprList,std::vector<gm2field::fixedProbeFrequency_t> fxprData,
//                                          std::vector<trolleyAnaEvent_t> Event,std::vector<trolleyAnaEvent_t> &EventCor);

int CorrectTRLYForDriftDuringMeasurement(int method,TF1 *fxprFit,
                                         std::vector<trolleyAnaEvent_t> Event,std::vector<trolleyAnaEvent_t> &EventCor); 

int CorrectTRLYForDriftDuringMeasurement(int method,std::vector<int> trlyList,std::vector<trolleyAnaEvent_t> trlyData,
                                         std::vector<trolleyAnaEvent_t> Event,std::vector<trolleyAnaEvent_t> &EventCor); 

int GetAverageTRLY(unsigned long long time,std::vector<int> trlyList,std::vector<trolleyAnaEvent_t> trlyData,double &mean,double &stdev);

int FindTrolleyEvent(unsigned long long time, std::vector<trolleyAnaEvent_t> trly, int TP);
int FindTrolleyEvent(unsigned long long time, std::vector<averageTrolleyAnaEvent_t> trly); 

int CorrectTRLYForDriftAcrossRuns(unsigned long long tStart,unsigned long long tStop,
                                  TGraph *gDrift,std::vector<trolleyAnaEvent_t> &trlyData,
                                  double &drift,double &drift_err); 

#endif 
