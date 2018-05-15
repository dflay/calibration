#ifndef FXPR_FUNCS_H
#define FXPR_FUNCS_H

// functions that use the fixed probe data structures  

#include <cstdlib> 
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>  
#include <cmath> 

#include "TStopwatch.h"
#include "TSpline.h"

#include "RootTreeStructs.h"
#include "gm2fieldMath.h"
#include "gm2fieldGraph.h"
#include "gm2fieldRootHelper.h"
#include "TemperatureSensor.h"

#include "Constants.h"
#include "fixedProbeEvent.h"
#include "plungingProbeAnaEvent.h"
#include "trolleyAnaEvent.h"

int FilterSingle(int nev,double T,std::vector<fixedProbeEvent_t> in,std::vector<fixedProbeEvent_t> &out);

int FindTransitionTimes(double thr,std::vector<fixedProbeEvent_t> fxprData,std::vector<double> &tLo,std::vector<double> &tHi);

int GetAverageFXPR(int method,unsigned long long time,std::vector<int> probe,
                   std::vector<gm2field::fixedProbeFrequency_t> fxprData,double &mean,double &stdev);

int GetAverageFXPR(int method,unsigned long long time,std::vector<int> probe,
                   std::vector<gm2field::fixedProbeFrequency_t> fxprData,
                   int NN,double *TIME,double *FREQ,double &mean,double &stdev);

int GetAverageFXPR(double time,std::vector<fixedProbeEvent_t> fxprData,double &mean,double &stdev); 

int GetAverageFXPRVectorsAlt(int method,
                             unsigned long long runEndPoint,unsigned long long runStartPoint,unsigned long long tStep,
                             std::vector<int> probeList,std::vector<gm2field::fixedProbeFrequency_t> fxprData,
                             std::vector<unsigned long long> &time,std::vector<double> &freq); 

int GetAverageFXPRVectors(int method,
                          unsigned long long t0,unsigned long long tStart,unsigned long long tStop,unsigned long long tStep,
                          std::vector<int> probeList,std::vector<gm2field::fixedProbeFrequency_t> fxprData,
                          std::vector<unsigned long long> &time,std::vector<double> &freq);  

int GetAverageFXPRVectorsNew(int method,
                             unsigned long long t0,unsigned long long tStart,unsigned long long tStop,unsigned long long tStep,
                             std::vector<int> probeList,std::vector<gm2field::fixedProbeFrequency_t> fxprData,
                             std::vector<unsigned long long> &time,std::vector<double> &freq);

int GetAverageFXPRVectorsNew(int method,
                             unsigned long long t0,unsigned long long tStart,unsigned long long tStop,unsigned long long tStep,
                             std::vector<int> probeList,std::vector<gm2field::fixedProbeFrequency_t> fxprData,
                             std::vector<fixedProbeEvent_t> &fxprDataAvg);

int FitFXPR(int lo,int hi,std::vector<fixedProbeEvent_t> fxprData,double *par);

double GetInterpolatedFXPRFreq(int method,int probe,unsigned long long time,std::vector<gm2field::fixedProbeFrequency_t> fxprData);
double GetInterpolatedFXPRFreq(int method,int probe,unsigned long long time,
                               std::vector<gm2field::fixedProbeFrequency_t> fxprData,
                               int NN,double *TIME,double *FREQ);  

#endif 
