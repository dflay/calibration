#ifndef CONSOLIDATE_H
#define CONSOLIDATE_H

// Custom functions to fill vectors of custom structs    

#include <cstdlib> 
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>  
#include <cmath> 

#include "RootTreeStructs.h"
#include "gm2fieldMath.h"
#include "gm2fieldGraph.h"
#include "gm2fieldRootHelper.h"
#include "TemperatureSensor.h"

#include "Constants.h"
#include "nmrAnaEvent.h"
#include "plungingProbeAnaEvent.h"
#include "trolleyAnaEvent.h"
#include "CustomImport.h"

int ApplyBlindingPP(double blind,std::vector<plungingProbeAnaEvent_t> &x); 
int CopyPlungingProbe(plungingProbeAnaEvent_t x,plungingProbeAnaEvent_t &y);  
int CopyPlungingProbe(std::vector<plungingProbeAnaEvent_t> x,std::vector<plungingProbeAnaEvent_t> &y);  
int ConsolidatePPData(std::vector<int> fxprList,std::vector<gm2field::fixedProbeFrequency_t> fxprData,
                      std::vector<gm2field::plungingProbeFrequency_t> ppData,
                      std::vector<gm2field::plungingProbeInfo_t> ppInfo,
                      std::vector<plungingProbeAnaEvent_t> &ppEvent);
int GetPlungingProbeData(int runNumber,int freqMethod,std::vector<plungingProbeAnaEvent_t> &ppEvent);
int ModifyPlungingProbeData(int method,plungingProbeAnaEvent_t &Data);
int FilterPlungingProbeData(std::vector<int> subRun,
                            std::vector<plungingProbeAnaEvent_t> x,
                            std::vector<plungingProbeAnaEvent_t> &y); 

int ApplyBlindingTRLY(double blind,std::vector<trolleyAnaEvent_t> &x); 
int CopyTrolleyProbe(std::vector<trolleyAnaEvent_t> x,std::vector<trolleyAnaEvent_t> &y);  
int ConsolidateTrolleyData(int startIndex,int method,std::vector<gm2field::trolleyTimeStamp_t> trlyTime,
                           std::vector<gm2field::trolleyPosition_t> trlyPos,
                           std::vector<gm2field::trolleyMonitor_t> trlyMon,
                           std::vector<gm2field::trolleyProbeFrequency_t> trlyFreq, 
                           std::vector<trolleyAnaEvent_t> &trlyEvent);
int GetTrolleyData(std::string date,int runNumber,int freqMethod,std::vector<trolleyAnaEvent_t> &trlyEvent);

#endif 
