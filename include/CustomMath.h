#ifndef CUSTOM_MATH_H
#define CUSTOM_MATH_H

// Test out grabbing PP data from the ROOT Tree,
// getting average field from fixed probes at a given point in time   

#include <cstdlib> 
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>  
#include <cmath> 

#include "TMath.h"
#include "TString.h"
#include "TSpline.h"

#include "RootTreeStructs.h"
#include "gm2fieldMath.h"
#include "gm2fieldGraph.h"
#include "gm2fieldRootHelper.h"
#include "TemperatureSensor.h"

#include "Constants.h"
#include "plungingProbeAnaEvent.h"
#include "trolleyAnaEvent.h"

int GetStats(std::string varName,int probe,std::vector<trolleyAnaEvent_t> data,double &mean,double &stdev);
int GetStats(std::string varName,std::vector<calibSwap_t> data,double &mean,double &stdev);  

int GetStats(std::vector<double> x,double &mean,double &stdev); 
int GetStats(int probeNumber,int method,std::vector<gm2field::fixedProbeFrequency_t> Data,double &mean,double &stdev); 
int GetStats(int probeNumber,int method,std::vector<gm2field::trolleyProbeFrequency_t> Data,double &mean,double &stdev); 
int GetStats(TString axis,std::vector<gm2field::psFeedback_t> Data,double &mean,double &stdev); 

double MultipoleFitFunc(double *x,double *par); 

#endif 
