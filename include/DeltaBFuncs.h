#ifndef DELTAB_FUNCS_H
#define DELTAB_FUNCS_H

// functions to compute DeltaB for the PP or trolley 

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
#include "gm2fieldUnits.h"
#include "TemperatureSensor.h"

#include "Constants.h"
#include "plungingProbeAnaEvent.h"
#include "trolleyAnaEvent.h"
#include "FXPRFuncs.h"

int CalculatePPDeltaB_ABA(bool correctDrift,
                          plungingProbeAnaEvent_t bare,plungingProbeAnaEvent_t grad,plungingProbeAnaEvent_t bare2,
                          double &deltaB,double &deltaB_err);

int CalculatePPDeltaB(bool correctDrift,int method,
                      TGraph *gDrift,
                      plungingProbeAnaEvent_t bare,plungingProbeAnaEvent_t grad,
                      double &deltaB,double &deltaB_err,
                      double &drift ,double &drift_err);

int CalculatePPDeltaB(bool correctDrift,int method,std::vector<int> trlyList,
                      std::vector<trolleyAnaEvent_t> trlyData,
                      plungingProbeAnaEvent_t bare,plungingProbeAnaEvent_t grad,
                      double &deltaB,double &deltaB_err,
                      double &drift ,double &drift_err);

int CalculateTRLYDeltaB_Stationary_ABA(bool correctDrift,int probe,
                                       std::vector<trolleyAnaEvent_t> bare,
                                       std::vector<trolleyAnaEvent_t> grad,
                                       std::vector<trolleyAnaEvent_t> bare2,
                                       double &deltaB,double &deltaB_err);

int CalculateTRLYDeltaB_Stationary(bool correctDrift,int method,int probe,
                                   TGraph *gDrift,
                                   std::vector<trolleyAnaEvent_t> bare,std::vector<trolleyAnaEvent_t> grad,
                                   double &deltaB,double &deltaB_err,
                                   double &drift ,double &drift_err);

int CalculateTRLYDeltaB_Stationary(bool correctDrift,int method,int probe,std::vector<int> trlyList,
                                   std::vector<trolleyAnaEvent_t> trlyData,
                                   std::vector<trolleyAnaEvent_t> bare,std::vector<trolleyAnaEvent_t> grad,
                                   double &deltaB,double &deltaB_err,
                                   double &drift ,double &drift_err);

TGraph *CalculateTRLYDeltaB_Moving(int method,int probe,std::vector<int> fxprList,std::vector<gm2field::fixedProbeFrequency_t> fxprData,
                                   std::vector<trolleyAnaEvent_t> bare,std::vector<trolleyAnaEvent_t> grad); 

#endif 
