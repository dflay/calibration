#ifndef TRLY_FUNCS_H
#define TRLY_FUNCS_H 

#include <cstdlib> 
#include <vector>
#include <string> 

#include "gm2fieldMath.h"

#include "nmr_meas.h"
#include "trolleyAnaEvent.h"
#include "fixedProbeEvent.h"

// useful trolley functions

// int GetTRLYStatsAtTime(bool UseTempCor,int probe,int nev,double fLO,std::vector<double> time,
//                        std::vector<trolleyAnaEvent_t> Data,std::vector<trolleySwapEvent_t> &Event); 

// int GetTRLYStats_sccToggle(int probe,int nev,std::vector<double> time,std::vector<trolleyAnaEvent_t> Data,
//                            std::vector<double> &TIME,std::vector<double> &MEAN,std::vector<double> &STDEV);

int GetTRLYStatsAtTime(bool UseTempCor,bool UseOscCor,int probe,int nev,double fLO,
                       std::vector<double> time,std::vector<averageFixedProbeEvent_t> fxpr,
                       std::vector<trolleyAnaEvent_t> Data,std::vector<calibSwap_t> &Event,double t0=0,double dsigdT=0);

int GetTRLYStatsAtTime_hybrid(bool UseTempCor,bool UseOscCor,int probe,int nev,double fLO,
                              std::vector<double> time,std::vector<averageFixedProbeEvent_t> fxpr,
                              std::vector<trolleyAnaEvent_t> Data,std::vector<calibSwap_t> &Event,
                              std::vector<double> t0,double dsigdT=0); 

int FindTransitionTimes(int step,double thr,std::vector<trolleyAnaEvent_t> Data,
                        std::vector< std::vector<double> > &timeLo,std::vector< std::vector<double> > &timeHi);
int FindTransitionTimes(int probe,int step,double thr,std::vector<trolleyAnaEvent_t> Data,
                        std::vector<double> &timeLo,std::vector<double> &timeHi);

int FilterSingle(std::string var,int probe,int nev,double T,std::vector<trolleyAnaEvent_t> in,std::vector<double> &x); 


int GetTRLYStats_sccToggle_old(bool useOscCor,int probe,int nev,std::vector<double> time,
                               std::vector<averageFixedProbeEvent_t> fxpr,std::vector<trolleyAnaEvent_t> Data,
                               std::vector<double> &TIME,std::vector<double> &MEAN,std::vector<double> &STDEV);  

int GetTRLYStatsAtTime_old(int probe,int nev,double fLO,std::vector<double> time,std::vector<trolleyAnaEvent_t> Data,
                           std::vector<double> &FREQ,std::vector<double> &FREQ_ERR,
                           std::vector<double> &TEMP,std::vector<double> &TEMP_ERR);

int GetTRLYStats_sccToggle(bool useOscCor,int probe,int nev,std::vector<double> time,
                           std::vector<averageFixedProbeEvent_t> fxpr,std::vector<trolleyAnaEvent_t> Data,
                           std::vector<double> &TIME,std::vector<double> &MEAN,std::vector<double> &STDEV,
                           double &mean_phi,double &stdev_phi,double t0=0);  

int GetTRLYStatsAtTime(int probe,int nev,double fLO,std::vector<double> time,std::vector<trolleyAnaEvent_t> Data,
                       std::vector<double> &FREQ,std::vector<double> &FREQ_ERR,
                       std::vector<double> &TEMP,std::vector<double> &TEMP_ERR);

#endif 
