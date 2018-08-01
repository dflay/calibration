#ifndef TRLY_FUNCS_H
#define TRLY_FUNCS_H 

#include <cstdlib> 
#include <vector>

#include "trolleyAnaEvent.h"

#include "gm2fieldMath.h"

// useful trolley functions

// int GetTrolleyProbePosition(int index,double *pos);
int GetTrolleyProbePositions(trolleyProbePosition_t &data);

int FindTransitionTimes(int step,double thr,std::vector<trolleyAnaEvent_t> Data,
                        std::vector< std::vector<double> > &timeLo,std::vector< std::vector<double> > &timeHi);
int FindTransitionTimes(int probe,int step,double thr,std::vector<trolleyAnaEvent_t> Data,
                        std::vector<double> &timeLo,std::vector<double> &timeHi);

int FilterSingle(int probe,int nev,double T,std::vector<trolleyAnaEvent_t> in,
                 std::vector<double> &freq,std::vector<double> &temp);
int FilterSingle(int probe,int nev,double T,std::vector<trolleyAnaEvent_t> in,
                 std::vector<double> &time,std::vector<double> &freq,std::vector<double> &temp);  
int FilterSingle(int probe,int nev,double T,std::vector<trolleyAnaEvent_t> in,std::vector<trolleyAnaEvent_t> &out); 

int GetTRLYStats_sccToggle(int probe,int nev,std::vector<double> time,std::vector<trolleyAnaEvent_t> Data,
                           std::vector<double> &TIME,std::vector<double> &MEAN,std::vector<double> &STDEV); 

int GetTRLYStatsAtTime(int probe,int nev,double fLO,std::vector<double> time,std::vector<trolleyAnaEvent_t> Data,
                       std::vector<double> &FREQ,std::vector<double> &FREQ_ERR,
                       std::vector<double> &TEMP,std::vector<double> &TEMP_ERR);

#endif 
