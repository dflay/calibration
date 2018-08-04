#ifndef CUSTOM_IMPORT_H
#define CUSTOM_IMPORT_H

// extra import functions 

#include <cstdlib> 
#include <vector>
#include <string> 

#include "trolleyAnaEvent.h"
#include "results.h"
#include "misalignment.h"
#include "nmrAnaEvent.h"
#include "perturbation.h"
#include "nmr_meas.h"
#include "grad_meas.h"
#include "deltab.h"
#include "blind.h"

int ImportNMRANAData(const char *inpath,std::vector<nmrAnaEvent_t> &Data);
int ImportDeltaBFileList_csv(const char *inpath,
                    std::vector<int> &x1,std::vector<std::string> &x2,
                    std::vector<double> &x3); 

int LoadTRLYSwapTimes(int probe,std::vector<double> &time); 
int LoadResultsProdData(const char *inpath,result_prod_t &data); 
int LoadMisalignmentData(const char *inpath,misalignment_t &data); 
int LoadCalibSwapData(const char *inpath,std::vector<calibSwap_t> &data);
int LoadImposedGradientData(const char *inpath,imposed_gradient_t &data); 
int LoadImposedGradientData(const char *inpath,std::vector<imposed_gradient_t> &data);
int LoadImposedAziGradData(const char *inpath,int probe,double &dBdz); 
int LoadTrolleyDeltaBData(const char *inpath,trolleyDeltaB_t &data); 
int LoadTrolleyPositionData(const char *inpath,trolleyProbePosition_t &data);
int LoadPerturbationData(const char *inpath,perturbation_t &pert);
int LoadFieldData(const char *inpath,nmr_meas_t &x); 
int LoadGradientData(const char *inpath,grad_meas_t &x);
int LoadGradientData(const char *inpath,std::vector<grad_meas_t> &x);
int LoadDeltaBData(const char *inpath,std::vector<deltab_t> &x);
int LoadDeltaBData_trlyXY(const char *inpath,int probe,std::vector<deltab_t> &x); 

int FindStartIndexTRLY(std::string date,int runNumber);
int LoadBlinding(blind_t *data); 
int SortRuns(std::vector<std::string> label,std::vector<int> allRuns,
             std::vector<int> &run,std::vector<int> &driftRun,std::vector<int> &index);
int SortRunsAlt(std::vector<std::string> label,std::vector<int> allRuns,
             std::vector<int> &run,std::vector<int> &index);

int ImportResults(std::string inpath,result_t &data); 

#endif 
