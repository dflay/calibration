#ifndef CUSTOM_IMPORT_H
#define CUSTOM_IMPORT_H

// extra import functions 

#include <cstdlib> 
#include <vector>
#include <string> 

#include "results.h"
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

int LoadPerturbationData(const char *inpath,perturbation_t &pert);
int LoadFieldData(const char *inpath,nmr_meas_t &x); 
int LoadGradientData(const char *inpath,std::vector<grad_meas_t> &x);
int LoadDeltaBData(const char *inpath,std::vector<deltab_t> &x);

int FindStartIndexTRLY(std::string date,int runNumber);
int LoadBlinding(blind_t *data); 
int SortRuns(std::vector<std::string> label,std::vector<int> allRuns,
             std::vector<int> &run,std::vector<int> &driftRun,std::vector<int> &index);

int ImportResults(std::string inpath,result_t &data); 
 

#endif 
