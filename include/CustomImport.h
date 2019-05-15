#ifndef CUSTOM_IMPORT_H
#define CUSTOM_IMPORT_H

// extra import functions 

#include <cstdlib>
#include <iostream>  
#include <vector>
#include <string> 

#include "gm2fieldFunc.h"
#include "gm2fieldImport.h"
#include "gm2fieldRootHelper.h"

#include "sccEvent.h"
#include "runSummary.h"
#include "imageResult.h"
#include "imageParameter.h"
#include "trolleyAnaEvent.h"
#include "plungingProbeAnaEvent.h"
#include "fixedProbeEvent.h"
#include "results.h"
#include "misalignment.h"
#include "nmrAnaEvent.h"
#include "NMRDAQEvent.h"
#include "perturbation.h"
#include "nmr_meas.h"
#include "grad_meas.h"
#include "deltab.h"
#include "blind.h"
#include "Constants.h"

// Reading trolley data 
int GetTrolleyData(std::string date,int run,int method,std::vector<trolleyAnaEvent_t> &trlyEvent,std::string version);

// Reading PP data
int GetPlungingProbeData(int run,int prMethod,int ppMethod,std::vector<plungingProbeAnaEvent_t> &data,
                         std::string version,std::string nmrAnaVersion,std::string cutData="UNKNOWN.json",bool useNMRANA=true); 
int ModifyPlungingProbeData(int method,plungingProbeAnaEvent_t &data,std::string nmrAnaVersion,std::string cutFile); 

// Reading SCC data 
int GetSurfaceCoilData(int run,std::vector<surfaceCoilEvent_t> &data,std::string version);  

// Everything else 
int SetDataFileParameters(std::string version,std::string &fileName,std::string &dataPath); 
int LoadImageResults(std::string inpath,std::vector<imageResult_t> &data); 
int LoadIMGTimes(std::string type,int trial,std::vector<double> &time); 
int LoadImageParameters(std::string inpath,std::string type,std::vector<imageParameter_t> &data);

int LoadRunSummaryData(const char *inpath,runSummary_t &x); 
int LoadNMRDAQEventData(const char *inpath,std::vector<NMRDAQEvent_t> &event); 
int ImportNMRANAData(const char *inpath,std::vector<nmrAnaEvent_t> &Data,std::string cutFile);
int ImportDeltaBFileList_csv(const char *inpath,
                    std::vector<int> &x1,std::vector<std::string> &x2,
                    std::vector<double> &x3); 

int LoadImagesData(const char *inpath,int probe,double &image,double &image_err); 

int LoadTimes(int probe,int runPeriod,std::string type,std::string dev,std::vector<double> &time); 
int LoadSCCTimes(int probe,int runPeriod,std::string dev,std::vector<double> &sccOff,std::vector<double> &sccOn);

int LoadResultsProdData(const char *inpath,result_prod_t &data); 
int LoadResultsProdFinalData(const char *inpath,result_prod_t &data); 

int LoadMisalignmentData(const char *inpath,misalignment_t &data); 
int LoadCalibSwapData(const char *inpath,std::vector<calibSwap_t> &data);

int LoadImposedGradientData(const char *inpath,imposed_gradient_t &data); 
int LoadImposedGradientData(const char *inpath,std::vector<imposed_gradient_t> &data);
int LoadImposedAziGradData(const char *inpath,int probe,double &dBdz,double &dBdz_err); 
int LoadImposedAziGradData_bak(const char *inpath,int probe,double &dBdz); 

int LoadTrolleyDeltaBData(const char *inpath,trolleyDeltaB_t &data); 
int LoadTrolleyPositionData(const char *inpath,trolleyProbePosition_t &data);

int LoadPerturbationData(const char *inpath,int probe,perturbation_t &pert);
int LoadPerturbationData_json(const char *inpath,perturbation_t &pert);
int LoadFieldData(const char *inpath,nmr_meas_t &x); 

int LoadGradientData(const char *inpath,grad_meas_t &x);
int LoadGradientData(const char *inpath,std::vector<grad_meas_t> &x);

int LoadDeltaBData(const char *inpath,std::vector<deltab_t> &x);
int LoadDeltaBData_trlyXY(const char *inpath,int probe,std::vector<deltab_t> &x); 
int LoadDeltaBData_trlyXYZ(const char *inpath,int probe,std::vector<deltab_t> &x); 

int FindStartIndexTRLY(std::string date,int runNumber);
int LoadBlinding(blind_t *data); 

int SortRuns(std::vector<std::string> label,std::vector<int> allRuns,
             std::vector<int> &run,std::vector<int> &driftRun,std::vector<int> &index);
int SortRunsAlt(std::vector<std::string> label,std::vector<int> allRuns,
             std::vector<int> &run,std::vector<int> &index);

int ImportResults(std::string inpath,result_t &data); 

// templated functions
// trolley  
//______________________________________________________________________________
template <typename T> int GetTrolleyMultipoles(int run,std::vector<T> &data,std::string version){
   int rc=0;
   std::string dirName    = "TreeGenTrolley";
   std::string treeName   = "trolley";
   std::string branchName = "FieldMultipole";
   std::string fileName,dataPath;
   rc = SetDataFileParameters(version,fileName,dataPath);
   rc = gm2fieldUtil::RootHelper::GetDataFromTree<T>(run,dirName,treeName,branchName,data,-1,-1,fileName,dataPath);
   return rc;
}
//______________________________________________________________________________
template <typename T> int GetTrolleyFrequencies(int run,std::vector<T> &data,std::string version){
   int rc=0;
   std::string dirName    = "TreeGenTrolley";
   std::string treeName   = "trolley";
   std::string branchName = "ProbeFrequency";
   std::string fileName,dataPath;
   rc = SetDataFileParameters(version,fileName,dataPath);
   rc = gm2fieldUtil::RootHelper::GetDataFromTree<T>(run,dirName,treeName,branchName,data,-1,-1,fileName,dataPath);
   return rc;
}
//______________________________________________________________________________
template <typename T> int GetTrolleyTimeStamps(int run,std::vector<T> &data,std::string version){
   int rc=0;
   std::string dirName    = "TreeGenTrolley";
   std::string treeName   = "trolley";
   std::string branchName = "TimeStamp";
   std::string fileName,dataPath;
   rc = SetDataFileParameters(version,fileName,dataPath);
   rc = gm2fieldUtil::RootHelper::GetDataFromTree<T>(run,dirName,treeName,branchName,data,-1,-1,fileName,dataPath);
   return rc;
}
//______________________________________________________________________________
template <typename T> int GetTrolleyPosition(int run,std::vector<T> &data,std::string version){
   int rc=0;
   std::string dirName    = "TreeGenTrolley";
   std::string treeName   = "trolley";
   std::string branchName = "Position";
   std::string fileName,dataPath;
   rc = SetDataFileParameters(version,fileName,dataPath);
   rc = gm2fieldUtil::RootHelper::GetDataFromTree<T>(run,dirName,treeName,branchName,data,-1,-1,fileName,dataPath);
   return rc;
}
//______________________________________________________________________________
template <typename T> int GetTrolleyMonitor(int run,std::vector<T> &data,std::string version){
   int rc=0;
   std::string dirName    = "TreeGenTrolley";
   std::string treeName   = "trolley";
   std::string branchName = "Monitor";
   std::string fileName,dataPath;
   rc = SetDataFileParameters(version,fileName,dataPath);
   rc = gm2fieldUtil::RootHelper::GetDataFromTree<T>(run,dirName,treeName,branchName,data,-1,-1,fileName,dataPath);
   return rc;
}
//______________________________________________________________________________
template <typename T> int GetTrolleyGalil(int run,std::vector<T> &data,std::string version){
   int rc=0;
   std::string dirName    = "TreeGenGalilTrolley";
   std::string treeName   = "tGalil";
   std::string branchName = "Trolley";
   std::string fileName,dataPath;
   rc = SetDataFileParameters(version,fileName,dataPath);
   rc = gm2fieldUtil::RootHelper::GetDataFromTree<T>(run,dirName,treeName,branchName,data,-1,-1,fileName,dataPath);
   return rc;
}
// plunging probe
//______________________________________________________________________________
template <typename T> int GetPlungingProbeFrequencies(int run,std::vector<T> &data,std::string version){
   int rc=0;
   std::string dirName    = "TreeGenPlungingProbe";
   std::string treeName   = "plungingProbe";
   std::string branchName = "Frequency";
   std::string fileName,dataPath;
   rc = SetDataFileParameters(version,fileName,dataPath);
   rc = gm2fieldUtil::RootHelper::GetDataFromTree<T>(run,dirName,treeName,branchName,data,-1,-1,fileName,dataPath);
   return rc;
}
//______________________________________________________________________________
template <typename T> int GetPlungingProbeInfo(int run,std::vector<T> &data,std::string version){
   int rc=0;
   std::string dirName    = "TreeGenPlungingProbe";
   std::string treeName   = "plungingProbe";
   std::string branchName = "Info";
   std::string fileName,dataPath;
   rc = SetDataFileParameters(version,fileName,dataPath);
   rc = gm2fieldUtil::RootHelper::GetDataFromTree<T>(run,dirName,treeName,branchName,data,-1,-1,fileName,dataPath);
   return rc;
}
// fixed probes 
//______________________________________________________________________________
template <typename T> int GetFixedProbeFrequencies(int run,std::vector<T> &data,std::string version){
   int rc=0;
   std::string dirName    = "TreeGenFixedProbe";
   std::string treeName   = "fixedProbe";
   std::string branchName = "Frequency";
   std::string fileName,dataPath;
   rc = SetDataFileParameters(version,fileName,dataPath);
   rc = gm2fieldUtil::RootHelper::GetDataFromTree<T>(run,dirName,treeName,branchName,data,-1,-1,fileName,dataPath);
   return rc;
}
// surface coils
//______________________________________________________________________________
template <typename T> int GetSurfaceCoil(int run,std::vector<T> &data,std::string version){
   int rc=0;
   std::string dirName    = "TreeGenSurfaceCoil";
   std::string treeName   = "surfaceCoils";
   std::string branchName = "data";
   std::string fileName,dataPath;
   rc = SetDataFileParameters(version,fileName,dataPath);
   rc = gm2fieldUtil::RootHelper::GetDataFromTree<T>(run,dirName,treeName,branchName,data,-1,-1,fileName,dataPath);
   return rc;
}

#endif 
