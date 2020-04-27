#ifndef CUSTOM_EXPORT_H
#define CUSTOM_EXPORT_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <string> 

#include "TString.h"

int PrintToFile_TRLY_dBz_enc(const char *outpath,int probe,int run,double phi,double phi_err); 

int PrintTRLYFrequencies(const char *outpath,int probe,std::vector<trolleyAnaEvent_t> data); 
int PrintTRLYPositions(const char *outpath,std::vector<int> probe,
                       std::vector<double> r,std::vector<double> dr,
                       std::vector<double> y,std::vector<double> dy,
                       std::vector<double> z,std::vector<double> dz);

int PrintToFile_sccTimes(const char *outpath,std::vector<double> off,std::vector<double> on); 

int PrintToFile_1dbl(const char *outpath,std::vector<double> x);
int PrintToFile_2dbl(const char *outpath,std::vector<double> x1,std::vector<double> x2);
int PrintToFile_4dbl(const char *outpath,std::vector<double> x1,std::vector<double> x2,
                     std::vector<double> x3,std::vector<double> x4); 

int PrintToFile(const char *outpath,std::vector<std::string> label,std::vector<double> x); 
int PrintToFile(const char *outpath,std::vector<std::string> label,std::vector<double> x1,std::vector<double> x2);
int PrintToFile(const char *outpath,std::vector<std::string> label,std::vector<double> x1,std::vector<double> x2,std::vector<double> x3);

int PrintToFile(const char *outpath,std::vector<std::string> label,
                std::vector<double> x1,std::vector<double> x2,
                std::vector<double> x3,std::vector<double> x4); 
int PrintToFile(const char *outpath,std::vector<std::string> label,
                std::vector<double> x1,std::vector<double> x2,
                std::vector<double> x3,std::vector<double> x4, 
                std::vector<double> x5,std::vector<double> x6); 

int PrintToFile(const char *outpath,std::string label,const int N,double *x,double *x_err);
int PrintToFile(const char *outpath,std::string label,const int N,double *x);
int PrintToFile(const char *outpath,std::string label,double *x,double *x_err);
int PrintToFile(const char *outpath,double *x);

#endif 
