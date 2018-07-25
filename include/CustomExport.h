#ifndef CUSTOM_EXPORT_H
#define CUSTOM_EXPORT_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <string> 

#include "TString.h"

int PrintTRLYPositions(const char *outpath,std::vector<int> probe,
                       std::vector<double> r,std::vector<double> dr,
                       std::vector<double> y,std::vector<double> dy,
                       std::vector<double> z,std::vector<double> dz);

int PrintToFile(const char *outpath,std::vector<std::string> label,
                std::vector<double> x1,std::vector<double> x2,
                std::vector<double> x3,std::vector<double> x4); 

int PrintToFile(const char *outpath,std::string label,const int N,double *x,double *x_err);
int PrintToFile(const char *outpath,std::string label,const int N,double *x);
int PrintToFile(const char *outpath,std::string label,double *x,double *x_err);
int PrintToFile(const char *outpath,double *x);

#endif 
