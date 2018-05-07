#ifndef CUSTOM_EXPORT_H
#define CUSTOM_EXPORT_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <string> 

#include "TString.h"

int PrintToFile(const char *outpath,std::string label,const int N,double *x,double *x_err);
int PrintToFile(const char *outpath,std::string label,const int N,double *x);
int PrintToFile(const char *outpath,std::string label,double *x,double *x_err);
int PrintToFile(const char *outpath,double *x);

#endif 
