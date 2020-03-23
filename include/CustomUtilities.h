#ifndef CUSTOM_UTILITIES_H
#define CUSTOM_UTILITIES_H

// a place to store useful miscellaneous functions
#include <cstdlib> 
#include <iostream>  
#include <ctime>
#include <sys/stat.h>
#include <dirent.h>
#include <string>
// #include <boost/filesystem.hpp>  

#include "TLine.h"

#include "date.h" 

// void mkdirTree(string sub, string dir);

enum cutType_util{
   kLowerBound = 0,
   kUpperBound = 1, 
   kRange      = 2
};
 
// int SplitString(std::string delim,std::string myStr,std::vector<std::string> &out); 
int PrintMessage(std::string label,std::string msg,int lineNumber=-1);
int GetDate(date_t &aDate);
int MakeDirectory(const char *path);
int GetFXPRCutTime(std::string inpath,int probe,int index,
                   unsigned long long &tMin,unsigned long long &tMax);  

std::string GetPath(std::string base,bool isBlind,std::string blindLabel,std::string date,bool isSyst=false,int systDirNum=0); 

TLine **GetLines(int color,double min,double max,std::vector<double> x);


#endif 
