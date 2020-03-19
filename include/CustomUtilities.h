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
 
int PrintMessage(std::string label,std::string msg,int lineNumber=-1);
int GetDate(date_t &aDate);
int MakeDirectory(const char *path);

std::string GetPath(std::string base,bool isBlind,std::string blindLabel,std::string date,bool isSyst=false,int systDirNum=0); 

TLine **GetLines(int color,double min,double max,std::vector<double> x);

unsigned long long GetFXPRCutTime(std::string inpath,int probe,int index,bool &isMax);  

#endif 
