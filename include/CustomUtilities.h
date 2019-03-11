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

#include "date.h" 

// void mkdirTree(string sub, string dir); 

int GetDate(date_t &aDate);
int MakeDirectory(const char *path);

std::string GetPath(std::string base,bool isBlind,std::string blindLabel,std::string date); 

#endif 
