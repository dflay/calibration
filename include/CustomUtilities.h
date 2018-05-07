#ifndef CUSTOM_UTILITIES_H
#define CUSTOM_UTILITIES_H

// a place to store useful miscellaneous functions
#include <cstdlib> 
#include <iostream>  
#include <ctime>
#include <sys/stat.h>
#include <dirent.h>

#include "date.h" 

int GetDate(date_t &aDate);
int MakeDirectory(const char *path); 

#endif 
