#ifndef ERROR_FUNCS_H
#define ERROR_FUNCS_H

#include <cstdlib>
#include <vector>
#include <string>

#include "TMath.h"

#include "deltab.h"
#include "grad_meas.h"

int GetDriftError(std::vector<deltab_t> pp,std::vector<deltab_t> trly,
                  std::vector<grad_meas_t> gradImposed,std::vector<grad_meas_t> gradShim,
                  std::vector<double> &errFXPR,std::vector<double> &errTRLY);

double GetDeltaDQ(std::string type,deltab_t pp,deltab_t trly,grad_meas_t gradImposed);  

#endif 
