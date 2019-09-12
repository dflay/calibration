#ifndef CALIB_FUNCS_H
#define CALIB_FUNCS_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>

#include "nmr_meas.h"
#include "grad_meas.h"
#include "perturbation.h"

double GetDeltaTerm(double sigma,double delta_b,double delta_s,double delta_p,
                    double delta_rd,double delta_d);

int GetOmegaP_err(perturbation_t pert,double T,double &err); 
int GetOmegaP_free(nmr_meas_t pp,perturbation_t pert,double &freq_free,double &freq_free_err);
int GetOmegaP_free(perturbation_t pert,double freq,double freqErr,double temp,double tempErr,double &freqFree,double &freqFreeErr);
int GetDiamagneticShielding(double sigma,double dsigma,double T,double &SIG,double &ERR);
int GetMagneticSusceptibility(double chi,double T,double &CHI);
int GetBulkMagneticSusceptibility(double eps,double deps,double T,double &delta_b,double &delta_b_err); 

#endif
