// Compare results 

#include <cstdlib> 
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>  
#include <cmath> 

#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSpline.h"

#include "RootTreeStructs.h"
#include "gm2fieldMath.h"
#include "gm2fieldGraph.h"
#include "gm2fieldRootHelper.h"
#include "gm2fieldUnits.h"
#include "TemperatureSensor.h"

#include "./include/Constants.h"
#include "./include/results.h"
#include "./include/perturbation.h"
#include "./include/deltab.h"
#include "./include/plungingProbeAnaEvent.h"
#include "./include/trolleyAnaEvent.h"
#include "./include/nmr_meas.h"
#include "./include/grad_meas.h"
#include "./include/date.h"

#include "./src/FXPRFuncs.C"
#include "./src/Consolidate.C"
#include "./src/CustomMath.C"
#include "./src/CustomGraph.C"
#include "./src/CustomImport.C"
#include "./src/CustomExport.C"
#include "./src/CustomAlgorithms.C"
#include "./src/DeltaBFuncs.C"
#include "./src/CalibFuncs.C"
#include "./src/ErrorFuncs.C"
#include "./src/CustomUtilities.C"

int CompareResults(){

   int rc=0;

   // date_t theDate; 
   // rc = GetDate(theDate); 

   result_t res1,res2; 
   
   std::string prefix = "./output/blinded/04-30-18";
   char path[200]; 
   // run 3077 
   sprintf(path,"%s/trly-shim-3077/fxpr-set-2/results.csv",prefix.c_str()); 
   std::string inpath = path; 
   ImportResults(inpath,res1);  

   // run 3084 
   sprintf(path,"%s/trly-shim-3084/fxpr-set-2/results.csv",prefix.c_str()); 
   inpath = path; 
   ImportResults(inpath,res2);  
   
   // average over these results 
   std::vector<double> diff,diff_err,trly,trly_err;
 
   diff.push_back(res1.diff_fxpr); 
   diff.push_back(res2.diff_fxpr); 
   diff_err.push_back(res1.diff_fxpr_err); 
   diff_err.push_back(res2.diff_fxpr_err); 

   double arg = res1.trly_fxpr + res1.driftShim_fxpr; 
   trly.push_back(arg); 
   arg = res2.trly_fxpr + res2.driftShim_fxpr; 
   trly.push_back(arg);
 
   trly_err.push_back(res1.trly_fxpr_err); 
   trly_err.push_back(res2.trly_fxpr_err); 

   double diff_stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(diff);                 // add in standard deviation 
   double diffAvg    = gm2fieldUtil::Math::GetMean<double>(diff); 
   double diffErr    = gm2fieldUtil::Math::GetMean<double>(diff_err);  // average over the error for each run  
   diffErr           = TMath::Sqrt( diffErr*diffErr + diff_stdev*diff_stdev); 

   double trly_stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(trly);                 // add in standard deviation 
   double trlyAvg    = gm2fieldUtil::Math::GetMean<double>(trly); 
   double trlyErr    = gm2fieldUtil::Math::GetMean<double>(trly_err);  // average over the error for each run  
   trlyErr           = TMath::Sqrt( trlyErr*trlyErr + trly_stdev*trly_stdev);  

   std::cout << "PP (uncorrected): " << std::endl;
   std::cout << Form("raw  = %.3lf +/- %.3lf Hz (%.3lf ppb)",res1.ppRaw,res1.ppRaw_err,res1.ppRaw_err/0.06179) << std::endl; 
   std::cout << Form("fxpr = %.3lf +/- %.3lf Hz (%.3lf ppb)",res1.ppRaw_fxpr,res1.ppRaw_fxpr_err,res1.ppRaw_fxpr_err/0.06179) << std::endl; 
   std::cout << "PP (free): " << std::endl;
   std::cout << Form("raw  = %.3lf +/- %.3lf Hz (%.3lf ppb)",res1.ppFree,res1.ppFree_err,res1.ppFree_err/0.06179) << std::endl; 
   std::cout << Form("fxpr = %.3lf +/- %.3lf Hz (%.3lf ppb)",res1.ppFree_fxpr,res1.ppFree_fxpr_err,res1.ppFree_fxpr_err/0.06179) << std::endl; 


   std::cout << "USING TRLY RUN 3077" << std::endl;
   std::cout << Form("raw   = %.3lf +/- %.3lf Hz (%.3lf ppb)",res1.trly     ,res1.trly_err     ,res1.trly_err/0.06179) << std::endl; 
   std::cout << Form("fxpr  = %.3lf +/- %.3lf Hz (%.3lf ppb)",trly[0],trly_err[0],trly_err[0]/0.06179) << std::endl;
   std::cout << Form("PP-TR = %.3lf +/- %.3lf Hz (%.3lf ppb)",res1.diff_fxpr,res1.diff_fxpr_err,res1.diff_fxpr_err/0.06179) << std::endl;

   std::cout << "USING TRLY RUN 3084" << std::endl;
   std::cout << Form("raw   = %.3lf +/- %.3lf Hz (%.3lf ppb)",res2.trly     ,res2.trly_err     ,res2.trly_err/0.06179) << std::endl; 
   std::cout << Form("fxpr  = %.3lf +/- %.3lf Hz (%.3lf ppb)",trly[1],trly_err[1],trly_err[1]/0.06179) << std::endl;
   std::cout << Form("PP-TR = %.3lf +/- %.3lf Hz (%.3lf ppb)",res2.diff_fxpr,res2.diff_fxpr_err,res2.diff_fxpr_err/0.06179) << std::endl;

   std::cout << "COMBINED RESULTS" << std::endl;
   std::cout << Form("trly  = %.3lf +/- %.3lf (%.3lf)",trlyAvg,trlyErr,trlyErr/0.06179) << std::endl;
   std::cout << Form("PP-TR = %.3lf +/- %.3lf (%.3lf)",diffAvg,diffErr,diffErr/0.06179) << std::endl;

   return 0;
}

