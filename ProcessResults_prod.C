// Compute final uncertainties 
// - Misalignment dx, dy, dz in quadrature 
// - Free-proton corrections   

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
#include "TText.h"

#include "RootTreeStructs.h"
#include "gm2fieldMath.h"
#include "gm2fieldGraph.h"
#include "gm2fieldRootHelper.h"
#include "gm2fieldUnits.h"
#include "TemperatureSensor.h"
#include "MovingAverage.h"

#include "./include/date.h"
#include "./include/Constants.h"
#include "./include/fixedProbeEvent.h"
#include "./include/trolleyAnaEvent.h"
#include "./include/nmr_meas.h"
#include "./include/perturbation.h" 

#include "./src/CustomUtilities.C"
#include "./src/CustomMath.C"
#include "./src/CustomImport.C"
#include "./src/CustomExport.C"
#include "./src/CustomAlgorithms.C"
#include "./src/CustomGraph.C"
#include "./src/OscFuncs.C"
#include "./src/InputManager.C"
#include "./src/FitFuncs.C"
#include "./src/TRLYFuncs.C"
#include "./src/CalibFuncs.C"

double gMarkerSize = 0.8;
double gXSize      = 0.05;  
double gYSize      = 0.06;  

int PrintResults(std::string outpath,result_prod_t result); 

int ProcessResults_prod(std::string configFile){

   int rc=0;
   int method = gm2fieldUtil::Constants::kPhaseDerivative;

   InputManager *inputMgr = new InputManager();
   inputMgr->Load(configFile);
   inputMgr->Print();

   std::string anaDate    = inputMgr->GetAnalysisDate();
   std::string blindLabel = inputMgr->GetBlindLabel();
   bool isBlind           = inputMgr->IsBlind();
   int probeNumber        = inputMgr->GetTrolleyProbe(); 
   int runPeriod          = inputMgr->GetRunPeriod(); 

   date_t theDate;
   GetDate(theDate);
  
   std::string outDir  = GetPath("output",isBlind,blindLabel,theDate.getDateString());
 
   // results  
   char outPath_result[500]; 
   sprintf(outPath_result,"%s/results_pr-%02d.csv",outDir.c_str(),probeNumber);

   result_prod_t result;
   rc = LoadResultsProdData(outPath_result,result); 

   char outPath_result_free[500]; 
   sprintf(outPath_result_free,"%s/results_free-prot_pr-%02d.csv",outDir.c_str(),probeNumber);

   result_prod_t result_free;
   rc = LoadResultsProdData(outPath_result_free,result_free); 

   double shot_err     = result.diffErr; 
   double shot_err_aba = result.diffErr_aba; 

   double shot_err_free     = result_free.diffErr; 
   double shot_err_free_aba = result_free.diffErr_aba; 

   // misalignments 
   char outPath_misalign[500]; 
   sprintf(outPath_misalign,"%s/misalignment_results_pr-%02d.csv",outDir.c_str(),probeNumber);

   // load misalignment results 
   misalignment_t mErr; 
   rc = LoadMisalignmentData(outPath_misalign,mErr); 

   double tot_misalign_err     = TMath::Sqrt(mErr.dB_x*mErr.dB_x + mErr.dB_y*mErr.dB_y + mErr.dB_z*mErr.dB_z);
   double tot_misalign_err_aba = TMath::Sqrt(mErr.dB_x_aba*mErr.dB_x_aba + mErr.dB_y_aba*mErr.dB_y_aba + mErr.dB_z_aba*mErr.dB_z_aba);

   // load in perturbation data (just in case we need it) 
   char inpath_pert[200];
   perturbation_t ppPert;
   sprintf(inpath_pert,"./input/perturbation/pp-pert_run-%d.json",runPeriod);
   LoadPerturbationData_json(inpath_pert,ppPert);

   // compute errors from free proton corrections 
   double freeProtErr=0;
   rc = GetOmegaP_err(ppPert,freeProtErr);

   char errStr[200],errStr_aba[200];
   sprintf(errStr    ,"%.3lf +/- %.3lf",shot_err    ,tot_misalign_err); 
   sprintf(errStr_aba,"%.3lf +/- %.3lf",shot_err_aba,tot_misalign_err_aba); 

   char errStr_free[200],errStr_free_aba[200];
   sprintf(errStr_free    ,"%.3lf +/- %.3lf +/- %.3lf",shot_err_free    ,tot_misalign_err,freeProtErr); 
   sprintf(errStr_free_aba,"%.3lf +/- %.3lf +/- %.3lf",shot_err_free_aba,tot_misalign_err_aba,freeProtErr); 

   result.mErr     = tot_misalign_err; 
   result.mErr_aba = tot_misalign_err_aba; 
   result.pErr     = 0.; 
   result.pErr_aba = 0.; 

   result_free.mErr     = tot_misalign_err; 
   result_free.mErr_aba = tot_misalign_err_aba; 
   result_free.pErr     = freeProtErr; 
   result_free.pErr_aba = freeProtErr; 
  
   std::cout << Form("======================= PROBE %02d RESULTS =======================",probeNumber) << std::endl;
   std::cout << "Bare" << std::endl;
   std::cout << Form("[RAW]: %.3lf +/- %s",result.diff    ,errStr)     << std::endl;
   std::cout << Form("[ABA]: %.3lf +/- %s",result.diff_aba,errStr_aba) << std::endl;
   std::cout << "Free proton" << std::endl;
   std::cout << Form("[RAW]: %.3lf +/- %s",result_free.diff    ,errStr_free)     << std::endl;
   std::cout << Form("[ABA]: %.3lf +/- %s",result_free.diff_aba,errStr_free_aba) << std::endl;

   // results  
   char outPath_final[500]; 
   sprintf(outPath_final,"%s/results_final_pr-%02d.csv",outDir.c_str(),probeNumber);
   std::string outpath_final = outPath_final;  
   rc = PrintResults(outpath_final,result);

   char outPath_final_free[500]; 
   sprintf(outPath_final_free,"%s/results_final_free-prot_pr-%02d.csv",outDir.c_str(),probeNumber);
   std::string outpath_final_free = outPath_final_free;  
   rc = PrintResults(outpath_final_free,result_free);

   return 0;
}
//______________________________________________________________________________
int PrintResults(std::string outpath,result_prod_t result){

   char header[500]; 
   // prepare the header 
   sprintf(header,"#type,diff,shot-err,misalign-err,free-prot-err"); 
   char myStr[1000]; 
   sprintf(myStr    ,"raw,%.3lf,%.3lf,%.3lf,%.3lf",result.diff,result.diffErr,result.mErr,result.pErr); 
   char myStr_aba[1000]; 
   sprintf(myStr_aba,"ABA,%.3lf,%.3lf,%.3lf,%.3lf",result.diff_aba,result.diffErr_aba,result.mErr_aba,result.pErr_aba); 

   std::ofstream outfile;
   outfile.open(outpath.c_str());
   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      outfile << header    << std::endl; 
      outfile << myStr     << std::endl;
      outfile << myStr_aba << std::endl;
      std::cout << "The data has been written to file: " << outpath << std::endl;
      outfile.close(); 
   }  
   return 0;
}
