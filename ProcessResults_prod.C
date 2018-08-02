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

#include "./src/InputManager.C"
#include "./src/FitFuncs.C"
#include "./src/FXPRFuncs.C"
#include "./src/TRLYFuncs.C"
#include "./src/CalibFuncs.C"
#include "./src/Consolidate.C"
#include "./src/CustomUtilities.C"
#include "./src/CustomMath.C"
#include "./src/CustomGraph.C"
#include "./src/CustomImport.C"
#include "./src/CustomExport.C"
#include "./src/CustomAlgorithms.C"

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

   std::string anaDate = inputMgr->GetAnalysisDate();
   bool isBlind        = inputMgr->IsBlind();
   bool isFreeProton   = inputMgr->GetFreeProtonStatus(); 
   int probeNumber     = inputMgr->GetTrolleyProbe(); 

   date_t theDate;
   GetDate(theDate);

   char plotDir[200];
   sprintf(plotDir,"./plots/%s",theDate.getDateString().c_str());
   rc = MakeDirectory(plotDir);

   char outDir[200];
   sprintf(outDir,"./output"); 
   if(isBlind)  sprintf(outDir,"%s/blinded"  ,outDir);
   if(!isBlind) sprintf(outDir,"%s/unblinded",outDir);
   sprintf(outDir,"%s/%s",outDir,theDate.getDateString().c_str()); 
   rc = MakeDirectory(outDir);

   // results  
   char outPath_result[500]; 
   sprintf(outPath_result,"%s/results_%s.csv",outDir,anaDate.c_str());

   // load results
   result_prod_t result;
   rc = LoadResultsProdData(outPath_result,result); 

   double shot_err     = result.diffErr; 
   double shot_err_aba = result.diffErr_aba; 

   // misalignments 
   char outPath_misalign[500]; 
   sprintf(outPath_misalign,"%s/misalignment_results_%s.csv",outDir,anaDate.c_str());

   // load misalignment results 
   misalignment_t mErr; 
   rc = LoadMisalignmentData(outPath_misalign,mErr); 

   double tot_misalign_err     = TMath::Sqrt(mErr.dB_x*mErr.dB_x + mErr.dB_y*mErr.dB_y + mErr.dB_z*mErr.dB_z);
   double tot_misalign_err_aba = TMath::Sqrt(mErr.dB_x_aba*mErr.dB_x_aba + mErr.dB_y_aba*mErr.dB_y_aba + mErr.dB_z_aba*mErr.dB_z_aba);

   // load in perturbation data (just in case we need it) 
   char inpath_pert[200];
   perturbation_t ppPert;
   sprintf(inpath_pert,"./input/perturbation/pp-pert.csv");
   LoadPerturbationData(inpath_pert,ppPert);

   // compute errors from free proton corrections if necessary 
   double freeProtErr=0;
   if(isFreeProton){
      rc        = GetOmegaP_err(ppPert,freeProtErr);
      std::cout << "Computing uncertainty from free-proton corrections" << std::endl;
   } 

   char errStr[200],errStr_aba[200];
   if(isFreeProton){
      sprintf(errStr    ,"%.3lf +/- %.3lf +/- %.3lf",shot_err    ,tot_misalign_err,freeProtErr); 
      sprintf(errStr_aba,"%.3lf +/- %.3lf +/- %.3lf",shot_err_aba,tot_misalign_err_aba,freeProtErr); 
   }else{
      sprintf(errStr    ,"%.3lf +/- %.3lf",shot_err    ,tot_misalign_err); 
      sprintf(errStr_aba,"%.3lf +/- %.3lf",shot_err_aba,tot_misalign_err_aba); 
   } 

   result.mErr     = tot_misalign_err; 
   result.mErr_aba = tot_misalign_err_aba; 
   result.pErr     = freeProtErr; 
   result.pErr_aba = freeProtErr; 
  
   std::cout << Form("======================= PROBE %02d RESULTS =======================",probeNumber) << std::endl;
   std::cout << Form("[RAW]: %.3lf +/- %s",result.diff    ,errStr)     << std::endl;
   std::cout << Form("[ABA]: %.3lf +/- %s",result.diff_aba,errStr_aba) << std::endl;

   // results  
   char outPath_final[500]; 
   sprintf(outPath_final,"%s/results_final_%s.csv",outDir,anaDate.c_str());
   std::string outpath_final = outPath_final;  

   rc = PrintResults(outpath_final,result);

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
