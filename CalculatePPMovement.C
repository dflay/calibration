// determine the distance to move the PP to align with 
// the trolley probe based upon Delta-B measurements.  

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

int CalculatePPMovement(std::string date,int probeNumber){

   int rc=0;
   int method = gm2fieldUtil::Constants::kPhaseDerivative;

   date_t theDate; 
   rc = GetDate(theDate); 

   char outdir[200],outpath[200],prefix[200];
   sprintf(outdir,"./output/delta-b/%02d-%02d-%02d",theDate.month,theDate.day,theDate.year-2000);  
   sprintf(prefix,"%s",outdir); 

   std::cout << "----------------------------------" << std::endl;
   std::cout << "CALCULATE DISTANCE TO MOVE PP TO ALIGN WITH TROLLEY" << std::endl;

   const int NA = 2;
   std::string gradName[NA] = {"rad","vert"; // ,"azi"};
 
   char inpath[200]; 
 
   // bring in DeltaB data
   std::vector<deltab_t> pp,trly;
   for(int i=0;i<NA;i++){
      sprintf(inpath,"%s/dB-pp_%s-grad_%s.csv",prefix,gradName[i].c_str(),date.c_str());
      LoadDeltaBData(inpath,pp); 
      sprintf(inpath,"%s/dB-trly_%s-grad_pr-%02d_%s.csv",prefix,gradName[i].c_str(),probeNumber,date.c_str());
      LoadDeltaBData(inpath,trly); 
   }

   // load in imposed gradients
   std::vector<grad_meas_t> gradient; 
   for(int i=0;i<NA;i++){
      sprintf(inpath,"%s/%s-grad_pr-%02d_%s.csv",prefix,gradName[i].c_str(),probeNumber,date.c_str());
      LoadGradientData(inpath,gradient);
   }

   if(rc!=0) return 1;

   // compute distance error 
   double dr[NA],dr_fxpr[NA],dr_trly[NA];
   double posErr[NA],posErrFXPR[NA],posErrTRLY[NA];
   for(int i=0;i<NA;i++){
      dr[i]         = (pp[i].dB-trly[i].dB)/gradient[i].grad;
      dr_fxpr[i]    = (pp[i].dB_fxpr-trly[i].dB_fxpr)/gradient[i].grad_fxpr;
      // dr_trly[i] = (pp[i].dB_trly-trly[i].dB_trly)/gradient[i].grad_trly;
      posErr[i]     = pp[i].dB - trly[i].dB;
      posErrFXPR[i] = pp[i].dB_fxpr - trly[i].dB_fxpr;
      posErrTRLY[i] = pp[i].dB_trly - trly[i].dB_trly;
   }

   // load in mm --> encoder count conversions
   char inpath_enc[200];
   sprintf(inpath_enc,"./input/calib-constants/pp-encoder.csv");
   std::vector<double> mmToEnc; 
   ImportPPEncoderCalib_csv(inpath_enc,mmToEnc);

   // convert to encoder counts needed to bring PP in alignment with trolley   
   double dr_enc[NA],dr_fxpr_enc[NA],dr_trly_enc[NA];
   // double sign[NA] = {-1.,-1., 1.};  // apply a sign correction to get the motion correct  
   double sign[NA] = {-1.,-1.};      // apply a sign correction to get the motion correct  
   for(int i=0;i<NA;i++){
      dr_enc[i]      = sign[i]*dr[i]*mmToEnc[i];        
      dr_fxpr_enc[i] = sign[i]*dr_fxpr[i]*mmToEnc[i]; 
      // dr_trly_enc[i] = sign[i]*dr_trly[i]*mmToEnc[i]; 
   }

   std::cout << "====================== RESULTS ======================" << std::endl; 
   std::cout << "--- Misalignment ---" << std::endl;
   std::cout << Form("x: %.3lf Hz grad = %.3lf Hz/mm dx = %.3lf mm",posErr[0],gradient[0].grad,dr[0]) << std::endl;
   std::cout << Form("y: %.3lf Hz grad = %.3lf Hz/mm dy = %.3lf mm",posErr[1],gradient[1].grad,dr[1]) << std::endl;
   // std::cout << Form("z: %.3lf Hz grad = %.3lf Hz/mm dz = %.3lf mm",posErr[2],gradient[2].grad,dr[2]) << std::endl;
   std::cout << "[Drift Cor FXPR]" << std::endl;
   std::cout << Form("x: %.3lf Hz grad = %.3lf Hz/mm dx = %.3lf mm",posErrFXPR[0],gradient[0].grad_fxpr,dr_fxpr[0]) << std::endl;
   std::cout << Form("y: %.3lf Hz grad = %.3lf Hz/mm dy = %.3lf mm",posErrFXPR[1],gradient[1].grad_fxpr,dr_fxpr[1]) << std::endl;
   // std::cout << Form("z: %.3lf Hz grad = %.3lf Hz/mm dz = %.3lf mm",posErrFXPR[2],gradient[2].grad_fxpr,dr_fxpr[2]) << std::endl;
   // std::cout << "[Drift Cor TRLY]" << std::endl;                                     
   // std::cout << Form("x: %.3lf Hz grad = %.3lf Hz/mm dx = %.3lf mm",posErrTRLY[0],dBdr[0].grad_trly,dr_trly[0]) << std::endl;
   // std::cout << Form("y: %.3lf Hz grad = %.3lf Hz/mm dy = %.3lf mm",posErrTRLY[1],dBdr[1].grad_trly,dr_trly[1]) << std::endl;
   // std::cout << Form("z: %.3lf Hz grad = %.3lf Hz/mm dz = %.3lf mm",posErrTRLY[2],dBdr[2].grad_trly,dr_trly[2]) << std::endl;
   std::cout << "--- Distance to Move ---" << std::endl;
   std::cout << Form("x: %.0lf enc counts",dr_enc[0]) << std::endl;
   std::cout << Form("y: %.0lf enc counts",dr_enc[1]) << std::endl;
   // std::cout << Form("z: %.0lf enc counts",dr_enc[2]) << std::endl;
   std::cout << "[Drift Cor FXPR]" << std::endl;
   std::cout << Form("x: %.0lf enc counts",dr_fxpr_enc[0]) << std::endl;
   std::cout << Form("y: %.0lf enc counts",dr_fxpr_enc[1]) << std::endl;
   // std::cout << Form("z: %.0lf enc counts",dr_fxpr_enc[2]) << std::endl;
   // std::cout << "[Drift Cor TRLY]" << std::endl;
   // std::cout << Form("x: %.0lf enc counts",dr_trly_enc[0]) << std::endl;
   // std::cout << Form("y: %.0lf enc counts",dr_trly_enc[1]) << std::endl;
   // std::cout << Form("z: %.0lf enc counts",dr_trly_enc[2]) << std::endl;
   std::cout << Form("-----------------------------------------------------------") << std::endl;

   return 0;
}

