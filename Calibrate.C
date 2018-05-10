// determine the calibration for a given trolley probe 

// TODO 
// 1. Add in temperature effects

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

#include "./src/InputManager.C"
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

int Calibrate(std::string configFile){

   std::cout << "----------------------------------" << std::endl;
   std::cout << "CALIBRATION CALCULATION" << std::endl;

   int rc=0;
   int method = gm2fieldUtil::Constants::kPhaseDerivative;

   std::string input_path = "./input/json/calibrate.json";
   InputManager *inputMgr = new InputManager();
   inputMgr->Load(configFile);
   inputMgr->Print();

   std::string date    = inputMgr->GetAnalysisDate();
   bool isBlind        = inputMgr->IsBlind();
   int probeNumber     = inputMgr->GetTrolleyProbe();

   std::vector<int> run;
   std::vector<std::string> runLabel;
   inputMgr->GetRunList(run); 
   inputMgr->GetRunLabels(runLabel);

   int shimRun_pp=0,shimRun_trly=0;
   int NRUN = run.size(); 
   for(int i=0;i<NRUN;i++){
      if(runLabel[i].compare("trly")==0) shimRun_trly = run[i];
      if(runLabel[i].compare("pp")==0  ) shimRun_pp   = run[i];
   } 

   date_t theDate; 
   rc = GetDate(theDate); 

   char outdir[200],outpath[200],prefix[200]; 
   if(isBlind)  sprintf(outdir ,"./output/blinded/%02d-%02d-%02d",theDate.month,theDate.day,theDate.year-2000);
   if(!isBlind) sprintf(outdir ,"./output/unblinded/%02d-%02d-%02d",theDate.month,theDate.day,theDate.year-2000);
   sprintf(prefix,"%s",outdir); 

   blind_t blind;
   ImportBlinding(blind);
   double blindErr_x  = blind.err_x; 
   double blindErr_y  = blind.err_y; 
   double blindErr_z  = blind.err_z; 
   double blindErr_pp = blind.err_pp; 
   double blindErr_tr = blind.err_tr; 

   double DR_Hz[3] = {0,0,0}; 
   double DR_Hz_blind[3] = {blindErr_x,blindErr_y,blindErr_z};
   if(isBlind){  
      for(int i=0;i<3;i++) DR_Hz[i] = DR_Hz_blind[i];
   }

   // Fixed Probe data for drift correction *between* shimmed runs
   bool trlyFirst = false;
   double driftCorShim[3],driftCorShimErr[3];
   std::vector<int> shimRun;
   // populate run vector in chornological order 
   if(shimRun_pp<shimRun_trly){
      shimRun.push_back(shimRun_pp); 
      shimRun.push_back(shimRun_trly); 
   }else{
      trlyFirst = true;
      shimRun.push_back(shimRun_trly); 
      shimRun.push_back(shimRun_pp); 
   }

   std::string fxprPath = "./input/probe-lists/fxpr-list.csv";
   std::vector<int> fxprList;
   gm2fieldUtil::Import::ImportData1<int>(fxprPath,"csv",fxprList);

   // get graph of fxpr data for drift correction *across* runs 
   double driftSign = -1.0; 
   std::vector<double> stats;
   TGraph *gDrift  = GetDriftTGraph(method,shimRun,fxprList,stats);
   double drift    = stats[2] - stats[0];  // difference of last anchor point and first anchor point of interpolation  
   double driftErr = TMath::Sqrt( stats[1]*stats[1] + stats[3]*stats[3] );  // in-quadrature sum of stdev for each anchor point  
   // double drift   = gDrift->Eval(timeBound[1]/1E+9) - gDrift->Eval(timeBound[0]/1E+9);
   if(trlyFirst) driftSign = 1;  // flip sign on the drift correction if the trolley run came first, since we calculate (PP-TRLY) 
  
   // fill the vector of drift corrections  
   driftCorShim[0] = 0;      // "raw" result 
   driftCorShim[1] = driftSign*drift;  // FXPR  
   driftCorShim[2] = driftSign*drift;  // TRLY -- take to be the same as FXPR  

   driftCorShimErr[0] = 0;  // "raw" result 
   driftCorShimErr[1] = driftErr;  // FXPR  
   driftCorShimErr[2] = driftErr;  // TRLY -- take to be the same as FXPR  

   std::string gradName[3] = {"rad","vert","azi"};
 
   char inpath[200]; 
 
   // bring in DeltaB data at final location  
   std::vector<deltab_t> pp,trly;
   for(int i=0;i<3;i++){
      sprintf(inpath,"%s/dB-pp_final-location_%s-grad_%s.csv",prefix,gradName[i].c_str(),date.c_str());
      LoadDeltaBData(inpath,pp); 
      sprintf(inpath,"%s/dB-trly_final-location_%s-grad_pr-%02d_%s.csv",prefix,gradName[i].c_str(),probeNumber,date.c_str());
      LoadDeltaBData(inpath,trly); 
   }

   // load in imposed gradients
   std::vector<grad_meas_t> gradient; 
   for(int i=0;i<3;i++){
      sprintf(inpath,"%s/%s-grad_pr-%02d_%s.csv",prefix,gradName[i].c_str(),probeNumber,date.c_str());
      LoadGradientData(inpath,gradient);
   }

   // load in gradients measured in shimmed field 
   std::vector<grad_meas_t> dBdr; 
   for(int i=0;i<3;i++){
      sprintf(inpath,"%s/%s-grad_final-location_pr-%02d_%s.csv",prefix,gradName[i].c_str(),probeNumber,date.c_str());
      rc = LoadGradientData(inpath,dBdr); 
   }

   if(rc!=0) return 1;

   // compute distance error 
   double dr[3],dr_fxpr[3],dr_trly[3];
   for(int i=0;i<3;i++){
      dr[i]      = (pp[i].dB-trly[i].dB)/gradient[i].grad;
      dr_fxpr[i] = (pp[i].dB_fxpr-trly[i].dB_fxpr)/gradient[i].grad_fxpr;
      dr_trly[i] = (pp[i].dB_trly-trly[i].dB_trly)/gradient[i].grad_trly;
   }

   // if(isBlind){
   //    for(int i=0;i<3;i++){
   //       dr[i]      = TMath::Sqrt( dr[i]*dr[i]          ); // + DR[i]*DR[i] ); 
   //       dr_fxpr[i] = TMath::Sqrt( dr_fxpr[i]*dr_fxpr[i]); // + DR[i]*DR[i] ); 
   //       dr_trly[i] = TMath::Sqrt( dr_trly[i]*dr_trly[i]); // + DR[i]*DR[i] ); 
   //    }
   // }

   // load in shimmed field data
   nmr_meas_t ppShim; 
   sprintf(inpath,"%s/pp_shimmed-field_%s.csv",prefix,date.c_str()); 
   LoadFieldData(inpath,ppShim);  

   nmr_meas_t trlyShim; 
   sprintf(inpath,"%s/trly_shimmed-field_pr-%02d_%s.csv",prefix,probeNumber,date.c_str()); 
   LoadFieldData(inpath,trlyShim);  

   double freqTr[3],freqTrErr[3];
   freqTr[0]    = trlyShim.freq; 
   freqTr[1]    = trlyShim.freq_fxpr; 
   freqTr[2]    = trlyShim.freq_trly; 

   freqTrErr[0] = trlyShim.freq_err; 
   freqTrErr[1] = trlyShim.freq_fxpr_err; 
   freqTrErr[2] = trlyShim.freq_trly_err; 

   if(isBlind){
      for(int i=0;i<3;i++) freqTrErr[i] = TMath::Sqrt( freqTrErr[i]*freqTrErr[i] + blindErr_tr*blindErr_tr );  
   }

   // load in perturbation data
   perturbation_t ppPert; 
   sprintf(inpath,"./input/perturbation/pp-pert.csv"); 
   LoadPerturbationData(inpath,ppPert);
  
   // find omega_p free (PP) 
   double freqFree[3],freqFreeErr[3];
   GetOmegaP_free(ppShim,ppPert,freqFree,freqFreeErr); 
  
   if(isBlind){
      for(int i=0;i<3;i++) freqFreeErr[i] = TMath::Sqrt( freqFreeErr[i]*freqFreeErr[i] + blindErr_pp*blindErr_pp );  
   }
 
   // compute the difference 
   double diff[3],diffErr[3]; 
   for(int i=0;i<3;i++){
      diff[i]    = freqFree[i] - (freqTr[i] + driftCorShim[i]); 
      diffErr[i] = TMath::Sqrt( freqFreeErr[i]*freqFreeErr[i] + freqTrErr[i]*freqTrErr[i] );  
   } 

   // temperature errors  
   double dT   = 0.;   // FIXME: this should be non-zero eventually  
   double dBdT = 10.*0.06179;

   // the last term here is the BLINDED values (if we're blinding)  
   double err[3],posErr[3],posErrFXPR[3],posErrTRLY[3];   
   double errTotSq=0,errTotSqFXPR=0,errTotSqTRLY=0;
   // fold in the error from the difference between PP and TRLY 
   // index is now (bare,fxpr,trly)
   errTotSq     += TMath::Power(diffErr[0],2.);
   errTotSqFXPR += TMath::Power(diffErr[1],2.);
   errTotSqTRLY += TMath::Power(diffErr[2],2.);
   // fold in temperature errors 
   errTotSq     += TMath::Power(dBdT*dT,2.);
   errTotSqFXPR += TMath::Power(dBdT*dT,2.);
   errTotSqTRLY += TMath::Power(dBdT*dT,2.);
   // adding in the error from dB/dx, dB/dy, dB/dz, dB/dT 
   // index is (x,y,z)
   for(int i=0;i<3;i++){
     // store the position error so we can look later
     posErr[i]     = TMath::Sqrt( TMath::Power(dr[i]*dBdr[i].grad,2.)           + TMath::Power(DR_Hz[i],2.) );  
     posErrFXPR[i] = TMath::Sqrt( TMath::Power(dr_fxpr[i]*dBdr[i].grad_fxpr,2.) + TMath::Power(DR_Hz[i],2.) ); 
     posErrTRLY[i] = TMath::Sqrt( TMath::Power(dr_trly[i]*dBdr[i].grad_trly,2.) + TMath::Power(DR_Hz[i],2.) );
     // fold into totals (in quadrature)  
     errTotSq     += TMath::Power(posErr[i]    ,2.);  
     errTotSqFXPR += TMath::Power(posErrFXPR[i],2.);  
     errTotSqTRLY += TMath::Power(posErrTRLY[i],2.); 
   }
   // now account for drift in dB/dr [imposed gradients]
   // index is (x,y,z) 
   std::vector<double> driftErrFXPR,driftErrTRLY;
   rc = GetDriftError(pp,trly,gradient,dBdr,driftErrFXPR,driftErrTRLY);
   for(int i=0;i<3;i++){
      errTotSqFXPR += TMath::Power(driftErrFXPR[i],2.);
      errTotSqTRLY += TMath::Power(driftErrTRLY[i],2.); 
   }

   // now fill the array    
   // add in error from drift across shimmed runs 
   err[0] = TMath::Sqrt( errTotSq ); 
   err[1] = TMath::Sqrt( errTotSqFXPR + driftCorShimErr[1]*driftCorShimErr[1] ); 
   err[2] = TMath::Sqrt( errTotSqTRLY + driftCorShimErr[2]*driftCorShimErr[2] );

   // print results to file  

   sprintf(outpath,"%s/results.csv",outdir);

   rc = MakeDirectory(outdir);

   const int NR = 3;
   double freqPPShim[NR]    = {ppShim.freq    ,ppShim.freq_fxpr    ,ppShim.freq_trly};
   double freqPPShimErr[NR] = {ppShim.freq_err,ppShim.freq_fxpr_err,ppShim.freq_trly_err};

   PrintToFile(outpath,"pp-raw"    ,NR,freqPPShim  ,freqPPShimErr  ); 
   PrintToFile(outpath,"pp-free"   ,NR,freqFree    ,freqFreeErr    ); 
   PrintToFile(outpath,"trly"      ,NR,freqTr      ,freqTrErr      ); 
   PrintToFile(outpath,"drift-shim",NR,driftCorShim,driftCorShimErr);
   PrintToFile(outpath,"diff"      ,NR,diff        ,err            );

   // print errors to file
   double DriftErrFXPR[NR],DriftErrTRLY[NR]; 
   for(int i=0;i<NR;i++){
      DriftErrFXPR[i] = driftErrFXPR[i];  
      DriftErrTRLY[i] = driftErrTRLY[i];  
   }

   // total drift error (to print to screen) 
   double driftErrFXPRTot=0,driftErrTRLYTot=0;
   for(int i=0;i<3;i++) driftErrFXPRTot += TMath::Power(DriftErrFXPR[i],2.);
   driftErrFXPRTot += TMath::Power(driftCorShimErr[1],2.); 
   driftErrFXPRTot  = TMath::Sqrt(driftErrFXPRTot); 

   for(int i=0;i<3;i++) driftErrTRLYTot += TMath::Power(DriftErrTRLY[i],2.); 
   driftErrTRLYTot += TMath::Power(driftCorShimErr[2],2.); 
   driftErrTRLYTot  = TMath::Sqrt(driftErrTRLYTot); 
 
   sprintf(outpath,"%s/results_errors.csv",outdir);  
   PrintToFile(outpath,"db-drift-fxpr",NR,DriftErrFXPR); 
   PrintToFile(outpath,"db-drift-trly",NR,DriftErrTRLY); 
   PrintToFile(outpath,"pos-raw"      ,NR,posErr); 
   PrintToFile(outpath,"pos-fxpr"     ,NR,posErrFXPR); 
   PrintToFile(outpath,"pos-trly"     ,NR,posErrTRLY); 

   // print to screen 

   std::cout << "========================= PP RESULTS ==========================" << std::endl;
   std::cout << Form("Type | No Drift Cor (Hz) | Drift Cor (Hz) [FXPR] | Drift Cor (Hz) [TRLY]") << std::endl;
   std::cout << Form("Raw  | %.3lf +/- %.3lf   | %.3lf +/- %.3lf       | %.3lf +/- %.3lf",
                             ppShim.freq,ppShim.freq_err,
                             ppShim.freq_fxpr,ppShim.freq_fxpr_err,
                             ppShim.freq_trly,ppShim.freq_trly_err) << std::endl;
   std::cout << Form("Free | %.3lf +/- %.3lf   | %.3lf +/- %.3lf       | %.3lf +/- %.3lf",
                             freqFree[0],freqFreeErr[0],
                             freqFree[1],freqFreeErr[1],
                             freqFree[2],freqFreeErr[2]) << std::endl;

   std::cout << "===================== CALIBRATION RESULTS =====================" << std::endl;
   std::cout << Form("Probe | No Drift Cor (Hz) | Drift Cor (Hz) [FXPR] | Drift Cor (Hz) [TRLY]") << std::endl;
   std::cout << Form("PP    | %.3lf +/- %.3lf | %.3lf +/- %.3lf | %.3lf +/- %.3lf",
                     freqFree[0],freqFreeErr[0],freqFree[1],freqFreeErr[1],freqFree[2],freqFreeErr[2]) << std::endl;
   std::cout << Form("TR%02d  | %.3lf +/- %.3lf | %.3lf +/- %.3lf | %.3lf +/- %.3lf",
                     probeNumber+1,freqTr[0],freqTrErr[0],freqTr[1],freqTrErr[1],freqTr[2],freqTrErr[2]) << std::endl;
   std::cout << Form("Drift (shim)  | %.3lf +/- %.3lf | %.3lf +/- %.3lf | %.3lf +/- %.3lf",
                     driftCorShim[0],driftCorShimErr[0],driftCorShim[1],driftCorShimErr[1],driftCorShim[2],driftCorShimErr[2]) << std::endl;
   std::cout << Form("-----------------------------------------------------------") << std::endl;
   std::cout << Form("PP-TR%02d | %.3lf +/- %.3lf (%.3lf ppb) | %.3lf +/- %.3lf (%.3lf ppb) | %.3lf +/- %.3lf (%.3lf ppb)",
                             probeNumber+1,diff[0],err[0],err[0]/0.06179,diff[1],err[1],err[1]/0.06179,diff[2],err[2],err[2]/0.06179) << std::endl;

   std::cout << "========================= ERROR BREAKDOWN ==========================" << std::endl;
   std::cout << "Absolute Calibration" << std::endl; 
   std::cout << Form("sigma     = %.3lf ppb",ppPert.sigma_err)     << std::endl; 
   std::cout << Form("chi       = %.3lf ppb",ppPert.chi_err)       << std::endl; 
   std::cout << Form("delta_m   = %.3lf ppb",ppPert.delta_m_err)   << std::endl; 
   std::cout << Form("delta_eps = %.3lf ppb",ppPert.delta_eps_err) << std::endl; 
   std::cout << Form("mag_img   = %.3lf ppb",ppPert.delta_mag_err) << std::endl; 
   std::cout << Form("-----------------------------------------------------------") << std::endl;
   std::cout << "Drift Errors " << std::endl;
   std::cout << "[FXPR drift]" << std::endl;
   std::cout << Form("x:   %.3lf Hz",driftErrFXPR[0]) << std::endl;
   std::cout << Form("y:   %.3lf Hz",driftErrFXPR[1]) << std::endl;
   std::cout << Form("z:   %.3lf Hz",driftErrFXPR[2]) << std::endl;
   std::cout << Form("tot: %.3lf Hz",driftErrFXPRTot) << std::endl; 
   std::cout << "[TRLY drift]" << std::endl;
   std::cout << Form("x:   %.3lf Hz",driftErrTRLY[0]) << std::endl;
   std::cout << Form("y:   %.3lf Hz",driftErrTRLY[1]) << std::endl;
   std::cout << Form("z:   %.3lf Hz",driftErrTRLY[2]) << std::endl;
   std::cout << Form("tot: %.3lf Hz",driftErrTRLYTot) << std::endl; 
   std::cout << Form("-----------------------------------------------------------") << std::endl;
   std::cout << "Misalignment" << std::endl;
   std::cout << "[No Drift Cor]" << std::endl;
   std::cout << Form("x: %.3lf Hz grad = %.3lf Hz/mm dx = %.3lf mm",posErr[0],dBdr[0].grad,dr[0]) << std::endl;
   std::cout << Form("y: %.3lf Hz grad = %.3lf Hz/mm dy = %.3lf mm",posErr[1],dBdr[1].grad,dr[1]) << std::endl;
   std::cout << Form("z: %.3lf Hz grad = %.3lf Hz/mm dz = %.3lf mm",posErr[2],dBdr[2].grad,dr[2]) << std::endl;
   std::cout << "[Drift Cor FXPR]" << std::endl;
   std::cout << Form("x: %.3lf Hz grad = %.3lf Hz/mm dx = %.3lf mm",posErrFXPR[0],dBdr[0].grad_fxpr,dr_fxpr[0]) << std::endl;
   std::cout << Form("y: %.3lf Hz grad = %.3lf Hz/mm dy = %.3lf mm",posErrFXPR[1],dBdr[1].grad_fxpr,dr_fxpr[1]) << std::endl;
   std::cout << Form("z: %.3lf Hz grad = %.3lf Hz/mm dz = %.3lf mm",posErrFXPR[2],dBdr[2].grad_fxpr,dr_fxpr[2]) << std::endl;
   std::cout << "[Drift Cor TRLY]" << std::endl;                                     
   std::cout << Form("x: %.3lf Hz grad = %.3lf Hz/mm dx = %.3lf mm",posErrTRLY[0],dBdr[0].grad_trly,dr_trly[0]) << std::endl;
   std::cout << Form("y: %.3lf Hz grad = %.3lf Hz/mm dy = %.3lf mm",posErrTRLY[1],dBdr[1].grad_trly,dr_trly[1]) << std::endl;
   std::cout << Form("z: %.3lf Hz grad = %.3lf Hz/mm dz = %.3lf mm",posErrTRLY[2],dBdr[2].grad_trly,dr_trly[2]) << std::endl;
   std::cout << Form("-----------------------------------------------------------") << std::endl;

   return 0;
}

