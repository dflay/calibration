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
// #include "./src/CalibFuncs.C"
#include "./src/FreeProton.C"

double gMarkerSize = 0.8;
double gXSize      = 0.05;  
double gYSize      = 0.06;  

int PrintResults(std::string outpath,result_prod_t result); 
int LoadSystematicUncertainties(const char *inpath,int probe,std::vector<std::string> type,double &err); 
int LoadSystematicUncertainty_final(int runPeriod,int probe,double &err); 

int ProcessResults_prod(std::string configFile){

   int rc=0;
   int method = gm2fieldUtil::Constants::kPhaseDerivative;

   InputManager *inputMgr = new InputManager();
   inputMgr->Load(configFile);
   inputMgr->Print();

   std::string ppID       = inputMgr->GetPPID(); 
   std::string blindLabel = inputMgr->GetBlindLabel();

   bool isBlind           = inputMgr->IsBlind();
   bool isMisalignCor     = inputMgr->GetMisalignCorStatus(); 

   int probeNumber        = inputMgr->GetTrolleyProbe(); 
   int runPeriod          = inputMgr->GetRunPeriod();

   // systematics 
   bool isSyst            = inputMgr->GetSystStatus();
   int systDirNum         = inputMgr->GetSystDirNum(); 

   // date_t theDate;
   // GetDate(theDate);
   std::string theDate = inputMgr->GetAnaDate();
  
   std::string outDir  = GetPath("output",isBlind,blindLabel,theDate,isSyst,systDirNum);
 
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
   double shot_err_opt = result.diffErr_aba;            // opt err is aba here -- opt only differs for misalignment   

   double shot_err_free     = result_free.diffErr; 
   double shot_err_free_aba = result_free.diffErr_aba; 
   double shot_err_free_opt = result_free.diffErr_aba;  // opt err is aba here  

   // misalignments 
   char outPath_misalign[500]; 
   // sprintf(outPath_misalign,"%s/misalignment_results_pr-%02d.csv",outDir.c_str(),probeNumber);

   // // load misalignment results 
   // misalignment_t mErr; 
   // rc = LoadMisalignmentData(outPath_misalign,mErr); 

   // double tot_misalign_err     = TMath::Sqrt(mErr.dB_x*mErr.dB_x + mErr.dB_y*mErr.dB_y + mErr.dB_z*mErr.dB_z);
   // double tot_misalign_err_aba = TMath::Sqrt(mErr.dB_x_aba*mErr.dB_x_aba + mErr.dB_y_aba*mErr.dB_y_aba + mErr.dB_z_aba*mErr.dB_z_aba);
   // double tot_misalign_err_opt = TMath::Sqrt(mErr.dB_x_opt*mErr.dB_x_opt + mErr.dB_y_opt*mErr.dB_y_opt + mErr.dB_z_opt*mErr.dB_z_opt);

   // misalignment correction
   // index 0 = raw, 1 = ABA, 2 = opt, 3 = bar-opt 
   std::vector<misalignCor_t> mCor;  
   sprintf(outPath_misalign,"%s/misalign-cor_pr-%02d.csv",outDir.c_str(),probeNumber);
   rc = LoadMisalignmentCorData(outPath_misalign,mCor);

   double tot_misalign_err         = mCor[0].err;
   double tot_misalign_err_aba     = mCor[1].err;
   double tot_misalign_err_opt     = mCor[2].err;
   double tot_misalign_err_bar_opt = mCor[3].err;

   char inpath[200];
   sprintf(inpath,"./input/perturbation/run-%d/%s.json",runPeriod,ppID.c_str());
   FreeProton *fp = new FreeProton(inpath); 

   // compute errors from free proton corrections
   double T0 = 25.0; // use T = 25 deg for the uncertainty calculation -- weak T dependence, so this is OK  
   double freeProtErr = fp->GetDelta_t_err(T0)*(0.06179/1E-9);  // converts to Hz 
   // rc = GetOmegaP_err(ppPert,freeProtErr);
   delete fp; 

   // apply the TRLY footprint correction 
   bool applyTrlyFP = inputMgr->GetTRLYFootprintStatus();  
   double trlyFP=0,trlyFPErr=0;
   if(applyTrlyFP){
      inputMgr->GetTRLYFootprint(trlyFP,trlyFPErr); 
   }
 
   // now load systematic uncertainty
   double systErr=0;
   // char inpath_syst[200];
   // sprintf(inpath_syst,"./input/json/run-%d/syst-err.json",runPeriod); 
   // std::vector<std::string> systType; 
   // systType.push_back("dB-tr"); 
   // systType.push_back("pp-freq-swap"); 
   // systType.push_back("pp-freq-dB"); 
   // systType.push_back("tr-freq-swap"); 
   // systType.push_back("tr-freq-dB");
   // rc = LoadSystematicUncertainties(inpath_syst,probeNumber-1,systType,systErr);  
   rc = LoadSystematicUncertainty_final(runPeriod,probeNumber,systErr);  

   char errStr[200],errStr_aba[200],errStr_opt[200],errStr_bar_opt[200];
   sprintf(errStr        ,"%.3lf",shot_err    ); 
   sprintf(errStr_aba    ,"%.3lf",shot_err_aba); 
   sprintf(errStr_opt    ,"%.3lf",shot_err_opt); 
   sprintf(errStr_bar_opt,"%.3lf",shot_err_opt); 

   char errStr_cor[200],errStr_cor_aba[200],errStr_cor_opt[200],errStr_cor_bar_opt[200];
   sprintf(errStr_cor        ,"%.3lf +/- %.3lf",shot_err    ,tot_misalign_err); 
   sprintf(errStr_cor_aba    ,"%.3lf +/- %.3lf",shot_err_aba,tot_misalign_err_aba); 
   sprintf(errStr_cor_opt    ,"%.3lf +/- %.3lf",shot_err_opt,tot_misalign_err_opt); 
   sprintf(errStr_cor_bar_opt,"%.3lf +/- %.3lf",shot_err_opt,tot_misalign_err_bar_opt); 

   char errStr_free[200],errStr_free_aba[200],errStr_free_opt[200];
   sprintf(errStr_free        ,"%.3lf +/- %.3lf",shot_err_free        ,freeProtErr); 
   sprintf(errStr_free_aba    ,"%.3lf +/- %.3lf",shot_err_free_aba    ,freeProtErr); 
   sprintf(errStr_free_opt    ,"%.3lf +/- %.3lf",shot_err_free_opt    ,freeProtErr);
 
   char errStr_cor_free[200],errStr_cor_free_aba[200],errStr_cor_free_opt[200],errStr_cor_free_bar_opt[200];
   sprintf(errStr_cor_free        ,"%.3lf +/- %.3lf +/- %.3lf",shot_err_free    ,tot_misalign_err        ,freeProtErr); 
   sprintf(errStr_cor_free_aba    ,"%.3lf +/- %.3lf +/- %.3lf",shot_err_free_aba,tot_misalign_err_aba    ,freeProtErr); 
   sprintf(errStr_cor_free_opt    ,"%.3lf +/- %.3lf +/- %.3lf",shot_err_free_opt,tot_misalign_err_opt    ,freeProtErr); 
   sprintf(errStr_cor_free_bar_opt,"%.3lf +/- %.3lf +/- %.3lf",shot_err_free_opt,tot_misalign_err_bar_opt,freeProtErr); 

   // apply TRLY footprint if necessary 
   if(applyTrlyFP){
      std::cout << "[ProcessResults_prod]: Applying TRLY footprint! " << std::endl;
      systErr               = TMath::Sqrt( systErr*systErr + trlyFPErr*trlyFPErr ); // update the systematic uncertainty
      result.diff                += trlyFP;
      result.diff_aba            += trlyFP;
      result.diff_opt            += trlyFP;
      result_free.diff           += trlyFP;
      result_free.diff_aba       += trlyFP;
      result_free.diff_opt       += trlyFP;
      result.diffCor             += trlyFP;
      result.diffCor_aba         += trlyFP;
      result.diffCor_opt         += trlyFP;
      result.diffCorBar_opt      += trlyFP;
      result_free.diffCor        += trlyFP;
      result_free.diffCor_aba    += trlyFP;
      result_free.diffCor_opt    += trlyFP;
      result_free.diffCorBar_opt += trlyFP;
   }

   result.systErr        = systErr; 
   result_free.systErr   = systErr;

   // update error strings 
   sprintf(errStr         ,"%s +/- %.3lf",errStr         ,systErr);  
   sprintf(errStr_aba     ,"%s +/- %.3lf",errStr_aba     ,systErr);  
   sprintf(errStr_opt     ,"%s +/- %.3lf",errStr_opt     ,systErr);  
   sprintf(errStr_free    ,"%s +/- %.3lf",errStr_free    ,systErr);  
   sprintf(errStr_free_aba,"%s +/- %.3lf",errStr_free_aba,systErr);  
   sprintf(errStr_free_opt,"%s +/- %.3lf",errStr_free_opt,systErr);  

   sprintf(errStr_cor         ,"%s +/- %.3lf",errStr_cor         ,systErr);  
   sprintf(errStr_cor_aba     ,"%s +/- %.3lf",errStr_cor_aba     ,systErr);  
   sprintf(errStr_cor_opt     ,"%s +/- %.3lf",errStr_cor_opt     ,systErr);  
   sprintf(errStr_cor_free    ,"%s +/- %.3lf",errStr_cor_free    ,systErr);  
   sprintf(errStr_cor_free_aba,"%s +/- %.3lf",errStr_cor_free_aba,systErr);  
   sprintf(errStr_cor_free_opt,"%s +/- %.3lf",errStr_cor_free_opt,systErr);  
   sprintf(errStr_cor_free_bar_opt,"%s +/- %.3lf",errStr_cor_free_bar_opt,systErr);  

   result.mErr         = tot_misalign_err; 
   result.mErr_aba     = tot_misalign_err_aba; 
   result.mErr_opt     = tot_misalign_err_opt; 
   result.mErr_bar_opt = tot_misalign_err_bar_opt; 
   result.pErr         = 0.; 
   result.pErr_aba     = 0.; 
   result.pErr_opt     = 0.; 

   result_free.mErr         = tot_misalign_err; 
   result_free.mErr_aba     = tot_misalign_err_aba; 
   result_free.mErr_opt     = tot_misalign_err_opt; 
   result_free.mErr_bar_opt = tot_misalign_err_bar_opt; 
   result_free.pErr         = freeProtErr; 
   result_free.pErr_aba     = freeProtErr; 
   result_free.pErr_opt     = freeProtErr; 
  
   std::cout << Form("******************************************************************") << std::endl;
   std::cout << Form("************************* PROBE %02d RESULTS ***********************",probeNumber) << std::endl;
   std::cout << "WITHOUT Misalignment correction" << std::endl;
   std::cout << "Bare" << std::endl;
   std::cout << Form("[RAW]: %.3lf +/- %s Hz",result.diff    ,errStr)     << std::endl;
   std::cout << Form("[ABA]: %.3lf +/- %s Hz",result.diff_aba,errStr_aba) << std::endl;
   std::cout << Form("[opt]: %.3lf +/- %s Hz",result.diff_opt,errStr_opt) << std::endl;
   std::cout << "Free proton" << std::endl;
   std::cout << Form("[RAW]: %.3lf +/- %s Hz",result_free.diff    ,errStr_free)     << std::endl;
   std::cout << Form("[ABA]: %.3lf +/- %s Hz",result_free.diff_aba,errStr_free_aba) << std::endl;
   std::cout << Form("[opt]: %.3lf +/- %s Hz",result_free.diff_opt,errStr_free_opt) << std::endl;
   std::cout << "WITH Misalignment correction" << std::endl;
   std::cout << "Bare" << std::endl;
   std::cout << Form("[RAW]:    %.3lf +/- %s Hz",result.diffCor       ,errStr_cor)         << std::endl;
   std::cout << Form("[ABA]:    %.3lf +/- %s Hz",result.diffCor_aba   ,errStr_cor_aba)     << std::endl;
   std::cout << Form("[opt]:    %.3lf +/- %s Hz",result.diffCor_opt   ,errStr_cor_opt)     << std::endl;
   std::cout << Form("[opt,bc]: %.3lf +/- %s Hz",result.diffCorBar_opt,errStr_cor_bar_opt) << std::endl;
   std::cout << "Free proton" << std::endl;
   std::cout << Form("[RAW]:    %.3lf +/- %s Hz",result_free.diffCor       ,errStr_cor_free)         << std::endl;
   std::cout << Form("[ABA]:    %.3lf +/- %s Hz",result_free.diffCor_aba   ,errStr_cor_free_aba)     << std::endl;
   std::cout << Form("[opt]:    %.3lf +/- %s Hz",result_free.diffCor_opt   ,errStr_cor_free_opt)     << std::endl;
   std::cout << Form("[opt,bc]: %.3lf +/- %s Hz",result_free.diffCorBar_opt,errStr_cor_free_bar_opt) << std::endl;
   std::cout << Form("******************************************************************") << std::endl;
   std::cout << Form("******************************************************************") << std::endl;

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
int LoadSystematicUncertainty_final(int runPeriod,int probe,double &err){
   // load final systematic uncertainties
   // uses the CSV manager 
   char prefix[200],inpath[200]; 
   sprintf(prefix,"./input/systematic-errors/run-%d",runPeriod);

   CSVManager *csvMgr = new CSVManager(); 
 
   // get frequency extraction uncertainties 
   bool headerStatus = true; 
   sprintf(inpath,"%s/freq-ext.csv",prefix);
   int rc = csvMgr->ReadFile(inpath,headerStatus);
   if(rc!=0) return rc; 

   std::vector<double> freq_misalignErr,freq_swapErr;
   rc = csvMgr->GetColumn_byName<double>("Tot_Corr_sys_freq",freq_misalignErr); 
   if(rc!=0) return rc; 
   rc = csvMgr->GetColumn_byName<double>("calibration_sys" ,freq_swapErr); 
   if(rc!=0) return rc; 

   csvMgr->ClearData(); 

   // now cut data 
   std::vector<double> dBtrErr,rapidSwapErr;
   sprintf(inpath,"%s/ana-cuts.csv",prefix);
   csvMgr->ReadFile(inpath,headerStatus); 
   rc = csvMgr->GetColumn_byName<double>("dB_tr"     ,dBtrErr); 
   if(rc!=0) return rc; 
   rc = csvMgr->GetColumn_byName<double>("rapid_swap",rapidSwapErr);
   if(rc!=0) return rc; 

   delete csvMgr;  
  
   // now take quadrature sum 
   int k = probe-1;
   double sum = TMath::Power(freq_misalignErr[k],2.) + TMath::Power(freq_swapErr[k],2.) 
              + TMath::Power(dBtrErr[k],2.)          + TMath::Power(rapidSwapErr[k],2.);
   err = TMath::Sqrt(sum);
   return 0;  
}
//______________________________________________________________________________
int LoadSystematicUncertainties(const char *inpath,int probe,std::vector<std::string> type,double &err){

  std::string path = inpath; 
  json data; 
  int rc = gm2fieldUtil::Import::ImportJSON(inpath,data);

  double arg=0,sum_sq=0; 
  const int N = 17;           // number of probes
  const int NT = type.size(); // number of types 
  for(int j=0;j<NT;j++){
     arg     = (double)data[type[j]][probe]; 
     sum_sq += arg*arg;
  }
  err = TMath::Sqrt(sum_sq); 
  return 0;
}
//______________________________________________________________________________
int PrintResults(std::string outpath,result_prod_t result){

   char header[500]; 
   // prepare the header 
   sprintf(header,"#type,diff,diff-mis-cor,shot-err,misalign-err,free-prot-err,syst"); 
   char myStr[1000]; 
   sprintf(myStr    ,"raw,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf",result.diff    ,result.diffCor    ,result.diffErr    ,result.mErr    ,result.pErr    ,result.systErr); 
   char myStr_aba[1000]; 
   sprintf(myStr_aba,"ABA,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf",result.diff_aba,result.diffCor_aba,result.diffErr_aba,result.mErr_aba,result.pErr_aba,result.systErr); 
   char myStr_opt[1000]; 
   sprintf(myStr_opt,"opt,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf",result.diff_opt,result.diffCor_opt,result.diffErr_opt,result.mErr_opt,result.pErr_opt,result.systErr); 
   char myStr_bar_opt[1000]; 
   sprintf(myStr_bar_opt,"optbc,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf",result.diff_opt,result.diffCorBar_opt,result.diffErr_opt,result.mErr_bar_opt,result.pErr_opt,result.systErr); 

   std::ofstream outfile;
   outfile.open(outpath.c_str());
   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      outfile << header        << std::endl; 
      outfile << myStr         << std::endl;
      outfile << myStr_aba     << std::endl;
      outfile << myStr_opt     << std::endl;
      outfile << myStr_bar_opt << std::endl;
      std::cout << "The data has been written to file: " << outpath << std::endl;
      outfile.close(); 
   }  
   return 0;
}
