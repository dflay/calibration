// Make tables of all results and print ROOT and JSON files that contain everything 

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
#include "gm2fieldImport.h"
#include "gm2fieldExport.h"
#include "TemperatureSensor.h"
#include "MovingAverage.h"

#include "./include/date.h"
#include "./include/Constants.h"
#include "./include/fixedProbeEvent.h"
#include "./include/trolleyAnaEvent.h"
#include "./include/nmr_meas.h"
#include "./include/grad_meas.h"
#include "./include/deltab.h"
#include "./include/results.h"
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

double gMarkerSize = 0.8;
double gXSize      = 0.05;  
double gYSize      = 0.06;  

int GetTemperatures(const char *outdir,int probeNumber,result_prod_t &res,result_prod_t &resFree); 
int GetFreeProtonCorrection(const char *outdir,int probeNumber,result_prod_t &resFree); 
int GetShimmedGradientFitPars(const char *outdir,int probeNumber,std::vector<grad_meas_t> &data); 
int GetMisalignmentCorByAxisData(const char *outdir,int probeNumber,std::vector<misalignCor_t> &data); 

int ApplyShimmedGradEstimates(int probe,std::string prodVersion,std::vector<grad_meas_t> &data); 
int PrintTableOfResults(const char *outpath,int probe,double diff,double diffFree,
                        double swapErr,double mErr,double pErr,double systErr);

int PrintTableOfResults_db(const char *outpath,int probe,
                           double dBx_pp,double dBx_pp_err,double dBx_tr,double dBx_tr_err,
                           double dBy_pp,double dBy_pp_err,double dBy_tr,double dBy_tr_err,
                           double dBz_pp,double dBz_pp_err,double dBz_tr,double dBz_tr_err); 

int PrintTableOfResults_grad(const char *outpath,int probe,double dBdx,double dBdx_err,
                             double dBdy,double dBdy_err,double dBdz,double dBdz_err);

int PrintTableOfResults_mis(const char *outpath,int probe,double dx,double dB_x,
                            double dy,double dB_y,double dz,double dB_z);

int GetImposedGradients(const char *inpath,std::vector<imposed_gradient_t> &imp_grad); 
int GetShimmedGradients(const char *prefix,int probe,std::string prodVersion,std::vector<std::string> gradName,
                        std::vector<grad_meas_t> &shim_grad);

int PrintResults(const char *outpath,std::vector<calib_result_t> data,bool toROOT=true,bool toJSON=true,bool toCSV=true); 

int PrintToFile_csv(const char *outpath,std::vector<calib_result_t> data); 

int PrintToFile_json(const char *outpath,std::vector<calib_result_t> data); 
int GetJSONObject(std::vector<calib_result_t> data,json &jData);

 int CollectData(std::vector<result_prod_t> r,std::vector<result_prod_t> rFree,
                std::vector<std::vector<imposed_gradient_t>> ig,std::vector<std::vector<grad_meas_t> > mg,
                std::vector<deltab_prod_t> dB_pp,std::vector<deltab_prod_t> dB_tr,
                std::vector<misalignment_t> mis,std::vector< std::vector<misalignCor_t> > mCor,
                std::vector< std::vector<misalignCor_t> > mCor_a,std::vector<calib_result_t> &data); 

int MakeTables_prod(int runPeriod,std::string theDate,int isSyst,int systDirNum){

   int rc=0;

   json configData; 
   char inpath_config[200];
   sprintf(inpath_config,"./input/json/run-%d/config.json",runPeriod); 
   std::string cPath = inpath_config; 
   rc = gm2fieldUtil::Import::ImportJSON(cPath,configData);

   // config data 
   std::string prodVersion = configData["prod-tag"]; 
   std::string blindLabel  = configData["blinding"]["label"]; 
   bool isBlind            = configData["blinding"]["enable"];  
   bool isMisCor           = (bool)( (int)configData["use-misalign-cor"] );

   std::string outDir = GetPath("output",isBlind,blindLabel,theDate,isSyst,systDirNum); 

   // char outDir[200];
   // if(isBlind)  sprintf(outDir,"./output/blinded/%s",blindLabel.c_str());
   // if(!isBlind) sprintf(outDir,"./output/unblinded");
   // sprintf(outDir,"%s/%s",outDir,theDate.c_str()); 

   char outPath[200]; 
   sprintf(outPath,"%s/calibData_%s",outDir.c_str(),theDate.c_str()); 

   result_prod_t result;
   std::vector<result_prod_t> res,resFree;

   std::vector<imposed_gradient_t> ig;   
   std::vector< std::vector<imposed_gradient_t> > impGrad;   

   deltab_prod_t dbPP,dbTR;  
   std::vector<deltab_prod_t> deltaB_pp,deltaB_trly;

   std::vector<grad_meas_t> grad; 
   std::vector< std::vector<grad_meas_t> > shimGrad;

   misalignment_t m; 
   std::vector<misalignment_t> misalign;

   std::vector<misalignCor_t> mc,mca; 
   std::vector< std::vector<misalignCor_t> > mCor,mCor_a; 

   char inpath[500]; 

   std::vector<std::string> gradName;
   gradName.push_back("rad");  
   gradName.push_back("vert");  
   gradName.push_back("azi");  

   std::vector<int> probe; 

   bool update = false;

   double syst=0;
   double sum_sq=0;
   double sumf=0,sum_sf=0,sum=0;
   double sumf_ppb=0,sum_ppb=0;
   double ppt=0,ppt_err=0; 
   double trt=0,trt_err=0; 

   int k=0; 
   int probeNumber = 0; 
   const int NTRLY = 17; 
   for(int i=0;i<NTRLY;i++){
      probeNumber = i+1; 
      std::cout << Form("LOADING PROBE %02d",probeNumber) << std::endl;
      // load raw result 
      sprintf(inpath,"%s/results_final_pr-%02d.csv",outDir.c_str(),probeNumber);
      rc = LoadResultsProdFinalData(inpath,result);
      if(rc!=0) continue;
      res.push_back(result); 
      // load free-proton result 
      sprintf(inpath,"%s/results_final_free-prot_pr-%02d.csv",outDir.c_str(),probeNumber);
      rc = LoadResultsProdFinalData(inpath,result);
      resFree.push_back(result);
      // load temperatures 
      rc = GetTemperatures(outDir.c_str(),probeNumber,res[i],resFree[i]);
      // load free proton correction 
      rc = GetFreeProtonCorrection(outDir.c_str(),probeNumber,resFree[i]);
      if(rc!=0) continue;
      // increment probe number
      probe.push_back(probeNumber);  
      // load OPTIMIZED Delta-B numbers
      sprintf(inpath,"%s/delta-b-opt_pr-%02d.csv",outDir.c_str(),probeNumber); 
      rc = LoadDeltaB_opt(inpath,dbPP,dbTR);
      deltaB_pp.push_back(dbPP);
      deltaB_trly.push_back(dbTR); 
      // load shimmed gradients 
      sprintf(inpath,"%s/shim-grad_opt_pr-%02d.csv",outDir.c_str(),probeNumber); 
      rc = LoadShimmedGrad_opt(inpath,grad);  
      // load shimmed field fit pars
      rc = GetShimmedGradientFitPars(outDir.c_str(),probeNumber,grad); 
      if(rc!=0) return 1;
      shimGrad.push_back(grad);
      // clean up 
      grad.clear(); 
      // get imposed gradients 
      sprintf(inpath,"%s/imposed-gradients_pr-%02d.csv",outDir.c_str(),probeNumber); 
      rc = GetImposedGradients(inpath,ig); 
      if(rc!=0) return 1;
      impGrad.push_back(ig); 
      // clean up 
      ig.clear(); 
      // load misalignment data 
      sprintf(inpath,"%s/misalignment_results_pr-%02d.csv",outDir.c_str(),probeNumber);
      rc = LoadMisalignmentData(inpath,m);
      if(rc!=0) return 1;
      misalign.push_back(m);
      // load misalignment CORRECTION data 
      sprintf(inpath,"%s/misalign-cor_pr-%02d.csv",outDir.c_str(),probeNumber); 
      rc = LoadMisalignmentCorData(inpath,mc); 
      if(rc!=0) return 1;
      mCor.push_back(mc);
      // misalignment corrections by axis
      rc = GetMisalignmentCorByAxisData(outDir.c_str(),probeNumber,mca); 
      if(rc!=0) return 1;
      mCor_a.push_back(mca); 
      // clean up 
      mc.clear();
      mca.clear();  
   }

   // now make tables

   // also print table of results to file 
   char outpath_res[200];
   sprintf(outpath_res,"%s/results_all-probes.csv",outDir.c_str());

   char outpath_res_free[200];
   sprintf(outpath_res_free,"%s/results_free-prot_all-probes.csv",outDir.c_str());

   char outpath_db[200];
   sprintf(outpath_db,"%s/results_db_all-probes.csv",outDir.c_str());

   char outpath_ig[200];
   sprintf(outpath_ig,"%s/results_imposed-gradients_all-probes.csv",outDir.c_str());

   char outpath_sg[200];
   sprintf(outpath_sg,"%s/results_shimmed-gradients_all-probes.csv",outDir.c_str());

   char outpath_mis[200];
   sprintf(outpath_mis,"%s/results_misalignments_all-probes.csv",outDir.c_str());

   const int NPROBES = probe.size(); 

   // PP-TRLY 
   std::cout << "===================================== PP-TRLY =====================================" << std::endl;
   for(int i=0;i<NPROBES;i++){
      sum_sq   = res[i].diffErr_aba*res[i].diffErr_aba + res[i].mErr_opt*res[i].mErr_opt + res[i].systErr*res[i].systErr;
      sum      = TMath::Sqrt(sum_sq);  
      sum_sq   = resFree[i].diffErr_aba*resFree[i].diffErr_aba + res[i].mErr_opt*res[i].mErr_opt + resFree[i].pErr*resFree[i].pErr 
                 + resFree[i].systErr*resFree[i].systErr;
      sumf     = TMath::Sqrt(sum_sq); 
      sumf_ppb = sumf/0.06179;  
      sum_sq   = resFree[i].diffErr_aba*resFree[i].diffErr_aba + res[i].mErr_opt*res[i].mErr_opt + resFree[i].systErr*resFree[i].systErr;
      sum_sf   = TMath::Sqrt(sum_sq);  
      std::cout << Form("%02d,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf",
                        probe[i],res[i].diff_aba,resFree[i].diff_aba,res[i].diffErr_aba,
                        resFree[i].diffErr_aba,res[i].mErr_opt,resFree[i].pErr,res[i].systErr,sum,sum_sf,sumf_ppb) << std::endl;
      rc = PrintTableOfResults(outpath_res,probe[i],res[i].diff_aba,resFree[i].diff_aba,res[i].diffErr_aba,res[i].mErr,resFree[i].pErr,res[i].systErr);
   }

   std::cout << "===================================== DELTA B =====================================" << std::endl;
   for(int i=0;i<NPROBES;i++){
      std::cout << Form("%02d,%.3lf ± %.3lf,%.3lf ± %.3lf,%.3lf ± %.3lf,%.3lf ± %.3lf,%.3lf ± %.3lf,%.3lf ± %.3lf",probe[i],
                        deltaB_pp[i].dB_opt[0],deltaB_pp[i].dB_opt_err[0],deltaB_trly[i].dB_opt[0],deltaB_trly[i].dB_opt_err[0],
                        deltaB_pp[i].dB_opt[1],deltaB_pp[i].dB_opt_err[1],deltaB_trly[i].dB_opt[1],deltaB_trly[i].dB_opt_err[1],
                        deltaB_pp[i].dB_opt[2],deltaB_pp[i].dB_opt_err[2],deltaB_trly[i].dB_opt[2],deltaB_trly[i].dB_opt_err[2]) << std::endl;

   rc = PrintTableOfResults_db(outpath_db,probe[i],
                               deltaB_pp[i].dB_opt[0],deltaB_pp[i].dB_opt_err[0],deltaB_trly[i].dB_opt[0],deltaB_trly[i].dB_opt_err[0],
                               deltaB_pp[i].dB_opt[1],deltaB_pp[i].dB_opt_err[1],deltaB_trly[i].dB_opt[1],deltaB_trly[i].dB_opt_err[1],
                               deltaB_pp[i].dB_opt[2],deltaB_pp[i].dB_opt_err[2],deltaB_trly[i].dB_opt[2],deltaB_trly[i].dB_opt_err[2]);
   }

   std::cout << "===================================== MISALIGNMENTS =====================================" << std::endl;
   for(int i=0;i<NPROBES;i++){
      sum_sq = misalign[i].dB_x_opt*misalign[i].dB_x_opt + misalign[i].dB_y_opt*misalign[i].dB_y_opt + misalign[i].dB_z_opt*misalign[i].dB_z_opt;
      sum    = TMath::Sqrt(sum_sq); 
      std::cout << Form("%02d,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf",probe[i],
                        misalign[i].dx_opt,misalign[i].dB_x_opt,
                        misalign[i].dy_opt,misalign[i].dB_y_opt,
                        misalign[i].dz_opt,misalign[i].dB_z_opt,
                        sum) << std::endl;

      rc = PrintTableOfResults_mis(outpath_mis,probe[i],
                                   misalign[i].dx_opt,misalign[i].dB_x_opt,
                                   misalign[i].dy_opt,misalign[i].dB_y_opt,
                                   misalign[i].dz_opt,misalign[i].dB_z_opt);
   }
   
   std::cout << "===================================== MISALIGNMENT CORRECTION =====================================" << std::endl;
   // note: we're printing the opt result! 
   for(int i=0;i<NPROBES;i++){
      std::cout << Form("%02d,%.3lf ± %.3lf",i+1,mCor[i][2].val,mCor[i][2].err) << std::endl;
   }
   
   std::cout << "===================================== SHIMMED GRADIENTS =====================================" << std::endl;
   for(int i=0;i<NPROBES;i++){
      std::cout << Form("%02d,%.3lf ± %.3lf,%.3lf ± %.3lf,%.3lf ± %.3lf",
                        probe[i],shimGrad[i][0].grad,shimGrad[i][0].grad_err,
                        shimGrad[i][1].grad,shimGrad[i][1].grad_err,
                        shimGrad[i][2].grad,shimGrad[i][2].grad_err) << std::endl;
      rc = PrintTableOfResults_grad(outpath_sg,probe[i],
                                    shimGrad[i][0].grad,shimGrad[i][0].grad_err,
	                            shimGrad[i][1].grad,shimGrad[i][1].grad_err,
	                            shimGrad[i][2].grad,shimGrad[i][2].grad_err); 
   }

   std::cout << "===================================== IMPOSED GRADIENTS =====================================" << std::endl;
   for(int i=0;i<NPROBES;i++){
      std::cout << Form("%02d,%.3lf,%.3lf,%.3lf",
                        probe[i],impGrad[i][0].grad,impGrad[i][1].grad,impGrad[i][2].grad) << std::endl;
      rc = PrintTableOfResults_grad(outpath_ig,probe[i],
                                    impGrad[i][0].grad,0.0,
	                            impGrad[i][1].grad,0.0,
	                            impGrad[i][2].grad,0.0);
   }

   std::vector<calib_result_t> data; 
   rc = CollectData(res,resFree,impGrad,shimGrad,deltaB_pp,deltaB_trly,misalign,mCor,mCor_a,data); 

   rc = PrintResults(outPath,data); 

   return 0;
}
//______________________________________________________________________________
int CollectData(std::vector<result_prod_t> r,std::vector<result_prod_t> rFree,
                std::vector<std::vector<imposed_gradient_t>> ig,std::vector<std::vector<grad_meas_t> > mg,
                std::vector<deltab_prod_t> dB_pp,std::vector<deltab_prod_t> dB_tr,
                std::vector<misalignment_t> mis,std::vector< std::vector<misalignCor_t> > mCor,
                std::vector< std::vector<misalignCor_t> > mCor_a,std::vector<calib_result_t> &data){

   // print the results to a ROOT file
   calib_result_t dataPt; 
   const int N = r.size();
   for(int i=0;i<N;i++){
      // no proton corrections, no misalignment 
      dataPt.calibCoeff               = r[i].diff; 
      dataPt.calibCoeffErr            = r[i].diffErr;  
      dataPt.calibCoeff_aba           = r[i].diff_aba; 
      dataPt.calibCoeffErr_aba        = r[i].diffErr_aba;
      dataPt.calibCoeff_opt           = r[i].diff_opt; 
      dataPt.calibCoeffErr_opt        = r[i].diffErr_opt;
      // no proton correction, with misalignment
      dataPt.calibCoeffCor            = r[i].diffCor; 
      dataPt.calibCoeffCorErr         = r[i].diffErr;         // note uncertainty is same as without mis cor!
      dataPt.calibCoeffCor_aba        = r[i].diffCor_aba;      
      dataPt.calibCoeffCorErr_aba     = r[i].diffErr_aba;     // note uncertainty is same as without mis cor!
      dataPt.calibCoeffCor_opt        = r[i].diffCor_opt;      
      dataPt.calibCoeffCorErr_opt     = r[i].diffErr_opt;     // note uncertainty is same as without mis cor!
      // with free proton, without misalignment 
      dataPt.calibCoeffFree           = rFree[i].diff;  
      dataPt.calibCoeffFreeErr        = rFree[i].diffErr;        
      dataPt.calibCoeffFree_aba       = rFree[i].diff_aba;    
      dataPt.calibCoeffFreeErr_aba    = rFree[i].diffErr_aba;   
      dataPt.calibCoeffFree_opt       = rFree[i].diff_opt;    
      dataPt.calibCoeffFreeErr_opt    = rFree[i].diffErr_opt;  
      // with free proton, with misalignment 
      dataPt.calibCoeffCorFree        = rFree[i].diffCor;  
      dataPt.calibCoeffCorFreeErr     = rFree[i].diffErr;     // note uncertainty is same as without mis cor! 
      dataPt.calibCoeffCorFree_aba    = rFree[i].diffCor_aba;  
      dataPt.calibCoeffCorFreeErr_aba = rFree[i].diffErr_aba; // note uncertainty is same as without mis cor!  
      dataPt.calibCoeffCorFree_opt    = rFree[i].diffCor_opt;  
      dataPt.calibCoeffCorFreeErr_opt = rFree[i].diffErr_opt; // note uncertainty is same as without mis cor! 
      dataPt.freeErr                  = rFree[i].pErr; 
      // free-proton related details.  
      dataPt.fpCor                 = rFree[i].freeProtCor; 
      dataPt.fpCorErr              = rFree[i].freeProtCorErr;
      dataPt.ppTemp                = rFree[i].ppTemp;  
      dataPt.ppTempErr             = rFree[i].ppTempErr;  
      dataPt.trTemp                = rFree[i].trTemp;  
      dataPt.trTempErr             = rFree[i].trTempErr;  
      // imposed gradient data 
      dataPt.dBdx_imp              = ig[i][0].grad;  
      dataPt.dBdy_imp              = ig[i][1].grad;
      dataPt.dBdz_imp              = ig[i][2].grad; 
      dataPt.dBdx_impErr           = ig[i][0].grad_err;  
      dataPt.dBdy_impErr           = ig[i][1].grad_err;
      dataPt.dBdz_impErr           = ig[i][2].grad_err; 
      // shimmed gradient data 
      dataPt.dBdx_shim             = mg[i][0].grad;  
      dataPt.dBdy_shim             = mg[i][1].grad;
      dataPt.dBdz_shim             = mg[i][2].grad; 
      dataPt.dBdx_shimErr          = mg[i][0].grad_err;  
      dataPt.dBdy_shimErr          = mg[i][1].grad_err;
      dataPt.dBdz_shimErr          = mg[i][2].grad_err; 
      // shimmed field fit parameters -- note REVERSED order of par index vs x,y,z 
      dataPt.shim_x_a              = mg[i][0].par[2]; 
      dataPt.shim_x_aErr           = mg[i][0].parErr[2]; 
      dataPt.shim_x_b              = mg[i][0].par[1]; 
      dataPt.shim_x_bErr           = mg[i][0].parErr[1]; 
      dataPt.shim_x_c              = mg[i][0].par[0]; 
      dataPt.shim_x_cErr           = mg[i][0].parErr[0]; 
      dataPt.shim_y_a              = mg[i][1].par[2]; 
      dataPt.shim_y_aErr           = mg[i][1].parErr[2]; 
      dataPt.shim_y_b              = mg[i][1].par[1]; 
      dataPt.shim_y_bErr           = mg[i][1].parErr[1]; 
      dataPt.shim_y_c              = mg[i][1].par[0]; 
      dataPt.shim_y_cErr           = mg[i][1].parErr[0]; 
      dataPt.shim_z_a              = mg[i][2].par[2]; 
      dataPt.shim_z_aErr           = mg[i][2].parErr[2]; 
      dataPt.shim_z_b              = mg[i][2].par[1]; 
      dataPt.shim_z_bErr           = mg[i][2].parErr[1]; 
      dataPt.shim_z_c              = mg[i][2].par[0]; 
      dataPt.shim_z_cErr           = mg[i][2].parErr[0]; 
      // delta-B data  
      dataPt.deltaB_pp_x           = dB_pp[i].dB_opt[0]; 
      dataPt.deltaB_pp_y           = dB_pp[i].dB_opt[1]; 
      dataPt.deltaB_pp_z           = dB_pp[i].dB_opt[2]; 
      dataPt.deltaB_pp_xErr        = dB_pp[i].dB_opt_err[0]; 
      dataPt.deltaB_pp_yErr        = dB_pp[i].dB_opt_err[1]; 
      dataPt.deltaB_pp_zErr        = dB_pp[i].dB_opt_err[2]; 
      dataPt.deltaB_tr_x           = dB_tr[i].dB_opt[0]; 
      dataPt.deltaB_tr_y           = dB_tr[i].dB_opt[1]; 
      dataPt.deltaB_tr_z           = dB_tr[i].dB_opt[2]; 
      dataPt.deltaB_tr_xErr        = dB_tr[i].dB_opt_err[0]; 
      dataPt.deltaB_tr_yErr        = dB_tr[i].dB_opt_err[1]; 
      dataPt.deltaB_tr_zErr        = dB_tr[i].dB_opt_err[2];
      // misalignment data.  note that this is in Hz 
      dataPt.dB_x                  = mis[i].dB_x_opt; 
      dataPt.dB_y                  = mis[i].dB_y_opt; 
      dataPt.dB_z                  = mis[i].dB_z_opt;
      dataPt.dB_r_tot              = TMath::Sqrt( dataPt.dB_x*dataPt.dB_x + dataPt.dB_y*dataPt.dB_y + dataPt.dB_z*dataPt.dB_z ); 
      // misalignment CORRECTION. always use the opt result.  this is in Hz 
      // this is a sum over all axes
      dataPt.misCor                = mCor[i][2].val; // index 2 = opt result     
      dataPt.misCor_err            = mCor[i][2].err; // index 2 = opt result 
      // now the second index is AXIS
      dataPt.misCor_x              = mCor_a[i][0].val;  
      dataPt.misCor_xErr           = mCor_a[i][0].err;
      dataPt.misCor_y              = mCor_a[i][1].val;
      dataPt.misCor_yErr           = mCor_a[i][1].err;
      dataPt.misCor_z              = mCor_a[i][2].val;
      dataPt.misCor_zErr           = mCor_a[i][2].err;
      dataPt.misCor_z_bar          = 0.;
      dataPt.misCor_zErr_bar       = 0.;
      // misalignment in mm; TODO: uncertainties? 
      dataPt.mis_x                 = mis[i].dx_opt;  
      dataPt.mis_xErr              = 0;  
      dataPt.mis_y                 = mis[i].dy_opt;  
      dataPt.mis_yErr              = 0;  
      dataPt.mis_z                 = mis[i].dz_opt;  
      dataPt.mis_zErr              = 0;  
      dataPt.mis_z_bar             = 0;  
      dataPt.mis_zErr_bar          = 0;  
      // systematic uncertainty (this is identical for raw and free results)  
      dataPt.systErr               = r[i].systErr;  
      // fill vector 
      data.push_back(dataPt);  
   }
   return 0;
}
//______________________________________________________________________________
int GetMisalignmentCorByAxisData(const char *outdir,int probeNumber,std::vector<misalignCor_t> &data){

   misalignCor_t inData;
   std::vector<double> v,ve;
   char inpath[200];
   sprintf(inpath,"%s/misalign-cor-opt_by-axis_pr-%02d.csv",outdir,probeNumber); 
   int rc = gm2fieldUtil::Import::ImportData2<double>(inpath,"csv",v,ve);

   std::string axis[3] = {"x","y","z"}; 

   // size is number of axes = 3
   const int N = v.size();
   for(int i=0;i<N;i++){
      inData.label = axis[i]; 
      inData.val   = v[i];
      inData.err   = ve[i];
      // std::cout << Form("probe %02d, misCor_%s = %.3lf +/- %.3lf",probeNumber,axis[i].c_str(),inData.val,inData.err) << std::endl;
      data.push_back(inData);
   }

   return 0; 

}
//______________________________________________________________________________
int GetFreeProtonCorrection(const char *outdir,int probeNumber,result_prod_t &resFree){

   char path[200];
   sprintf(path,"%s/pp-swap_fpc_pr-%02d.csv",outdir,probeNumber);
   std::string inpath = path; 

   std::vector<double> fpc,fpce;
   int rc = gm2fieldUtil::Import::ImportData2<double>(inpath,"csv",fpc,fpce);
   if(rc!=0) return 1;

   // WARNING: we use the error on the mean, since that's effectively the weighted 
   // average on the uncertainties that we use as weights.  Those uncertainties 
   // are more meaningful than the scatter in the data here
   double mean=0,err=0; 
   std::vector<double> w; 
   const int N = fpc.size();
   for(int i=0;i<N;i++) w.push_back(1./fpce[i]); 
 
   rc    = gm2fieldUtil::Math::GetWeightedMean<double>(fpc,w,mean,err); 
   // stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(fpc); 

   // for systematic uncertainties, we take the MEAN of all uncertainties; 
   // averaging over a systematic value doesn't make it better.  
   double meanErr = gm2fieldUtil::Math::GetMean<double>(fpce);   

   resFree.freeProtCor    = mean;     
   resFree.freeProtCorErr = meanErr;    

   return 0; 
}
//______________________________________________________________________________
int GetShimmedGradientFitPars(const char *outdir,int probeNumber,std::vector<grad_meas_t> &data){

  std::vector<double> xPar,xParErr;
  std::vector<double> yPar,yParErr;
  std::vector<double> zPar,zParErr;

  char inpath[200]; 
  sprintf(inpath,"%s/shimmed-grad-x_pars_pr-%02d.csv",outdir,probeNumber);
  LoadShimmedFieldFitPars(inpath,xPar,xParErr); 
  sprintf(inpath,"%s/shimmed-grad-y_pars_pr-%02d.csv",outdir,probeNumber);
  LoadShimmedFieldFitPars(inpath,yPar,yParErr); 
  sprintf(inpath,"%s/shimmed-grad-z_pars_pr-%02d.csv",outdir,probeNumber);
  LoadShimmedFieldFitPars(inpath,zPar,zParErr); 

  // note the dimension of data is 3 -- x, y, z
  // in case we have different dimensions for each axis, we do them separately  
  int NPAR = xPar.size();
  for(int i=0;i<NPAR;i++){ 
     data[0].par[i]    = xPar[i];
     data[0].parErr[i] = xParErr[i];
  }
  NPAR = yPar.size();
  for(int i=0;i<NPAR;i++){ 
     data[1].par[i]    = yPar[i];
     data[1].parErr[i] = yParErr[i];
  }
  NPAR = zPar.size();
  for(int i=0;i<NPAR;i++){ 
     data[2].par[i]    = zPar[i];
     data[2].parErr[i] = zParErr[i];
  }

  return 0;

}
//______________________________________________________________________________
int GetTemperatures(const char *outdir,int probeNumber,result_prod_t &res,result_prod_t &resFree){
   // load PP and TRLY temperatures
   int rc=0;

   char inpath[200]; 
   std::vector<calibSwap_t> pp,tr;
   sprintf(inpath,"%s/pp-swap-data_pr-%02d.csv",outdir,probeNumber);
   rc = LoadCalibSwapData(inpath,pp);
   if(rc!=0) return 1;
   sprintf(inpath,"%s/trly-swap-data_pr-%02d.csv",outdir,probeNumber);
   rc = LoadCalibSwapData(inpath,tr);
   if(rc!=0) return 1;

   // get temps 
   const int N = pp.size();
   std::vector<double> x,y;
   for(int i=0;i<N;i++){
      x.push_back(pp[i].temp);
      y.push_back(tr[i].temp);
   }

   // PP
   double ppt     = gm2fieldUtil::Math::GetMean<double>(x); 
   double ppt_err = gm2fieldUtil::Math::GetStandardDeviation<double>(x); 
   // TRLY
   double trt     = gm2fieldUtil::Math::GetMean<double>(y); 
   double trt_err = gm2fieldUtil::Math::GetStandardDeviation<double>(y);

   // apply temperatures 
   res.ppTemp        = ppt;  
   res.ppTempErr     = ppt_err;  
   res.trTemp        = trt;  
   res.trTempErr     = trt_err;  
   resFree.ppTemp    = ppt;  
   resFree.ppTempErr = ppt_err;  
   resFree.trTemp    = trt;  
   resFree.trTempErr = trt_err;  

   return 0; 
   
}
//______________________________________________________________________________
int PrintResults(const char *filename,std::vector<calib_result_t> data,bool toROOT,bool toJSON,bool toCSV){
   // print results to multiple formats 

   char outpath_root[200],outpath_json[200],outpath_csv[200];
   sprintf(outpath_root,"%s.root",filename); 
   sprintf(outpath_json,"%s.json",filename); 
   sprintf(outpath_csv ,"%s.csv" ,filename); 

   // ROOT file parameters 
   gm2fieldUtil::rootData_t rd;
   rd.fileName      = outpath_root;
   rd.treeName      = "CAL";
   rd.branchName    = "B";
   rd.leafStructure = calib_result_str; 

   int rc=0;
   if(toROOT) rc = gm2fieldUtil::Export::PrintToROOTFile<calib_result_t>(rd,data); 
   if(toJSON) rc = PrintToFile_json(outpath_json,data); 
   if(toCSV)  rc = PrintToFile_csv(outpath_csv,data); 
   
   return rc;
}
//______________________________________________________________________________
int PrintToFile_csv(const char *outpath,std::vector<calib_result_t> data){

   char outStr[1000]; 
   const int N = data.size();

   std::string header = "Probe,calibCoeff,calibCoeffErr,calibCoeff_cor,calibCoeffErr_cor,";
   header            += "calibCoeffFree,calibCoeffFreeErr,calibCoeffFree_cor,calibCoeffFreeErr_cor,";
   header            += "deltaB_tr_x,deltaB_tr_xErr,deltaB_tr_y,deltaB_tr_yErr,deltaB_tr_z,deltaB_tr_zErr,";
   header            += "deltaB_pp_x,deltaB_pp_xErr,deltaB_pp_y,deltaB_pp_yErr,deltaB_pp_z,deltaB_pp_zErr,";
   header            += "dBdx_imp,dBdx_impErr,dBdy_imp,dBdy_impErr,dBdz_imp,dBdz_impErr,";
   header            += "dBdx_shim,dBdx_shimErr,dBdy_shim,dBdy_shimErr,dBdz_shim,dBdz_shimErr,";
   header            += "mis_x,mis_xErr,mis_y,mis_yErr,mis_z,mis_zErr,"; 
   header            += "misCor_x,misCor_xErr,misCor_y,misCor_yErr,misCor_z,misCor_zErr,";
   header            += "shim_x_a,shim_x_aErr,shim_x_b,shim_x_bErr,shim_x_c,shim_x_cErr,";
   header            += "shim_y_a,shim_y_aErr,shim_y_b,shim_y_bErr,shim_y_c,shim_y_cErr,";
   header            += "shim_z_a,shim_z_aErr,shim_z_b,shim_z_bErr,shim_z_c,shim_z_cErr,";
   header            += "mis_z_bar,mis_zErr_bar,misCor_z_bar,misCor_zErr_bar";

   std::ofstream outfile; 
   outfile.open(outpath);
   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      outfile << header << std::endl;
      for(int i=0;i<N;i++){
	 sprintf(outStr,"%02d,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf",
	       // data[i].calibCoeff            ,data[i].calibCoeffErr        ,
	       // data[i].calibCoeff_aba        ,data[i].calibCoeffErr_aba    ,
               i+1,
	       data[i].calibCoeff_opt        ,data[i].calibCoeffErr_opt    ,
	       // data[i].calibCoeffCor         ,data[i].calibCoeffCorErr     ,
	       // data[i].calibCoeffCor_aba     ,data[i].calibCoeffCorErr_aba ,
	       data[i].calibCoeffCor_opt     ,data[i].calibCoeffCorErr_opt ,
	       // data[i].calibCoeffFree        ,data[i].calibCoeffFreeErr        ,
	       // data[i].calibCoeffFree_aba    ,data[i].calibCoeffFreeErr_aba    ,
	       data[i].calibCoeffFree_opt    ,data[i].calibCoeffFreeErr_opt    ,
	       // data[i].calibCoeffCorFree     ,data[i].calibCoeffCorFreeErr     ,
	       // data[i].calibCoeffCorFree_aba ,data[i].calibCoeffCorFreeErr_aba ,
	       data[i].calibCoeffCorFree_opt ,data[i].calibCoeffCorFreeErr_opt ,

	       data[i].deltaB_tr_x           ,data[i].deltaB_tr_xErr,
	       data[i].deltaB_tr_y           ,data[i].deltaB_tr_yErr,
	       data[i].deltaB_tr_z           ,data[i].deltaB_tr_zErr,

	       data[i].deltaB_pp_x           ,data[i].deltaB_pp_xErr,
	       data[i].deltaB_pp_y           ,data[i].deltaB_pp_yErr,
	       data[i].deltaB_pp_z           ,data[i].deltaB_pp_zErr,

	       data[i].dBdx_imp              ,data[i].dBdx_impErr,
	       data[i].dBdy_imp              ,data[i].dBdy_impErr,
	       data[i].dBdz_imp              ,data[i].dBdz_impErr,

	       data[i].dBdx_shim             ,data[i].dBdx_shimErr,
	       data[i].dBdy_shim             ,data[i].dBdy_shimErr,
	       data[i].dBdz_shim             ,data[i].dBdz_shimErr,

	       data[i].mis_x                 ,data[i].mis_xErr,
	       data[i].mis_y                 ,data[i].mis_yErr,
	       data[i].mis_z                 ,data[i].mis_zErr,

	       data[i].misCor_x              ,data[i].misCor_xErr,
	       data[i].misCor_y              ,data[i].misCor_yErr,
	       data[i].misCor_z              ,data[i].misCor_zErr,

	       data[i].shim_x_a              ,data[i].shim_x_aErr,
	       data[i].shim_x_b              ,data[i].shim_x_bErr,
	       data[i].shim_x_c              ,data[i].shim_x_cErr,

	       data[i].shim_y_a              ,data[i].shim_y_aErr,
	       data[i].shim_y_b              ,data[i].shim_y_bErr,
	       data[i].shim_y_c              ,data[i].shim_y_cErr,

	       data[i].shim_z_a              ,data[i].shim_z_aErr,
	       data[i].shim_z_b              ,data[i].shim_z_bErr,
	       data[i].shim_z_c              ,data[i].shim_z_cErr,

	       data[i].mis_z_bar             ,data[i].mis_zErr_bar,
	       data[i].misCor_zErr_bar       ,data[i].misCor_zErr_bar); 

	 outfile << outStr << std::endl;
      }
      outfile.close();
      std::cout << "The data has been written to the file: " << outpath << std::endl;
   }
   return 0;

}
//______________________________________________________________________________
int PrintToFile_json(const char *outpath,std::vector<calib_result_t> data){

   json jData;
   int rc = GetJSONObject(data,jData);

   std::ofstream outfile;
   outfile.open(outpath); 

   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
   }else{
      outfile << std::setw(5) << jData << std::endl;
      outfile.close();
   }
   return 0;
}
//______________________________________________________________________________
int GetJSONObject(std::vector<calib_result_t> data,json &jData){

   const int N = data.size(); 

   for(int i=0;i<N;i++){
      // without free proton, without misalignment 
      jData["calibCoeff"][i]               = data[i].calibCoeff; 
      jData["calibCoeff_aba"][i]           = data[i].calibCoeff_aba; 
      jData["calibCoeff_opt"][i]           = data[i].calibCoeff_opt; 
      jData["calibCoeffErr"][i]            = data[i].calibCoeffErr; 
      jData["calibCoeffErr_aba"][i]        = data[i].calibCoeffErr_aba;
      jData["calibCoeffErr_opt"][i]        = data[i].calibCoeffErr_opt;
      // without free proton, with misalignment
      jData["calibCoeffCor"][i]            = data[i].calibCoeffCor; 
      jData["calibCoeffCor_aba"][i]        = data[i].calibCoeffCor_aba; 
      jData["calibCoeffCor_opt"][i]        = data[i].calibCoeffCor_opt; 
      jData["calibCoeffCorErr"][i]         = data[i].calibCoeffCorErr; 
      jData["calibCoeffCorErr_aba"][i]     = data[i].calibCoeffCorErr_aba;
      jData["calibCoeffCorErr_opt"][i]     = data[i].calibCoeffCorErr_opt;
      // with free proton, without misalignment 
      jData["calibCoeffFree"][i]           = data[i].calibCoeffFree; 
      jData["calibCoeffFree_aba"][i]       = data[i].calibCoeffFree_aba; 
      jData["calibCoeffFree_opt"][i]       = data[i].calibCoeffFree_opt; 
      jData["calibCoeffFreeErr"][i]        = data[i].calibCoeffFreeErr; 
      jData["calibCoeffFreeErr_aba"][i]    = data[i].calibCoeffFreeErr_aba;
      jData["calibCoeffFreeErr_opt"][i]    = data[i].calibCoeffFreeErr_opt;
      // with free proton, with misalignment 
      jData["calibCoeffCorFree"][i]        = data[i].calibCoeffCorFree; 
      jData["calibCoeffCorFree_aba"][i]    = data[i].calibCoeffCorFree_aba; 
      jData["calibCoeffCorFree_opt"][i]    = data[i].calibCoeffCorFree_opt; 
      jData["calibCoeffCorFreeErr"][i]     = data[i].calibCoeffCorFreeErr; 
      jData["calibCoeffCorFreeErr_aba"][i] = data[i].calibCoeffCorFreeErr_aba;
      jData["calibCoeffCorFreeErr_opt"][i] = data[i].calibCoeffCorFreeErr_opt;
      jData["freeErr"][i]               = data[i].freeErr;
      // free-proton related numbers
      jData["fpCor"][i]                 = data[i].fpCor;  
      jData["fpCorErr"][i]              = data[i].fpCorErr;  
      jData["ppTemp"][i]                = data[i].ppTemp;  
      jData["ppTempErr"][i]             = data[i].ppTempErr;  
      jData["trTemp"][i]                = data[i].trTemp;  
      jData["trTempErr"][i]             = data[i].trTempErr;  
      // imposed gradient data 
      jData["dBdx_imp"][i]              = data[i].dBdx_imp;      
      jData["dBdy_imp"][i]              = data[i].dBdy_imp;      
      jData["dBdz_imp"][i]              = data[i].dBdz_imp;      
      jData["dBdx_impErr"][i]           = data[i].dBdx_impErr;      
      jData["dBdy_impErr"][i]           = data[i].dBdy_impErr;      
      jData["dBdz_impErr"][i]           = data[i].dBdz_impErr;      
      // shimmed gradient data 
      jData["dBdx_shim"][i]             = data[i].dBdx_shim;      
      jData["dBdy_shim"][i]             = data[i].dBdy_shim;      
      jData["dBdz_shim"][i]             = data[i].dBdz_shim;      
      jData["dBdx_shimErr"][i]          = data[i].dBdx_shimErr;      
      jData["dBdy_shimErr"][i]          = data[i].dBdy_shimErr;      
      jData["dBdz_shimErr"][i]          = data[i].dBdz_shimErr;     
      // auxilary shimmed gradient data 
      jData["shim_x_a"][i]              = data[i].shim_x_a; 
      jData["shim_x_aErr"][i]           = data[i].shim_x_aErr; 
      jData["shim_x_b"][i]              = data[i].shim_x_b; 
      jData["shim_x_bErr"][i]           = data[i].shim_x_bErr; 
      jData["shim_x_c"][i]              = data[i].shim_x_c; 
      jData["shim_x_cErr"][i]           = data[i].shim_x_cErr; 
      jData["shim_y_a"][i]              = data[i].shim_y_a; 
      jData["shim_y_aErr"][i]           = data[i].shim_y_aErr; 
      jData["shim_y_b"][i]              = data[i].shim_y_b; 
      jData["shim_y_bErr"][i]           = data[i].shim_y_bErr; 
      jData["shim_y_c"][i]              = data[i].shim_y_c; 
      jData["shim_y_cErr"][i]           = data[i].shim_y_cErr;
      jData["shim_z_a"][i]              = data[i].shim_z_a; 
      jData["shim_z_aErr"][i]           = data[i].shim_z_aErr; 
      jData["shim_z_b"][i]              = data[i].shim_z_b; 
      jData["shim_z_bErr"][i]           = data[i].shim_z_bErr; 
      jData["shim_z_c"][i]              = data[i].shim_z_c; 
      jData["shim_z_cErr"][i]           = data[i].shim_z_cErr; 
      // delta-B data  
      jData["deltaB_pp_x"][i]           = data[i].deltaB_pp_x; 
      jData["deltaB_pp_y"][i]           = data[i].deltaB_pp_y; 
      jData["deltaB_pp_z"][i]           = data[i].deltaB_pp_z; 
      jData["deltaB_pp_xErr"][i]        = data[i].deltaB_pp_xErr; 
      jData["deltaB_pp_yErr"][i]        = data[i].deltaB_pp_yErr; 
      jData["deltaB_pp_zErr"][i]        = data[i].deltaB_pp_zErr; 
      jData["deltaB_tr_x"][i]           = data[i].deltaB_tr_x; 
      jData["deltaB_tr_y"][i]           = data[i].deltaB_tr_y; 
      jData["deltaB_tr_z"][i]           = data[i].deltaB_tr_z; 
      jData["deltaB_tr_xErr"][i]        = data[i].deltaB_tr_xErr; 
      jData["deltaB_tr_yErr"][i]        = data[i].deltaB_tr_yErr; 
      jData["deltaB_tr_zErr"][i]        = data[i].deltaB_tr_zErr; 
      // misalignment data 
      jData["dB_x"][i]                  = data[i].dB_x;
      jData["dB_y"][i]                  = data[i].dB_y;
      jData["dB_z"][i]                  = data[i].dB_z;
      jData["dB_r_tot"][i]              = data[i].dB_r_tot;
      // misalignment CORRECTION data 
      jData["misCor"][i]                = data[i].misCor;   
      jData["misCor_err"][i]            = data[i].misCor_err;
      jData["misCor_x"][i]              = data[i].misCor_x;   
      jData["misCor_xErr"][i]           = data[i].misCor_xErr;
      jData["misCor_y"][i]              = data[i].misCor_y;   
      jData["misCor_yErr"][i]           = data[i].misCor_yErr;
      jData["misCor_z"][i]              = data[i].misCor_z;   
      jData["misCor_zErr"][i]           = data[i].misCor_zErr;
      jData["misCor_z_bar"][i]          = data[i].misCor_z_bar;   
      jData["misCor_zErr_bar"][i]       = data[i].misCor_zErr_bar;
      // auxilary misalignment numbers.  These are in mm!  
      jData["mis_x"][i]                 = data[i].mis_x;  
      jData["mis_xErr"][i]              = data[i].mis_xErr;  
      jData["mis_y"][i]                 = data[i].mis_y;  
      jData["mis_yErr"][i]              = data[i].mis_yErr;  
      jData["mis_z"][i]                 = data[i].mis_z;  
      jData["mis_zErr"][i]              = data[i].mis_zErr;  
      jData["mis_z_bar"][i]             = data[i].mis_z_bar;  
      jData["mis_zErr_bar"][i]          = data[i].mis_zErr_bar;  
      // systematic uncertainty 
      jData["systErr"][i]               = data[i].systErr;  
   }
   return 0; 
}
//______________________________________________________________________________
int GetImposedGradients(const char *inpath,std::vector<imposed_gradient_t> &imp_grad){

   // gather imposed gradients; the file contains values for x,y,z 

   imposed_gradient_t x; 
   std::string sa,sg;

   std::ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      std::cout << "Cannot open the file: " << inpath << std::endl;
      return 1;
   }else{
      while( !infile.eof() ){
	 std::getline(infile,sa,',');
	 std::getline(infile,sg); 
         x.grad = std::atof( sg.c_str() ); 
	 imp_grad.push_back(x); 
      }
      infile.close();
      imp_grad.pop_back();
   }

   return 0;
}
//______________________________________________________________________________
int GetShimmedGradients(const char *prefix,int probe,std::string prodVersion,std::vector<std::string> gradName,
                        std::vector<grad_meas_t> &shim_grad){
   // load in shimmed gradients
   // we do this for N runs!

   const int NAXES = gradName.size();

   char inpath[500];
   grad_meas_t meas;
   for(int i=0;i<NAXES;i++){
      sprintf(inpath,"%s/%s-grad_pr-%02d.csv",prefix,gradName[i].c_str(),probe);
      LoadGradientData(inpath,meas);
      // push back into vector 
      shim_grad.push_back(meas);
   }

   // apply estimates if necessary
   int rc = ApplyShimmedGradEstimates(probe,prodVersion,shim_grad);

   return 0;
}
//______________________________________________________________________________
int ApplyShimmedGradEstimates(int probe,std::string prodVersion,std::vector<grad_meas_t> &data){
   // load estimates of shimmed gradient and apply to the data 
   // for use when we don't know the shimmed gradient 
   const int NAXES = 3;
   char axis[NAXES] = {'x','y','z'};

   std::string spr,sx,sy,sz;
   double ipr,ig,ix,iy,iz; 

   char inpath[200]; 
   sprintf(inpath,"./input/gradients/shim-grad-estimates_%s.csv",prodVersion.c_str());

   std::ifstream infile;
   infile.open(inpath); 
   if( infile.fail() ){
      std::cout << "Cannot open the file: " << inpath << std::endl;
      return 1;
   }else{
      while( !infile.eof() ){
	 std::getline(infile,spr,','); 
	 std::getline(infile,sx,','); 
	 std::getline(infile,sy,','); 
	 std::getline(infile,sz); 
	 ipr = std::atoi( spr.c_str() );
	 ix  = std::atof( sx.c_str() );
	 iy  = std::atof( sy.c_str() );
	 iz  = std::atof( sz.c_str() );
	 if(probe==ipr){
	    for(int i=0;i<NAXES;i++){
	       if(data[i].grad==0){
		  // bad axis, replace data 
		  if(i==0) ig = ix;
		  if(i==1) ig = iy;
		  if(i==2) ig = iz;
		  data[i].grad     = ig;
		  data[i].grad_err = fabs(ig);
		  std::cout << Form("WARNING for probe %02d: gradient estimate of %.3lf +/- %.3lf Hz applied for axis %c",
                                    probe,ig,fabs(ig),axis[i]) << std::endl;
	       }
            }
         }
      }
   }
   return 0;
}
//______________________________________________________________________________
int PrintTableOfResults(const char *outpath,int probe,double diff,double diffFree,
                        double swapErr,double mErr,double pErr,double systErr){
   // printing the final results for ALL probes to a single file.  All values in Hz 
   // diff    = PP-TRLY 
   // swapErr = point-to-point error due to the multiple swaps (RMS) 
   // mErr    = misalignment error (dB/dx)*dx 
   // pErr    = free-proton error 
   // totErr  = in-quadrature sum of all errors  
   char myStr[500]; 
   sprintf(myStr,"%02d,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf",probe,diff,diffFree,swapErr,mErr,pErr,systErr);  
   std::ofstream outfile;
   outfile.open(outpath,std::ios::app);
   if( outfile.fail() ){ 
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1; 
   }else{
      outfile << myStr << std::endl;
      outfile.close(); 
   }
   return 1; 
}
//______________________________________________________________________________
int PrintTableOfResults_db(const char *outpath,int probe,
                           double dBx_pp,double dBx_pp_err,double dBx_tr,double dBx_tr_err,
                           double dBy_pp,double dBy_pp_err,double dBy_tr,double dBy_tr_err,
                           double dBz_pp,double dBz_pp_err,double dBz_tr,double dBz_tr_err){

   // printing the final Delta-B results for ALL probes to a single file.  All values in Hz 
   char myStr[500]; 
   sprintf(myStr,"%02d,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf",
           probe,dBx_pp,dBx_pp_err,dBx_tr,dBx_tr_err,
                 dBy_pp,dBy_pp_err,dBy_tr,dBy_tr_err,
                 dBz_pp,dBz_pp_err,dBz_tr,dBz_tr_err); 
   
   std::ofstream outfile;
   outfile.open(outpath,std::ios::app);
   if( outfile.fail() ){ 
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1; 
   }else{
      outfile << myStr << std::endl;
      outfile.close(); 
   }
   return 1; 
}
//______________________________________________________________________________
int PrintTableOfResults_grad(const char *outpath,int probe,double dBdx,double dBdx_err,
                             double dBdy,double dBdy_err,double dBdz,double dBdz_err){

   // printing the final Delta-B results for ALL probes to a single file.  All values in Hz 
   char myStr[500]; 
   sprintf(myStr,"%02d,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf",probe,dBdx,dBdx_err,dBdy,dBdy_err,dBdz,dBdz_err); 
   
   std::ofstream outfile;
   outfile.open(outpath,std::ios::app);
   if( outfile.fail() ){ 
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1; 
   }else{
      outfile << myStr << std::endl;
      outfile.close(); 
   }
   return 1; 
}
//______________________________________________________________________________
int PrintTableOfResults_mis(const char *outpath,int probe,double dx,double dB_x,
                            double dy,double dB_y,double dz,double dB_z){

   // printing the final Delta-B results for ALL probes to a single file.  All values in Hz 
   char myStr[500]; 
   sprintf(myStr,"%02d,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf",probe,dx,dB_x,dy,dB_y,dz,dB_z); 
   
   std::ofstream outfile;
   outfile.open(outpath,std::ios::app);
   if( outfile.fail() ){ 
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1; 
   }else{
      outfile << myStr << std::endl;
      outfile.close(); 
   }
   return 1; 
}
