// Make tables of all results   
// TODO: Remove need to update shimmed gradient numbers with estimates    

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
#include "./src/CalibFuncs.C"

double gMarkerSize = 0.8;
double gXSize      = 0.05;  
double gYSize      = 0.06;  

int ApplyShimmedGradEstimates(int probe,std::string prodVersion,std::vector<grad_meas_t> &data); 
int PrintTableOfResults(const char *outpath,int probe,double diff,double diffFree,
                        double swapErr,double mErr,double pErr);

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

int PrintResults(const char *outpath,std::vector<calib_result_t> data,bool toROOT=true,bool toJSON=true); 

int PrintToFile_json(const char *outpath,std::vector<calib_result_t> data); 
int GetJSONObject(std::vector<calib_result_t> data,json &jData);

int CollectData(std::vector<result_prod_t> r,std::vector<result_prod_t> rFree,
                std::vector<std::vector<imposed_gradient_t>> ig,std::vector<std::vector<grad_meas_t> > mg,
                std::vector<deltab_prod_t> dB_pp,std::vector<deltab_prod_t> dB_tr,
                std::vector<misalignment_t> mis,std::vector< std::vector<misalignCor_t> > mCor,std::vector<calib_result_t> &data); 

// int MakeTables_prod(int runPeriod,bool isBlind,std::string blindLabel,std::string theDate){
int MakeTables_prod(int runPeriod,std::string theDate){

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

   char outDir[200];
   if(isBlind)  sprintf(outDir,"./output/blinded/%s",blindLabel.c_str());
   if(!isBlind) sprintf(outDir,"./output/unblinded");
   sprintf(outDir,"%s/%s",outDir,theDate.c_str()); 

   char outPath[200]; 
   sprintf(outPath,"%s/calibData_%s",outDir,theDate.c_str()); 

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

   std::vector<misalignCor_t> mc; 
   std::vector< std::vector<misalignCor_t> > mCor; 

   char inpath[500]; 

   std::vector<std::string> gradName;
   gradName.push_back("rad");  
   gradName.push_back("vert");  
   gradName.push_back("azi");  

   std::vector<int> probe; 

   bool update = false;

   double sum_sq=0;
   double sumf=0,sum_sf=0,sum=0;
   double sumf_ppb=0,sum_ppb=0;

   int k=0; 
   int probeNumber = 0; 
   const int NTRLY = 17; 
   for(int i=0;i<NTRLY;i++){
      probeNumber = i+1; 
      std::cout << Form("LOADING PROBE %02d",probeNumber) << std::endl;
      // load raw result 
      sprintf(inpath,"%s/results_final_pr-%02d.csv",outDir,probeNumber);
      rc = LoadResultsProdFinalData(inpath,result);
      if(rc!=0) continue;
      probe.push_back(probeNumber);  
      res.push_back(result); 
      // load free-proton result 
      sprintf(inpath,"%s/results_final_free-prot_pr-%02d.csv",outDir,probeNumber);
      rc = LoadResultsProdFinalData(inpath,result);
      resFree.push_back(result); 
      // load OPTIMIZED Delta-B numbers
      sprintf(inpath,"%s/delta-b-opt_pr-%02d.csv",outDir,probeNumber); 
      rc = LoadDeltaB_opt(inpath,dbPP,dbTR);
      deltaB_pp.push_back(dbPP);
      deltaB_trly.push_back(dbTR); 
      // load shimmed gradients 
      rc = GetShimmedGradients(outDir,probeNumber,prodVersion,gradName,grad); 
      shimGrad.push_back(grad);
      // clean up 
      grad.clear(); 
      // get imposed gradients 
      sprintf(inpath,"%s/imposed-gradients_pr-%02d.csv",outDir,probeNumber); 
      rc = GetImposedGradients(inpath,ig); 
      if(rc!=0) return 1;
      impGrad.push_back(ig); 
      // clean up 
      ig.clear(); 
      // load misalignment data 
      sprintf(inpath,"%s/misalignment_results_pr-%02d.csv",outDir,probeNumber);
      rc = LoadMisalignmentData(inpath,m);
      if(rc!=0) return 1;
      misalign.push_back(m);
      // load misalignment CORRECTION data 
      sprintf(inpath,"%s/misalign-cor_pr-%02d.csv",outDir,probeNumber); 
      rc = LoadMisalignmentCorData(inpath,mc); 
      mCor.push_back(mc);
      // clean up 
      mc.clear(); 
   }

   // now make tables

   // also print table of results to file 
   char outpath_res[200];
   sprintf(outpath_res,"%s/results_all-probes.csv",outDir);

   char outpath_res_free[200];
   sprintf(outpath_res_free,"%s/results_free-prot_all-probes.csv",outDir);

   char outpath_db[200];
   sprintf(outpath_db,"%s/results_db_all-probes.csv",outDir);

   char outpath_ig[200];
   sprintf(outpath_ig,"%s/results_imposed-gradients_all-probes.csv",outDir);

   char outpath_sg[200];
   sprintf(outpath_sg,"%s/results_shimmed-gradients_all-probes.csv",outDir);

   char outpath_mis[200];
   sprintf(outpath_mis,"%s/results_misalignments_all-probes.csv",outDir);

   const int NPROBES = probe.size(); 

   // PP-TRLY 
   std::cout << "===================================== PP-TRLY =====================================" << std::endl;
   for(int i=0;i<NPROBES;i++){
      sum_sq   = res[i].diffErr_aba*res[i].diffErr_aba + res[i].mErr_opt*res[i].mErr_opt;
      sum      = TMath::Sqrt(sum_sq);  
      sum_sq   = resFree[i].diffErr_aba*resFree[i].diffErr_aba + res[i].mErr_opt*res[i].mErr_opt + resFree[i].pErr*resFree[i].pErr;
      sumf     = TMath::Sqrt(sum_sq); 
      sumf_ppb = sumf/0.06179;  
      sum_sq   = resFree[i].diffErr_aba*resFree[i].diffErr_aba + res[i].mErr_opt*res[i].mErr_opt;
      sum_sf   = TMath::Sqrt(sum_sq);  
      std::cout << Form("%02d,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf",
                        probe[i],res[i].diff_aba,resFree[i].diff_aba,res[i].diffErr_aba,
                        resFree[i].diffErr_aba,res[i].mErr_opt,resFree[i].pErr,sum,sum_sf,sumf_ppb) << std::endl;
      rc = PrintTableOfResults(outpath_res,probe[i],res[i].diff_aba,resFree[i].diff_aba,res[i].diffErr_aba,res[i].mErr,resFree[i].pErr);
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
   rc = CollectData(res,resFree,impGrad,shimGrad,deltaB_pp,deltaB_trly,misalign,mCor,data); 

   rc = PrintResults(outPath,data); 

   return 0;
}
//______________________________________________________________________________
int CollectData(std::vector<result_prod_t> r,std::vector<result_prod_t> rFree,
                std::vector<std::vector<imposed_gradient_t>> ig,std::vector<std::vector<grad_meas_t> > mg,
                std::vector<deltab_prod_t> dB_pp,std::vector<deltab_prod_t> dB_tr,
                std::vector<misalignment_t> mis,std::vector< std::vector<misalignCor_t> > mCor,std::vector<calib_result_t> &data){

   // print the results to a ROOT file
   calib_result_t dataPt; 
   const int N = r.size();
   for(int i=0;i<N;i++){
      // no proton corrections 
      dataPt.calibCoeff            = r[i].diff; 
      dataPt.calibCoeffErr         = r[i].diffErr;  
      dataPt.calibCoeff_aba        = r[i].diff_aba; 
      dataPt.calibCoeffErr_aba     = r[i].diffErr_aba;
      // free proton corrections applied 
      dataPt.calibCoeffFree        = rFree[i].diff;  
      dataPt.calibCoeffFreeErr     = rFree[i].diffErr;  
      dataPt.calibCoeffFree_aba    = rFree[i].diff_aba;  
      dataPt.calibCoeffFreeErr_aba = rFree[i].diffErr_aba; 
      dataPt.freeErr               = rFree[i].pErr;  
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
      dataPt.dx                    = mis[i].dB_x_opt; 
      dataPt.dy                    = mis[i].dB_y_opt; 
      dataPt.dz                    = mis[i].dB_z_opt;
      dataPt.dr                    = TMath::Sqrt( dataPt.dx*dataPt.dx + dataPt.dy*dataPt.dy + dataPt.dz*dataPt.dz ); 
      // misalignment CORRECTION. always use the opt result.  this is in Hz  
      dataPt.misCor                = mCor[i][2].val;  
      dataPt.misCor_err            = mCor[i][2].err;  
      // fill vector 
      data.push_back(dataPt);  
   }
   return 0;
}
//______________________________________________________________________________
int PrintResults(const char *filename,std::vector<calib_result_t> data,bool toROOT,bool toJSON){
   // ROOT file parameters 
   char outpath_root[200],outpath_json[200];
   sprintf(outpath_root,"%s.root",filename); 
   sprintf(outpath_json,"%s.json",filename); 

   gm2fieldUtil::rootData_t rd;
   rd.fileName      = outpath_root;
   rd.treeName      = "CAL";
   rd.branchName    = "B";
   rd.leafStructure = calib_result_str; 

   int rc=0;
   if(toROOT) rc = gm2fieldUtil::Export::PrintToROOTFile<calib_result_t>(rd,data); 
   if(toJSON) rc = PrintToFile_json(outpath_json,data); 
   
   return rc;
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
      jData["calibCoeff"][i]            = data[i].calibCoeff; 
      jData["calibCoeff_aba"][i]        = data[i].calibCoeff_aba; 
      jData["calibCoeffErr"][i]         = data[i].calibCoeffErr; 
      jData["calibCoeffErr_aba"][i]     = data[i].calibCoeffErr_aba;
      jData["calibCoeffFree"][i]        = data[i].calibCoeffFree; 
      jData["calibCoeffFree_aba"][i]    = data[i].calibCoeffFree_aba; 
      jData["calibCoeffFreeErr"][i]     = data[i].calibCoeffFreeErr; 
      jData["calibCoeffFreeErr_aba"][i] = data[i].calibCoeffFreeErr_aba;
      jData["freeErr"][i]               = data[i].freeErr; 
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
      jData["dx"][i]                    = data[i].dx;
      jData["dy"][i]                    = data[i].dy;
      jData["dz"][i]                    = data[i].dz;
      jData["dr"][i]                    = data[i].dr;
      // misalignment CORRECTION data 
      jData["misCor"][i]                = data[i].misCor;   
      jData["misCor_err"][i]            = data[i].misCor_err;   
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
                        double swapErr,double mErr,double pErr){
   // printing the final results for ALL probes to a single file.  All values in Hz 
   // diff    = PP-TRLY 
   // swapErr = point-to-point error due to the multiple swaps (RMS) 
   // mErr    = misalignment error (dB/dx)*dx 
   // pErr    = free-proton error 
   // totErr  = in-quadrature sum of all errors  
   char myStr[500]; 
   sprintf(myStr,"%02d,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf",probe,diff,diffFree,swapErr,mErr,pErr);  
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
