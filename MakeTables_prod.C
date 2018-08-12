// Make tables of all results    

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
#include "./include/grad_meas.h"
#include "./include/deltab.h"
#include "./include/results.h"
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

int GetImposedGradients(const char *inpath,std::vector<imposed_gradient_t> &imp_grad); 
int GetShimmedGradients(const char *prefix,int probe,std::vector<std::string> gradName,
                        std::vector<grad_meas_t> &shim_grad);

int MakeTables_prod(bool isBlind){

   int rc=0;
   int method = gm2fieldUtil::Constants::kPhaseDerivative;

   date_t theDate;
   GetDate(theDate);

   char outDir[200];
   sprintf(outDir,"./output"); 
   if(isBlind)  sprintf(outDir,"%s/blinded"  ,outDir);
   if(!isBlind) sprintf(outDir,"%s/unblinded",outDir);
   sprintf(outDir,"%s/%s",outDir,theDate.getDateString().c_str()); 

   result_prod_t result;
   std::vector<result_prod_t> res,resFree;

   std::vector<imposed_gradient_t> ig;   
   std::vector< std::vector<imposed_gradient_t> > impGrad;   

   // deltab_prod_t DB; 
   std::vector<deltab_t> db; 
   std::vector< std::vector<deltab_t> > deltaB_pp,deltaB_trly; 

   std::vector<grad_meas_t> grad; 
   std::vector< std::vector<grad_meas_t> > shimGrad;

   misalignment_t m; 
   std::vector<misalignment_t> misalign;

   char inpath[500]; 

   std::vector<std::string> gradName;
   gradName.push_back("rad");  
   gradName.push_back("vert");  
   gradName.push_back("azi");  

   std::vector<int> probe; 

   bool update = false;

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
      // load PP Delta-B  
      for(int j=0;j<3;j++){
	 sprintf(inpath,"%s/dB-pp_final-location_%s-grad_pr-%02d.csv",outDir,gradName[j].c_str(),probeNumber); 
	 rc = LoadDeltaBData(inpath,db);
	 if(db[j].dB_fxpr==0){
	    std::cout << Form("--> No PP ABA data for probe %02d, axis %d.  Using raw data",probeNumber,j) << std::endl; 
	    db[j].dB_fxpr     = db[j].dB;
	    db[j].dB_fxpr_err = db[j].dB_err;
	    update = true;
         }
      }
      // clean up the PP delta B container, and store in new vector 
      deltaB_pp.push_back(db);
      db.clear();
      // load TRLY Delta-B values
      // rad and vert already computed offline
      sprintf(inpath,"./input/delta-b/trly_xy_07-18.csv");
      LoadDeltaBData_trlyXY(inpath,probeNumber,db);
      // z axis 
      sprintf(inpath,"%s/dB-trly_final-location_azi-grad_pr-%02d.csv",outDir,probeNumber);
      LoadDeltaBData(inpath,db);  
      for(int j=0;j<3;j++){
	 if(db[j].dB_fxpr==0){
	    std::cout << Form("--> No TRLY ABA data for probe %02d, axis %d.  Using raw data",probeNumber,j) << std::endl; 
	    std::cout << db[j].dB << " " << db[j].dB_err << std::endl;
	    db[j].dB_fxpr     = db[j].dB;
	    db[j].dB_fxpr_err = db[j].dB_err;
	    update = true;
         }
      }
      // clean up the TRLY delta B container, and store in new vector  
      deltaB_trly.push_back(db);
      db.clear(); 
      // load shimmed gradients 
      rc = GetShimmedGradients(outDir,probeNumber,gradName,grad); 
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
      std::cout << misalign.size() << std::endl;
      if(update){
	 // need to recompute misalignment numbers
	 misalign[k].dx_aba   = (deltaB_pp[k][0].dB_fxpr-deltaB_trly[k][0].dB_fxpr)/impGrad[k][0].grad; 
	 misalign[k].dy_aba   = (deltaB_pp[k][1].dB_fxpr-deltaB_trly[k][1].dB_fxpr)/impGrad[k][1].grad; 
	 misalign[k].dz_aba   = (deltaB_pp[k][2].dB_fxpr-deltaB_trly[k][2].dB_fxpr)/impGrad[k][2].grad; 
	 misalign[k].dB_x_aba = misalign[k].dx_aba*shimGrad[k][0].grad; 
	 misalign[k].dB_y_aba = misalign[k].dy_aba*shimGrad[k][1].grad; 
	 misalign[k].dB_z_aba = misalign[k].dz_aba*shimGrad[k][2].grad; 
	 std::cout << "--> Misalignments recalculated." << std::endl;
      }
      // reset 
      update = false; 
      k++;
   }

   // now make tables

   const int NPROBES = probe.size(); 

   double sum=0,sum_sq=0;

   // PP-TRLY 
   std::cout << "===================================== PP-TRLY =====================================" << std::endl;
   for(int i=0;i<NPROBES;i++){
      sum_sq = res[i].diffErr_aba*res[i].diffErr_aba + res[i].mErr_aba*res[i].mErr_aba + resFree[i].pErr*resFree[i].pErr;
      sum    = TMath::Sqrt(sum_sq)/0.06179;  
      std::cout << Form("%02d,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf",
                        probe[i],res[i].diff_aba,resFree[i].diff_aba,res[i].diffErr_aba,res[i].mErr_aba,resFree[i].pErr,sum) << std::endl;
   }

   std::cout << "===================================== DELTA B =====================================" << std::endl;
   for(int i=0;i<NPROBES;i++){
      std::cout << Form("%02d,%.3lf ± %.3lf,%.3lf ± %.3lf,%.3lf ± %.3lf,%.3lf ± %.3lf,%.3lf ± %.3lf,%.3lf ± %.3lf",probe[i],
                        deltaB_pp[i][0].dB_fxpr,deltaB_pp[i][0].dB_fxpr_err,deltaB_trly[i][0].dB_fxpr,deltaB_trly[i][0].dB_fxpr_err,
                        deltaB_pp[i][1].dB_fxpr,deltaB_pp[i][1].dB_fxpr_err,deltaB_trly[i][1].dB_fxpr,deltaB_trly[i][1].dB_fxpr_err,
                        deltaB_pp[i][2].dB_fxpr,deltaB_pp[i][2].dB_fxpr_err,deltaB_trly[i][2].dB_fxpr,deltaB_trly[i][2].dB_fxpr_err) << std::endl;
   }

   std::cout << "===================================== MISALIGNMENTS =====================================" << std::endl;
   for(int i=0;i<NPROBES;i++){
      sum_sq = misalign[i].dB_x_aba*misalign[i].dB_x_aba + misalign[i].dB_y_aba*misalign[i].dB_y_aba + misalign[i].dB_z_aba*misalign[i].dB_z_aba;
      sum    = TMath::Sqrt(sum_sq); 
      std::cout << Form("%02d,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf",probe[i],
                        misalign[i].dx_aba,misalign[i].dB_x_aba,
                        misalign[i].dy_aba,misalign[i].dB_y_aba,
                        misalign[i].dz_aba,misalign[i].dB_z_aba,
                        sum) << std::endl;
   }
   
   std::cout << "===================================== SHIMMED GRADIENTS =====================================" << std::endl;
   for(int i=0;i<NPROBES;i++){
      std::cout << Form("%02d,%.3lf ± %.3lf,%.3lf ± %.3lf,%.3lf ± %.3lf",
                        probe[i],shimGrad[i][0].grad,shimGrad[i][0].grad_err,
                        shimGrad[i][1].grad,shimGrad[i][1].grad_err,
                        shimGrad[i][2].grad,shimGrad[i][2].grad_err) << std::endl;
   }

   std::cout << "===================================== IMPOSED GRADIENTS =====================================" << std::endl;
   for(int i=0;i<NPROBES;i++){
      std::cout << Form("%02d,%.3lf,%.3lf,%.3lf",
                        probe[i],impGrad[i][0].grad,impGrad[i][1].grad,impGrad[i][2].grad) << std::endl;
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
int GetShimmedGradients(const char *prefix,int probe,std::vector<std::string> gradName,
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
   return 0;
}
