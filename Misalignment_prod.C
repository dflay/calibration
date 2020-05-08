// Compare PP and TRLY Delta-B measurements, imposed and shimmed gradients to obtain 
// misalignment uncertainties   

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
// #include "./src/CalibFuncs.C"
#include "./src/TRLYFuncs.C"

double gMarkerSize = 0.8;
double gXSize      = 0.05;  
double gYSize      = 0.06;  

double GetQ0Error(double DB,double DB_err,double dBidq,double dBidq_err); 
double GetGradientError(double key,std::vector<double> pos,std::vector<double> err); 
double GetMisalignmentError(double dBsdx,double dBsdx_err,double dB,double dB_err,
                            double dBidx,double dBidx_err);
 
int GetShimmedGradients(const char *prefix,int runPeriod,int probe,
                        std::vector<std::string> gradName,std::vector<grad_meas_t> &shim_grad_avg);

int GetImposedGradients(const char *indir,int probeNumber,int fitDim,std::vector<double> coord,
                        std::vector<double> &impG,std::vector<double> &impGerr); 

int GetDeltaB_optimal(deltab_t data,double &DB,double &DB_err); 
int PrintShimmedGradToFile(const char *outpath,std::vector<grad_meas_t> data); 

int Misalignment_prod(std::string configFile){

   char msg[200]; 
   
   int rc=0;
   int method = gm2fieldUtil::Constants::kPhaseDerivative;

   InputManager *inputMgr = new InputManager();
   inputMgr->Load(configFile);
   inputMgr->Print();

   std::string fitFunc    = "pol2";  // FIXME: should add this to the input manager... 
   std::string blindLabel = inputMgr->GetBlindLabel(); 
   std::string prodTag    = inputMgr->GetProductionTag();
 
   bool isBlind           = inputMgr->IsBlind();
   bool useOscCor         = inputMgr->GetOscCorStatus(); 
   bool useShimGradAlt    = inputMgr->GetShimGradAltStatus(); 
   int probe              = inputMgr->GetTrolleyProbe();
   int runPeriod          = inputMgr->GetRunPeriod(); 
   int impGradFitDim      = inputMgr->GetImpGradFitDimension(); 
   // systematics 
   bool isSyst            = inputMgr->GetSystStatus();
   int systDirNum         = inputMgr->GetSystDirNum(); 

   // date_t theDate;
   // GetDate(theDate);
   std::string theDate = inputMgr->GetAnaDate();

   std::string plotDir = GetPath("plots" ,isBlind,blindLabel,theDate,isSyst,systDirNum);
   std::string outDir  = GetPath("output",isBlind,blindLabel,theDate,isSyst,systDirNum);

   char outPath[500],outPath_fxpr[500],outPath_ig[500],outPath_dB[500]; 
   sprintf(outPath     ,"%s/misalignment_results_pr-%02d.csv"     ,outDir.c_str(),probe);
   sprintf(outPath_fxpr,"%s/misalignment_results_fxpr_pr-%02d.csv",outDir.c_str(),probe);
   sprintf(outPath_ig  ,"%s/imposed-gradients_pr-%02d.csv"        ,outDir.c_str(),probe);
   sprintf(outPath_dB  ,"%s/delta-b-opt_pr-%02d.csv"              ,outDir.c_str(),probe); // optimized Delta-B values 
   std::string outpath    = outPath; 

   std::vector<std::string> gradName;
   gradName.push_back("rad");
   gradName.push_back("vert");
   gradName.push_back("azi");

   std::vector<std::string> axis; 
   axis.push_back("r"); 
   axis.push_back("y"); 
   axis.push_back("z"); 

   const int NG = gradName.size(); 

   // load PP Delta-B values
   char inpath[500];  
   std::vector<deltab_t> pp,trly;
   for(int i=0;i<NG;i++){
      sprintf(inpath,"%s/dB-pp_final-location_%s-grad_pr-%02d.csv",outDir.c_str(),gradName[i].c_str(),probe);
      LoadDeltaBData(inpath,pp);
   }

   // now get TRLY Delta-B 
   for(int i=0;i<NG;i++){ 
      sprintf(inpath,"%s/dB-trly_final-location_%s-grad_pr-%02d.csv",outDir.c_str(),gradName[i].c_str(),probe); 
      LoadDeltaBData(inpath,trly);
   }

   // std::cout << "TRLY DB VALS" << std::endl;
   // const int NN = trly.size();
   // for(int i=0;i<NN;i++){
   //    std::cout << trly[i].dB << std::endl;
   // }

   // gather trolley probe coordinates
   trolleyProbePosition_t trlyProbePos;
   rc = GetTrolleyProbePositions(trlyProbePos); 

   double rad_coord  = trlyProbePos.r[probe-1]; 
   double vert_coord = trlyProbePos.y[probe-1]; 
   double azi_coord  = trlyProbePos.phi[probe-1]; 

   // load shimmed gradients.  We took NRUNS, so we'll average over them
   // FIXME: 5/14/19: I don't see what I meant by averaging.  
   std::vector<grad_meas_t> shim_grad;
   rc = GetShimmedGradients(outDir.c_str(),runPeriod,probe,gradName,shim_grad);

   CheckShimGrad(probe,shim_grad); 

   // load imposed gradients
   std::vector<double> coord;
   coord.push_back(rad_coord);  
   coord.push_back(vert_coord);  
   coord.push_back(azi_coord);  
   std::vector<double> imposedGrad,imposedGrad_err;
   rc = GetImposedGradients(outDir.c_str(),probe,impGradFitDim,coord,imposedGrad,imposedGrad_err); 
   
   // check to make sure gradients are non-zero 
   const int NI = imposedGrad.size();
   for(int i=0;i<NI;i++){
      if(imposedGrad[i]==0){
	 imposedGrad[i] = 1E+3;  // make something crazy big
	 sprintf(msg,"[Misalignment_prod]: Invalid imposed gradient for axis %d!  Using 1000 Hz.",i);
         Logger::PrintMessage(Logger::kERR,"default",msg,'a');
      } 
   }

   // load barcode correction if necessary 
   char inpath_bc[200]; 
   sprintf(inpath_bc,"./input/barcode/run-%d/barcode-cor.csv",runPeriod); 
   std::vector<double> dz_bc,dz_bc_err;
   CSVManager *csvMgr = new CSVManager(); 
   csvMgr->ReadFile(inpath_bc,true);
   csvMgr->GetColumn_byName("dz_bar"    ,dz_bc); 
   csvMgr->GetColumn_byName("dz_bar_err",dz_bc_err);
   delete csvMgr;  

   double DZ_BC     = dz_bc[probe-1]; 
   double DZ_BC_ERR = dz_bc_err[probe-1]; 

   // do we flip the sign on the dB difference? 
   double q0_sf = 1.; 
   if(useShimGradAlt) q0_sf = -1.; 

   // compute misalignment error in mm and in Hz 
   double cor=0,cor_err_sq=0;                    // for eventually determining the size of a correction term (raw) 
   double cor_aba=0,cor_err_aba_sq=0;            // for eventually determining the size of a correction term (ABA) 
   double cor_opt=0,cor_err_opt_sq=0;            // for eventually determining the size of a correction term (opt) 
   double cor_bar_opt=0,cor_err_bar_opt_sq=0;    // for eventually determining the size of a correction term (opt) (barcode on z) 
   double arg=0,arg_aba=0,arg_opt=0,ddB=0,ddB_err=0,merr=0,merr_aba=0,merr_opt=0;
   double ddB_opt=0,ddB_opt_err=0,pp_dB_opt=0,pp_dB_opt_err=0;
   double trly_dB_opt=0,trly_dB_opt_err=0;
   double shim_grad_opt=0,shim_grad_alt=0,shim_grad_alt_err=0;
   double dqe=0,dqe_aba=0,dqe_opt=0,dqe_bar_opt=0; 
   double dq_arg=0,dq_arg_aba=0,dq_arg_opt=0;                         // misalignment in mm (optimized)  
   double dq_arg_bar_opt=0,arg_bar_opt=0;                             // misalignment in mm (optimized)  
   double merr_bar_opt=0;
   const int NAXES = gradName.size(); 
   std::vector<double> PP_DB,PP_DB_ERR;                         // PP Delta-B (optimized)    
   std::vector<double> TRLY_DB,TRLY_DB_ERR;                     // TRLY Delta-B (optimized)
   std::vector<double> dq,dq_aba,dq_opt,dq_bar_opt;             // misalignment in mm 
   std::vector<double> dqERR,dqERR_aba,dqERR_opt,dqERR_bar_opt; // error on misalignment in mm 
   std::vector<double> mERR,mERR_aba,mERR_opt;
   std::vector<double> vc,vce;                       // misalignment correction as a function of axis 
   std::vector<double> vcb,vcbe;                     // misalignment correction as a function of axis (w/ barcode for z) 
   for(int i=0;i<NAXES;i++){
      // raw 
      ddB     = pp[i].dB - trly[i].dB;   // difference in Delta-B
      ddB_err = TMath::Sqrt(pp[i].dB_err*pp[i].dB_err + trly[i].dB_err*trly[i].dB_err); 
      // misalignment in mm
      dq_arg  = q0_sf*ddB/imposedGrad[i];      // dividing by imposedGrad gives error in mm 
      dq.push_back(dq_arg);
      dqe     = GetQ0Error(ddB,ddB_err,imposedGrad[i],imposedGrad_err[i]); 
      dqERR.push_back(dqe);
      // correction in Hz
      arg     = shim_grad[i].grad*dq[i]; 
      merr    = GetMisalignmentError(shim_grad[i].grad,shim_grad[i].grad_err,ddB,ddB_err,imposedGrad[i],imposedGrad_err[i]);
      mERR.push_back(merr); 
      // ABA drift corrected 
      ddB        = pp[i].dB_fxpr - trly[i].dB_fxpr; 
      ddB_err    = TMath::Sqrt(pp[i].dB_fxpr_err*pp[i].dB_fxpr_err + trly[i].dB_fxpr_err*trly[i].dB_fxpr_err); 
      // misalignment in mm
      dq_arg_aba = q0_sf*ddB/imposedGrad[i]; 
      dq_aba.push_back(dq_arg_aba);
      dqe_aba    = GetQ0Error(ddB,ddB_err,imposedGrad[i],imposedGrad_err[i]); 
      dqERR_aba.push_back(dqe_aba);
      // correction in Hz
      arg_aba    = shim_grad[i].grad_fxpr*dq_aba[i]; 
      merr_aba   = GetMisalignmentError(shim_grad[i].grad_fxpr,shim_grad[i].grad_fxpr_err,ddB,ddB_err,imposedGrad[i],imposedGrad_err[i]);
      mERR_aba.push_back(merr_aba); 
      // OPTIMIZED: Use ABA where possible, otherwise use raw
      // find Delta-B values
      GetDeltaB_optimal(pp[i]  ,pp_dB_opt  ,pp_dB_opt_err  );
      GetDeltaB_optimal(trly[i],trly_dB_opt,trly_dB_opt_err);
      // difference in Delta-B 
      ddB_opt     = pp_dB_opt - trly_dB_opt; 
      ddB_opt_err = TMath::Sqrt(pp_dB_opt_err*pp_dB_opt_err + trly_dB_opt_err*trly_dB_opt_err);
      PP_DB.push_back(pp_dB_opt);  
      TRLY_DB.push_back(trly_dB_opt);  
      PP_DB_ERR.push_back(pp_dB_opt_err);  
      TRLY_DB_ERR.push_back(trly_dB_opt_err);  
      // misalignment in mm   
      dq_arg_opt  = q0_sf*ddB_opt/imposedGrad[i];
      dq_opt.push_back(dq_arg_opt);
      dqe_opt     = GetQ0Error(ddB_opt,ddB_opt_err,imposedGrad[i],imposedGrad_err[i]); 
      dqERR_opt.push_back(dqe_opt);
      // correction term   
      arg_opt     = shim_grad[i].grad_fxpr*dq_arg_opt;  
      merr_opt    = GetMisalignmentError(shim_grad[i].grad_fxpr,shim_grad[i].grad_fxpr_err,ddB_opt,ddB_opt_err,imposedGrad[i],imposedGrad_err[i]);  
      // now do a barcode correction for z; only for opt since that's what matters 
      // first make a copy of necessary values 
      dq_arg_bar_opt = dq_arg_opt; 
      dqe_bar_opt    = dqe_opt; 
      arg_bar_opt    = arg_opt; 
      merr_bar_opt   = merr_opt;
      if(i==2){
         // apply a barcode correction and add add'l uncertainty in quadrature 
	 std::cout << Form("[Misalignment_prod]: Applying barcode correction = %.3lf +/- %.3lf mm",DZ_BC,DZ_BC_ERR) << std::endl;
         dq_arg_bar_opt += DZ_BC; 
         dqe_bar_opt     = TMath::Sqrt( dqe_bar_opt*dqe_bar_opt + DZ_BC_ERR*DZ_BC_ERR );  
	 shim_grad_opt   = shim_grad[i].grad_fxpr; 
	 arg_bar_opt     = shim_grad_opt*dq_arg_bar_opt; 
         merr_bar_opt    = TMath::Sqrt( merr_bar_opt*merr_bar_opt + shim_grad_opt*shim_grad_opt*DZ_BC_ERR*DZ_BC_ERR );  
      }
      dq_bar_opt.push_back(dq_arg_bar_opt);
      dqERR_bar_opt.push_back(dqe_bar_opt);
      // now compute a correction size if we're to correct for this misalignment 
      cor                 += arg; 
      cor_err_sq          += merr*merr;  
      cor_aba             += arg_aba; 
      cor_err_aba_sq      += merr_aba*merr_aba;  
      cor_opt             += arg_opt; 
      cor_err_opt_sq      += merr_opt*merr_opt; 
      cor_bar_opt         += arg_bar_opt; 
      cor_err_bar_opt_sq  += merr_bar_opt*merr_bar_opt; 
      // store each axis result
      vc.push_back(arg_opt);
      vce.push_back(merr_opt);
      // barcode applied for z  
      vcb.push_back(arg_bar_opt);
      vcbe.push_back(merr_bar_opt);  
   } 

   // error on correction term 
   double cor_err         = TMath::Sqrt(cor_err_sq); 
   double cor_aba_err     = TMath::Sqrt(cor_err_aba_sq); 
   double cor_opt_err     = TMath::Sqrt(cor_err_opt_sq); 
   double cor_bar_opt_err = TMath::Sqrt(cor_err_bar_opt_sq); 

   std::cout << Form("------------------------- RESULTS FOR TRLY PROBE %02d -------------------------",probe) << std::endl;
   for(int i=0;i<NAXES;i++){
      std::cout << Form("axis: %s",axis[i].c_str()) << std::endl;
      std::cout << Form("PP Delta-B") << std::endl;
      std::cout << Form("raw = %.3lf +/- %.3lf Hz, ABA = %.3lf +/- %.3lf Hz",
                        pp[i].dB,pp[i].dB_err,pp[i].dB_fxpr,pp[i].dB_fxpr_err) << std::endl;
      std::cout << Form("TRLY Delta-B") << std::endl;
      std::cout << Form("raw = %.3lf +/- %.3lf Hz, ABA = %.3lf +/- %.3lf Hz",
                        trly[i].dB,trly[i].dB_err,trly[i].dB_fxpr,trly[i].dB_fxpr_err) << std::endl;
      std::cout << Form("GRADIENTS") << std::endl; 
      std::cout << Form("imposed gradient = %.3lf +/- %.3lf Hz/mm",imposedGrad[i]   ,imposedGrad_err[i])    << std::endl; 
      std::cout << Form("shimmed gradient = %.3lf +/- %.3lf Hz/mm",shim_grad[i].grad,shim_grad[i].grad_err) << std::endl; 
      std::cout << Form("MISALIGNMENT") << std::endl;
      std::cout << Form("raw: %.3lf mm (cor = %.3lf +/- %.3lf Hz)",dq[i]    ,vc[i],vce[i]) << std::endl;
      std::cout << Form("ABA: %.3lf mm (cor = %.3lf +/- %.3lf Hz)",dq_aba[i],vc[i],vce[i]) << std::endl;
      std::cout << Form("opt: %.3lf mm (cor = %.3lf +/- %.3lf Hz)",dq_opt[i],vc[i],vce[i]) << std::endl;
      std::cout << Form("MISALIGNMENT (BARCODE APPLIED FOR Z)") << std::endl;
      std::cout << Form("opt: %.3lf mm (cor = %.3lf +/- %.3lf Hz)",dq_bar_opt[i] ,vcb[i],vcbe[i]) << std::endl;
      // std::cout << Form("raw: %.3lf mm (%.3lf +/- %.3lf Hz)",dq[i]     ,ERR[i]     ,mERR[i]     ) << std::endl;
      // std::cout << Form("ABA: %.3lf mm (%.3lf +/- %.3lf Hz)",dq_fxpr[i],ERR_fxpr[i],mERR_fxpr[i]) << std::endl;
      // std::cout << Form("opt: %.3lf mm (%.3lf +/- %.3lf Hz)",dq_opt[i] ,ERR_opt[i] ,mERR_opt[i] ) << std::endl;
      std::cout << "--------------" << std::endl;
   }
   std::cout << "CORRECTION TERM (sum over all axes)" << std::endl;
   std::cout << Form("raw     = %.3lf +/- %.3lf Hz",cor        ,cor_err)         << std::endl;
   std::cout << Form("ABA     = %.3lf +/- %.3lf Hz",cor_aba    ,cor_aba_err)     << std::endl;
   std::cout << Form("opt     = %.3lf +/- %.3lf Hz",cor_opt    ,cor_opt_err)     << std::endl;
   std::cout << Form("opt(bc) = %.3lf +/- %.3lf Hz",cor_bar_opt,cor_bar_opt_err) << std::endl;
   std::cout << "-------------------------------------------------------------"  << std::endl;

   std::vector<double> COR,COR_ERR;
   COR.push_back(cor);
   COR.push_back(cor_aba);
   COR.push_back(cor_opt);
   COR.push_back(cor_bar_opt);
   COR_ERR.push_back(cor_err);
   COR_ERR.push_back(cor_aba_err);
   COR_ERR.push_back(cor_opt_err);
   COR_ERR.push_back(cor_bar_opt_err);
   std::vector<std::string> label; 
   label.push_back("raw"); 
   label.push_back("ABA"); 
   label.push_back("opt"); 
   label.push_back("optbc"); 

   rc = PrintToFile(outPath   ,axis,dq,dqERR,dq_aba,dqERR_aba,dq_opt,dqERR_opt,dq_bar_opt,dqERR_bar_opt);
   rc = PrintToFile(outPath_ig,axis,imposedGrad); 
   rc = PrintToFile(outPath_dB,axis,PP_DB,PP_DB_ERR,TRLY_DB,TRLY_DB_ERR); 

   // print what we've found since these are vetted to be the right combination of raw and ABA results  
   char outpath_cor[200];
   sprintf(outpath_cor,"%s/misalign-cor_pr-%02d.csv",outDir.c_str(),probe); 
   rc = PrintToFile(outpath_cor,label,COR,COR_ERR);
   // individual axis results [opt] -- includes barcode version!
   sprintf(outpath_cor,"%s/misalign-cor-opt_by-axis_pr-%02d.csv",outDir.c_str(),probe); 
   rc = PrintToFile_4dbl(outpath_cor,vc,vce,vcb,vcbe);
   // updates shimmed field gradients (replaced missing values with worst gradient from other axes)  
   char outpath_sg[200];
   sprintf(outpath_sg,"%s/shim-grad_opt_pr-%02d.csv",outDir.c_str(),probe); 
   rc = PrintShimmedGradToFile(outpath_sg,shim_grad); 

   return 0;
}
//______________________________________________________________________________
double GetQ0Error(double DB,double DB_err,double dBidq,double dBidq_err){
   // quadrature sum of uncertainties 
   double q0=0,T1=0,T2=0;
   if(DB!=0){
      T1 = DB_err/DB;
   }
   if(dBidq!=0){ 
      T2 = dBidq_err/dBidq;
      q0 = DB/dBidq;
   } 
   double err = q0*TMath::Sqrt( TMath::Power(T1,2.) + TMath::Power(T2,2.) ); 
   return err;
}
//______________________________________________________________________________
double GetMisalignmentError(double dBsdx,double dBsdx_err,double dB,double dB_err,
                            double dBidx,double dBidx_err){
   // dBsdx = shimmed field gradient; dBidx = imposed gradient
   // still YET a new way!
   // this is a propagation of errors, writing F = A*(C/D) and taking derivatives.
   bool noShimGrad=false;
   double T1=0,T2=0,T3=0; 
   if(dB!=0){
      T1 = dB_err/dB;
   } 
   if(dBidx!=0){
      T2 = dBidx_err/dBidx;
   } 
   if(dBsdx!=0){
      T3 = dBsdx_err/dBsdx; 
   }else{
      noShimGrad = true;
   }
   // compute the uncertainty 
   double F=0,sum_sq=0,DF=0,q0=0;
   if(noShimGrad){
      // no shimmed gradient!  Need a reasonable estimate 
      q0     = dB/dBidx;
      F      = 0; 
      DF     = TMath::Sqrt( q0*q0*dBsdx_err*dBsdx_err + q0*q0*dBidx_err*dBidx_err + dB_err*dB_err ); // all quantities in Hz 
   }else{
      // shimmed gradient exists; use usual quadrature sum of fractional uncertainties 
      F      = dBsdx*(dB/dBidx);  
      sum_sq = TMath::Power(T1,2.) + TMath::Power(T2,2.) + TMath::Power(T3,2.); 
      DF     = TMath::Abs(F)*TMath::Sqrt(sum_sq);
   } 
   return DF; 
   // // new way
   // // dx err 
   // double dx     = dB/dBidx; 
   // double sum_sq = TMath::Power(dB_err/dB,2.) + TMath::Power(dBidx_err/dBidx,2.); 
   // double dx_err = dx*TMath::Sqrt(sum_sq); 
   // // now put it all together
   // double F      = dBsdx*dx;  
   // sum_sq        = TMath::Power(dx_err/dx,2.) + TMath::Power(dBsdx_err/dBsdx,2.); 
   // double DF     = F*TMath::Sqrt(sum_sq);
   // return DF; 
   // // old way 
   // double dx = dB/dBidx; 
   // double T1 = dBsdx*dx; 
   // double a1 = dB_err/dBidx; 
   // double a2 = (-1.*dB*dBidx_err)/(dBidx*dBidx);
   // double T2 = dBsdx*(a1+a2); 
   // double arg_sq = T1*T1 + T2*T2; 
   // double err = TMath::Sqrt(arg_sq); 
   // return err;  
}
//______________________________________________________________________________
int GetShimmedGradients(const char *prefix,int runPeriod,int probe,
                        std::vector<std::string> gradName,std::vector<grad_meas_t> &shim_grad_avg){
   // load in shimmed gradients
   // we do this for N runs!
 
   const int NAXES = gradName.size(); 

   std::vector<double> grad,grad_aba; 
   double mean=0,stdev=0;
   double mean_aba=0,stdev_aba=0;
 
   char inpath[500]; 
   grad_meas_t meas; 
   std::vector<grad_meas_t> shim_grad; 
   for(int i=0;i<NAXES;i++){
      sprintf(inpath,"%s/%s-grad_pr-%02d.csv",prefix,gradName[i].c_str(),probe);
      // // determine which file to look at 
      // if( runPeriod==1 && i==2 && (probe==9||probe==11) ){
      //    // need an estimate for probes 9 and 11 for z axis
      //    std::cout << Form("[Misalignment::GetShimmedGradients]: WARNING! Will use ESTIMATE for probe %d, axis %d",probe,i) << std::endl; 
      //    sprintf(inpath,"%s/%s-grad_est_pr-%02d.csv",prefix,gradName[i].c_str(),probe);
      // }else{
      //    sprintf(inpath,"%s/%s-grad_pr-%02d.csv",prefix,gradName[i].c_str(),probe);
      // }
      LoadGradientData(inpath,meas);
      // push back into vector 
      shim_grad_avg.push_back(meas); 
   }
   return 0;
}
//______________________________________________________________________________
double GetGradientError(double key,std::vector<double> x,std::vector<double> err){
   // find the error on the gradient by searching the position vector x for the key value 
   // and using the corresponding unceratinty from the err vector; take average if we 
   // find an in-between coordinate that's bounded by the indices found by the binary search  
   int lo=0,hi=0;
   gm2fieldUtil::Algorithm::BinarySearch<double>(x,key,lo,hi);

   const int N = x.size(); 
   if(lo<0) lo = 0;
   if(hi>N) hi = N-1;  

   double ERR=0; 
   if( x[lo]==key ){
      ERR = err[lo]; 
   }else if( x[hi]==key ){
      ERR = err[hi];
   }else{
      ERR = 0.5*(err[lo] + err[hi]);
   }
   return ERR; 
}
//______________________________________________________________________________
int GetDeltaB_optimal(deltab_t data,double &DB,double &DB_err){
   // determine which values, raw or ABA, to use 
   if(data.dB_fxpr!=0){
      DB     = data.dB_fxpr;
      DB_err = data.dB_fxpr_err;
   }else{
      DB     = data.dB;
      DB_err = data.dB_err;
   }
   return 0;
}
//______________________________________________________________________________
int PrintShimmedGradToFile(const char *outpath,std::vector<grad_meas_t> data){
   // print gradient results to file 
   const int N = data.size();
   char axis[3] = {'r','y','z'};
   char outStr[200]; 

   std::ofstream outfile;
   outfile.open(outpath);
   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      for(int i=0;i<N;i++){
	 sprintf(outStr,"%c,%.3lf,%.3lf,%.3lf,%.3lf",axis[i],data[i].grad,data[i].grad_err,data[i].grad_fxpr,data[i].grad_fxpr_err);
         outfile << outStr << std::endl;  
      }
      outfile.close();
   }
   return 0;
}
//______________________________________________________________________________
int GetImposedGradients(const char *indir,int probeNumber,int fitDim,std::vector<double> coord,
                        std::vector<double> &impG,std::vector<double> &impGerr){

   std::string path;
   char inpath_igrad[500];
   std::vector<double> xPar,xParErr;
   std::vector<double> yPar,yParErr;
   std::vector<double> RAD,VER,radGrad_err,vertGrad_err;  
   std::vector<double> igx,igxErr;
   std::vector<double> igy,igyErr;

   double rad_grad,rad_grad_err;
   double vert_grad,vert_grad_err;

   TF1 *fitRadGrad  = NULL;
   TF1 *fitVertGrad = NULL;

   int rc=0;
   if(fitDim==1){
      // radial 
      sprintf(inpath_igrad,"%s/imposed-grad-x_fit-pars.csv",indir);
      rc = LoadFitPars(inpath_igrad,xPar,xParErr);
      // vertical 
      sprintf(inpath_igrad,"%s/imposed-grad-y_fit-pars.csv",indir);
      rc = LoadFitPars(inpath_igrad,yPar,yParErr);
      // load fit errors for rad and vert gradients
      // radial
      sprintf(inpath_igrad,"%s/imposed-grad-x_fit-err.csv",indir);
      path = inpath_igrad;  
      rc = gm2fieldUtil::Import::ImportData2<double,double>(path,"csv",VER,radGrad_err);
      // vertical 
      sprintf(inpath_igrad,"%s/imposed-grad-y_fit-err.csv",indir);
      path = inpath_igrad;  
      rc = gm2fieldUtil::Import::ImportData2<double,double>(path,"csv",RAD,vertGrad_err);
   }else if(fitDim==2){
      // just load the 2D fit results giving the imposed gradient and errors in vectors
      // radial
      sprintf(inpath_igrad,"%s/imposed-grad-x_2D.csv",indir);
      path = inpath_igrad; 
      rc = gm2fieldUtil::Import::ImportData2<double,double>(path,"csv",igx,igxErr);
      // vertical 
      sprintf(inpath_igrad,"%s/imposed-grad-y_2D.csv",indir);
      path = inpath_igrad; 
      rc = gm2fieldUtil::Import::ImportData2<double,double>(path,"csv",igy,igyErr);
   }else{
      std::cout << "[Misalignment_prod::GetImposedGradients]: Invalid fit dimension!" << std::endl;
      return 1;
   }
   
   const int NPX = xPar.size();
   const int NPY = yPar.size();

   double rad_coord  = coord[0]; 
   double vert_coord = coord[1]; 
   double azi_coord  = coord[2]; 

   if(fitDim==1){
      // build fits
      fitRadGrad  = new TF1("fitRadGrad",MyFitFunc_pol2_simple,-35,35,NPX); 
      for(int i=0;i<NPX;i++){
	 fitRadGrad->SetParameter(i,xPar[i]); 
	 fitRadGrad->SetParError(i,xParErr[i]); 
      }

      fitVertGrad  = new TF1("fitVertGrad",MyFitFunc_pol2_simple,-35,35,NPY); 
      for(int i=0;i<NPY;i++){
	 fitVertGrad->SetParameter(i,yPar[i]); 
	 fitVertGrad->SetParError(i,yParErr[i]); 
      }
      // evaluate fit at probe coordinate to get gradients (Hz/mm)  
      rad_grad  = fitRadGrad->Eval(vert_coord);
      vert_grad = fitVertGrad->Eval(rad_coord);
      // get the appropriate error -- note the switched coordinates!  
      rad_grad_err  = GetGradientError(vert_coord,VER,radGrad_err);  
      vert_grad_err = GetGradientError(rad_coord ,RAD,vertGrad_err); 
   }else if(fitDim==2){
      rad_grad      = igx[probeNumber-1]; 
      rad_grad_err  = igxErr[probeNumber-1]; 
      vert_grad     = igy[probeNumber-1]; 
      vert_grad_err = igyErr[probeNumber-1]; 
   }else{
      std::cout << "[Misalignment_prod::GetImposedGradients]: Invalid fit dimension!" << std::endl;
      return 1;
   }

   // we do z differently 
   double azi_grad=0,azi_grad_err=0;
   sprintf(inpath_igrad,"%s/imposed-grad_z_pr-%02d.csv",indir,probeNumber);
   rc = LoadImposedAziGradData(inpath_igrad,azi_grad,azi_grad_err);

   // fill output vectors
   impG.push_back(rad_grad); 
   impG.push_back(vert_grad); 
   impG.push_back(azi_grad); 

   impGerr.push_back(rad_grad_err); 
   impGerr.push_back(vert_grad_err); 
   impGerr.push_back(azi_grad_err); 

   return 0;
}
