// Compare PP and TRLY measurements across multiple swaps back and forth  

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
#include "./src/CustomImport.C"
#include "./src/CustomExport.C"
#include "./src/CustomMath.C"
#include "./src/CustomGraph.C"
#include "./src/CustomAlgorithms.C"
#include "./src/OscFuncs.C"
#include "./src/InputManager.C"
#include "./src/FitFuncs.C"
// #include "./src/CalibFuncs.C"
#include "./src/TRLYFuncs.C"

double gMarkerSize = 0.8;
double gXSize      = 0.05;  
double gYSize      = 0.06;  

int PrintResults(std::string outpath,std::vector<double> trial,
                 std::vector<double> x1,std::vector<double> dx1,
                 std::vector<double> x2,std::vector<double> dx2);

int PrintResults(std::string outpath,std::vector<double> x,std::vector<double> dx);
int PrintResults(std::string outpath,std::vector<double> x1,std::vector<double> x2,std::vector<double> dx); 


int CalibrateMultiSwap_prod(std::string configFile){

   char msg[200]; 

   int rc=0;
   int method = gm2fieldUtil::Constants::kPhaseDerivative;

   InputManager *inputMgr = new InputManager();
   inputMgr->Load(configFile);
   inputMgr->Print();

   std::string runDate    = inputMgr->GetRunDate();
   std::string blindLabel = inputMgr->GetBlindLabel(); 
   bool isBlind           = inputMgr->IsBlind();   
   bool useTimeWeight     = inputMgr->GetTimeWeightStatus();
   bool isMisalignCor     = inputMgr->GetMisalignCorStatus();
   int probeNumber        = inputMgr->GetTrolleyProbe();
   // systematics 
   bool isSyst            = inputMgr->GetSystStatus();
   int  systDirNum        = inputMgr->GetSystDirNum();  

   // if(isMisalignCor){
      sprintf(msg,"[CalibrateMultiSwap_prod]: Will use misalignment correction"); 
      Logger::PrintMessage(Logger::kINFO,"default",msg,'a');
   // }

   // date_t theDate;
   // GetDate(theDate);
   std::string theDate = inputMgr->GetAnaDate();

   std::string plotDir = GetPath("plots" ,isBlind,blindLabel,theDate,isSyst,systDirNum);
   std::string outDir  = GetPath("output",isBlind,blindLabel,theDate,isSyst,systDirNum);

   char outPath_swap[500]; 
   sprintf(outPath_swap,"%s/results_swap-data_pr-%02d.csv",outDir.c_str(),probeNumber);
   std::string outpath_swap = outPath_swap; 
 
   char outPath_result[500]; 
   sprintf(outPath_result,"%s/results_pr-%02d.csv",outDir.c_str(),probeNumber);
   std::string outpath_result = outPath_result;  

   sprintf(outPath_swap,"%s/results_free-prot_swap-data_pr-%02d.csv",outDir.c_str(),probeNumber);
   std::string outpath_swap_free = outPath_swap; 
 
   sprintf(outPath_result,"%s/results_free-prot_pr-%02d.csv",outDir.c_str(),probeNumber);
   std::string outpath_result_free = outPath_result;  

   // load the PP data
   std::vector<calibSwap_t> ppData,ppData_free;
   std::vector<double> ppTime,ppFreq,ppFreqErr,ppFreqFree,ppFreqFreeErr,ppTemp; 
   char ppPath[200]; 
   sprintf(ppPath,"%s/pp-swap-data_pr-%02d.csv",outDir.c_str(),probeNumber); 
   LoadCalibSwapData(ppPath,ppData); 
   sprintf(ppPath,"%s/pp-swap-data_free-prot_pr-%02d.csv",outDir.c_str(),probeNumber); 
   LoadCalibSwapData(ppPath,ppData_free); 
   const int NPP = ppData.size();
   for(int i=0;i<NPP;i++){
      ppTime.push_back(ppData[i].time); 
      ppFreq.push_back(ppData[i].freq); 
      ppFreqErr.push_back(ppData[i].freqErr); 
      ppFreqFree.push_back(ppData_free[i].freq); 
      ppFreqFreeErr.push_back(ppData_free[i].freqErr);
      ppTemp.push_back(ppData_free[i].temp);  
   } 

   // mean PP temp 
   double pp_temp_mean  = gm2fieldUtil::Math::GetMean<double>(ppTemp); 
   double pp_temp_stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(ppTemp); 

   // load the trly data
   std::vector<calibSwap_t> trlyData; 
   std::vector<double> trTime,trFreq,trFreqErr,trTemp,trTempErr; 
   char trlyPath[200]; 
   sprintf(trlyPath,"%s/trly-swap-data_pr-%02d.csv",outDir.c_str(),probeNumber); 
   LoadCalibSwapData(trlyPath,trlyData); 
   const int NTR = trlyData.size();
   for(int i=0;i<NTR;i++){
      trTime.push_back(trlyData[i].time); 
      trFreq.push_back(trlyData[i].freq); 
      trFreqErr.push_back(trlyData[i].freqErr); 
      trTemp.push_back(trlyData[i].temp); 
   } 

   // load the trly data (temp corrected)
   std::vector<calibSwap_t> trlyData_tc; 
   std::vector<double> trTime_tc,trFreq_tc,trFreqErr_tc,trTemp_tc,trTempErr_tc; 
   sprintf(trlyPath,"%s/trly-swap-data_temp-cor_pr-%02d.csv",outDir.c_str(),probeNumber); 
   LoadCalibSwapData(trlyPath,trlyData_tc); 
   const int NTR_tc = trlyData_tc.size();
   for(int i=0;i<NTR_tc;i++){
      trTime_tc.push_back(trlyData_tc[i].time); 
      trFreq_tc.push_back(trlyData_tc[i].freq); 
      trFreqErr_tc.push_back(trlyData_tc[i].freqErr); 
      trTemp_tc.push_back(trlyData_tc[i].temp); 
   } 

   // mean trly temp 
   double tr_temp_mean  = gm2fieldUtil::Math::GetMean<double>(trTemp); 
   double tr_temp_stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(trTemp); 
 
   if(NPP!=NTR){
      std::cout << "WARNING: Number of swaps for PP and TRLY do not match!  Exiting..." << std::endl;
      std::cout << "PP swaps:   " << NPP << std::endl;
      std::cout << "TRLY swaps: " << NTR << std::endl;
      if(NPP>NTR){
	 std::cout << "Removing the last PP swap" << std::endl;
	 ppTime.pop_back();
	 ppFreq.pop_back();
	 ppFreqErr.pop_back();
	 ppFreqFree.pop_back();
	 ppFreqFreeErr.pop_back();
      }else{
	 std::cout << "Removing the last TRLY swap" << std::endl;
	 trTime.pop_back();
	 trFreq.pop_back();
	 trFreqErr.pop_back();
      }
   }
   
   // misalignment correction
   // index 0 = raw, 1 = ABA, 2 = opt, 3 = opt w/barcode for z  
   std::vector<misalignCor_t> mCor;
   char outPath_misalign[200]; 
   sprintf(outPath_misalign,"%s/misalign-cor_pr-%02d.csv",outDir.c_str(),probeNumber);
   rc = LoadMisalignmentCorData(outPath_misalign,mCor);
 
   // determine if the PP or TRLY came first 
   bool ppFirst = false;
   if( ppTime[0]<trTime[0] ) ppFirst = true; 

   // compute difference
   double err=0;
   double mean=0,mean_free=0; 
   double mean_aba=0,mean_free_aba=0; 
   double stdev=0,stdev_free=0; 
   double stdev_aba=0,stdev_free_aba=0; 
   double stdev_opt=0,stdev_free_opt=0; 

   // raw result
   std::vector<double> diff,diffErr; 
   rc = GetDifference(ppFreq,ppFreqErr,trFreq,trFreqErr,diff,diffErr);
   rc = GetWeightedAverageStats(diff,diffErr,mean,err,stdev); 

   // free proton result: uses free proton PP, temp corrected TRLY 
   std::vector<double> diff_free,diffErr_free; 
   rc = GetDifference(ppFreqFree,ppFreqFreeErr,trFreq_tc,trFreqErr_tc,diff_free,diffErr_free);
   rc = GetWeightedAverageStats(diff_free,diffErr_free,mean_free,err,stdev_free); 

   // ABA: raw PP, TRLY 
   std::vector<double> diff_aba,diffErr_aba; 
   if(ppFirst){
      // PP came first! A = PP, B = TRLY 
      rc = GetDifference_ABA_final(useTimeWeight,ppTime,ppFreq,ppFreqErr,trTime,trFreq,trFreqErr,diff_aba,diffErr_aba);
   }else{
      // TRLY came first! A = TRLY, B = PP 
      rc = GetDifference_ABA_final(useTimeWeight,trTime,trFreq,trFreqErr,ppTime,ppFreq,ppFreqErr,diff_aba,diffErr_aba);
      // must resverse sign of results since we want PP-TRLY 
      for(int i=0;i<NPP;i++) diff_aba[i] *= -1.; 
   } 

   rc = GetWeightedAverageStats(diff_aba,diffErr_aba,mean_aba,err,stdev_aba); 
   stdev_opt = stdev_aba; 

   // ABA: free proton PP, temp corrected TRLY
   std::vector<double> diff_free_aba,diffErr_free_aba;
   if(ppFirst){ 
      // PP came first! A = PP, B = TRLY 
      rc = GetDifference_ABA_final(useTimeWeight,ppTime,ppFreqFree,ppFreqFreeErr,trTime_tc,trFreq_tc,trFreqErr_tc,diff_free_aba,diffErr_free_aba);
   }else{
      // TRLY came first! A = TRLY, B = PP 
      rc = GetDifference_ABA_final(useTimeWeight,trTime_tc,trFreq_tc,trFreqErr_tc,ppTime,ppFreqFree,ppFreqFreeErr,diff_free_aba,diffErr_free_aba);
      // must reverse sign of results since we want PP-TRLY 
      for(int i=0;i<NPP;i++) diff_free_aba[i] *= -1.; 
   }

   rc = GetWeightedAverageStats(diff_free_aba,diffErr_free_aba,mean_free_aba,err,stdev_free_aba); 
   stdev_free_opt = stdev_free_aba; 

   // create an opt result, which is the same for ABA since we always have swaps 
   double mean_opt      = mean_aba; 
   double mean_free_opt = mean_free_aba; 

   // shielded proton
   double mean_cor              = mean          + mCor[0].val; // raw  
   double mean_cor_aba          = mean_aba      + mCor[1].val; // ABA 
   double mean_cor_opt          = mean_opt      + mCor[2].val; // opt
   double mean_cor_bar_opt      = mean_opt      + mCor[3].val; // opt + barcode
   // free proton result
   double mean_cor_free         = mean_free     + mCor[0].val; // raw  
   double mean_cor_free_aba     = mean_free_aba + mCor[1].val; // ABA 
   double mean_cor_free_opt     = mean_free_opt + mCor[2].val; // opt
   double mean_cor_free_bar_opt = mean_free_opt + mCor[3].val; // opt + barcode

   // if(isMisalignCor){
   //    // apply the misalignment correction
   //    mean          += mCor[0].val;  
   //    mean_aba      += mCor[1].val;  
   //    mean_opt      += mCor[2].val;  
   //    mean_free     += mCor[0].val;  
   //    mean_free_aba += mCor[1].val;  
   //    mean_free_opt += mCor[2].val; 
   // }

   std::cout << Form("======================= TRLY PROBE %02d RESULTS =======================",probeNumber) << std::endl;
   std::cout << "WITHOUT Misalignment corrections" << std::endl; 
   std::cout << "Bare Result" << std::endl;
   std::cout << Form("[RAW] mean = %.3lf +/- %.3lf Hz (%.3lf ppb)",mean    ,stdev,stdev/0.06179)         << std::endl;
   std::cout << Form("[ABA] mean = %.3lf +/- %.3lf Hz (%.3lf ppb)",mean_aba,stdev_aba,stdev_aba/0.06179) << std::endl;
   std::cout << Form("[opt] mean = %.3lf +/- %.3lf Hz (%.3lf ppb)",mean_opt,stdev_opt,stdev_opt/0.06179) << std::endl;
   std::cout << "Free-Proton Result (no free-proton errors yet)" << std::endl;
   std::cout << Form("[RAW] mean = %.3lf +/- %.3lf Hz (%.3lf ppb)",mean_free    ,stdev_free,stdev_free/0.06179)         << std::endl;
   std::cout << Form("[ABA] mean = %.3lf +/- %.3lf Hz (%.3lf ppb)",mean_free_aba,stdev_free_aba,stdev_free_aba/0.06179) << std::endl;
   std::cout << Form("[opt] mean = %.3lf +/- %.3lf Hz (%.3lf ppb)",mean_free_opt,stdev_free_opt,stdev_free_opt/0.06179) << std::endl;
   std::cout << "WITH Misalignment corrections (note: stat errors only)" << std::endl; 
   std::cout << "Bare Result" << std::endl;
   std::cout << Form("[RAW]    mean = %.3lf +/- %.3lf Hz (%.3lf ppb)",mean_cor        ,stdev,stdev/0.06179)         << std::endl;
   std::cout << Form("[ABA]    mean = %.3lf +/- %.3lf Hz (%.3lf ppb)",mean_cor_aba    ,stdev_aba,stdev_aba/0.06179) << std::endl;
   std::cout << Form("[opt]    mean = %.3lf +/- %.3lf Hz (%.3lf ppb)",mean_cor_opt    ,stdev_opt,stdev_opt/0.06179) << std::endl;
   std::cout << Form("[opt,bc] mean = %.3lf +/- %.3lf Hz (%.3lf ppb)",mean_cor_bar_opt,stdev_opt,stdev_opt/0.06179) << std::endl;
   std::cout << "Free-Proton Result (no free-proton errors yet)" << std::endl;
   std::cout << Form("[RAW]    mean = %.3lf +/- %.3lf Hz (%.3lf ppb)",mean_cor_free        ,stdev_free,stdev_free/0.06179)         << std::endl;
   std::cout << Form("[ABA]    mean = %.3lf +/- %.3lf Hz (%.3lf ppb)",mean_cor_free_aba    ,stdev_free_aba,stdev_free_aba/0.06179) << std::endl;
   std::cout << Form("[opt]    mean = %.3lf +/- %.3lf Hz (%.3lf ppb)",mean_cor_free_opt    ,stdev_free_opt,stdev_free_opt/0.06179) << std::endl;
   std::cout << Form("[opt,bc] mean = %.3lf +/- %.3lf Hz (%.3lf ppb)",mean_cor_free_bar_opt,stdev_free_opt,stdev_free_opt/0.06179) << std::endl;

   std::vector<double> trial,trial_aba,trial_opt; 
   const int ND = diff.size();
   for(int i=0;i<ND;i++)  trial.push_back(i+1); 
   const int NDA = diff_aba.size();
   for(int i=0;i<NDA;i++) trial_aba.push_back(i+1); 

   // make plots 
   TGraphErrors *gPP            = gm2fieldUtil::Graph::GetTGraphErrors(ppTime   ,ppFreq       ,ppFreqErr       );
   TGraphErrors *gPP_free       = gm2fieldUtil::Graph::GetTGraphErrors(ppTime   ,ppFreqFree   ,ppFreqFreeErr   );
   TGraphErrors *gTR            = gm2fieldUtil::Graph::GetTGraphErrors(trTime   ,trFreq       ,trFreqErr       );
   TGraphErrors *gDiff          = gm2fieldUtil::Graph::GetTGraphErrors(trial    ,diff         ,diffErr         );
   TGraphErrors *gDiff_aba      = gm2fieldUtil::Graph::GetTGraphErrors(trial_aba,diff_aba     ,diffErr_aba     );
   TGraphErrors *gDiff_free     = gm2fieldUtil::Graph::GetTGraphErrors(trial    ,diff_free    ,diffErr_free    );
   TGraphErrors *gDiff_free_aba = gm2fieldUtil::Graph::GetTGraphErrors(trial_aba,diff_free_aba,diffErr_free_aba);

   gm2fieldUtil::Graph::SetGraphParameters(gPP           ,20,kBlue   );
   gm2fieldUtil::Graph::SetGraphParameters(gPP_free      ,21,kBlue   );
   gm2fieldUtil::Graph::SetGraphParameters(gTR           ,20,kRed    );
   gm2fieldUtil::Graph::SetGraphParameters(gDiff         ,20,kBlack  );
   gm2fieldUtil::Graph::SetGraphParameters(gDiff_aba     ,20,kGreen+2);
   gm2fieldUtil::Graph::SetGraphParameters(gDiff_free    ,20,kBlack  );
   gm2fieldUtil::Graph::SetGraphParameters(gDiff_free_aba,20,kGreen+2);

   // add an extra entry to the ABA result to get the same length as the raw result 
   trial_aba.push_back(-1);
   diff_aba.push_back(0); 
   diffErr_aba.push_back(0); 
   diff_free_aba.push_back(0); 
   diffErr_free_aba.push_back(0); 

   // make vectors for main results 

   // storing temperatures -- avg over all trials 
   std::vector<double> result,resultErr,result_cor;  
   result.push_back(mean); 
   result.push_back(mean_aba); 
   result.push_back(mean_opt); 
   result.push_back(mean_opt); 
   resultErr.push_back(stdev); 
   resultErr.push_back(stdev_aba); 
   resultErr.push_back(stdev_opt); 
   resultErr.push_back(stdev_opt); 
   result_cor.push_back(mean_cor); 
   result_cor.push_back(mean_cor_aba); 
   result_cor.push_back(mean_cor_opt); 
   result_cor.push_back(mean_cor_bar_opt); 
   
   std::vector<double> result_free,resultErr_free,result_cor_free; 
   result_free.push_back(mean_free); 
   result_free.push_back(mean_free_aba); 
   result_free.push_back(mean_free_opt); 
   result_free.push_back(mean_free_opt); 
   resultErr_free.push_back(stdev_free); 
   resultErr_free.push_back(stdev_free_aba); 
   resultErr_free.push_back(stdev_free_opt); 
   resultErr_free.push_back(stdev_free_opt); 
   result_cor_free.push_back(mean_cor_free); 
   result_cor_free.push_back(mean_cor_free_aba); 
   result_cor_free.push_back(mean_cor_free_opt); 
   result_cor_free.push_back(mean_cor_free_bar_opt); 

   // print to file 
   rc = PrintResults(outpath_swap  ,trial,diff,diffErr,diff_aba,diffErr_aba); 
   rc = PrintResults(outpath_result,result,result_cor,resultErr); 

   rc = PrintResults(outpath_swap_free  ,trial,diff_free,diffErr_free,diff_free_aba,diffErr_free_aba); 
   rc = PrintResults(outpath_result_free,result_free,result_cor_free,resultErr_free); 

   TMultiGraph *mgp = new TMultiGraph(); 
   mgp->Add(gPP     ,"lp"); 
   // mgp->Add(gPP_free,"lp"); 
   mgp->Add(gTR     ,"lp"); 

   TMultiGraph *mgd = new TMultiGraph();
   mgd->Add(gDiff    ,"lp");
   mgd->Add(gDiff_aba,"lp");

   TMultiGraph *mgd_free = new TMultiGraph();
   mgd_free->Add(gDiff_free    ,"lp");
   mgd_free->Add(gDiff_free_aba,"lp");

   TLegend *L = new TLegend(0.6,0.6,0.8,0.8);
   L->AddEntry(gDiff    ,"Raw"       ,"p");    
   L->AddEntry(gDiff_aba,"ABA Method","p");   

   TCanvas *c1 = new TCanvas("c1","PP, TRLY plots",1200,800);
   c1->Divide(1,3);

   c1->cd(1);
   mgp->Draw("a"); 
   // gm2fieldUtil::Graph::SetGraphLabels(mgp,"PP (Blue), TRLY (Red) Data","",Form("Frequency-%.3lf MHz (Hz)",mean_freq/1E+6));
   gm2fieldUtil::Graph::SetGraphLabels(mgp,"PP (Blue), TRLY (Red) Data","","Frequency (Hz)");
   gm2fieldUtil::Graph::UseTimeDisplay(mgp);
   gm2fieldUtil::Graph::SetGraphLabelSizes(mgp,gXSize,gYSize);  
   mgp->Draw("a");
   c1->Update(); 

   c1->cd(2);
   mgd->Draw("ap"); 
   gm2fieldUtil::Graph::SetGraphLabels(mgd,"PP-TRLY Data (black = raw, green = ABA)","Trial","Frequency Difference (Hz)");
   gm2fieldUtil::Graph::SetGraphLabelSizes(mgd,gXSize,gYSize);  
   mgd->Draw("ap");
   c1->Update();
 
   c1->cd(3);
   mgd_free->Draw("ap"); 
   gm2fieldUtil::Graph::SetGraphLabels(mgd_free,"PP(free)-TRLY Data (black = raw, green = ABA)","Trial","Frequency Difference (Hz)");
   gm2fieldUtil::Graph::SetGraphLabelSizes(mgd_free,gXSize,gYSize);  
   mgd_free->Draw("ap");
   c1->Update(); 

   // save the plot  
   c1->cd(); 
   TString calibPlots = Form("%s/pp-trly_calib-data_pr-%02d_%s.png",plotDir.c_str(),probeNumber,runDate.c_str());
   c1->Print(calibPlots);

   return 0;
}
//______________________________________________________________________________
int PrintResults(std::string outpath,std::vector<double> trial,
                 std::vector<double> x1,std::vector<double> dx1,
                 std::vector<double> x2,std::vector<double> dx2){

   const int N = trial.size(); 
   char myStr[1000],header[500]; 
   // prepare the header 
   sprintf(header,"#trial,diff,err,diff_ABA,err_ABA");  

   std::ofstream outfile;
   outfile.open(outpath.c_str());
   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      outfile << header << std::endl; 
      for(int i=0;i<N;i++){
	 sprintf(myStr,"%02d,%.3lf,%.3lf,%.3lf,%.3lf",
                 (int)trial[i],x1[i],dx1[i],x2[i],dx2[i]);
	 outfile << myStr << std::endl;
      }
      std::cout << "The data has been written to file: " << outpath << std::endl;
   }  
   return 0;
}
//______________________________________________________________________________
int PrintResults(std::string outpath,
                std::vector<double> x,std::vector<double> dx){

   const int N = x.size(); 
   char myStr[1000],header[500]; 
   // prepare the header 
   sprintf(header,"#diff,err;first-row-is-raw;second-is-ABA");  

   std::ofstream outfile;
   outfile.open(outpath.c_str());
   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      outfile << header << std::endl; 
      for(int i=0;i<N;i++){
	 sprintf(myStr,"%.3lf,%.3lf",x[i],dx[i]);
	 outfile << myStr << std::endl;
      }
      std::cout << "The data has been written to file: " << outpath << std::endl;
   }  
   return 0;
}
//______________________________________________________________________________
int PrintResults(std::string outpath,
                std::vector<double> x1,std::vector<double> x2,std::vector<double> dx){

   const int N = x1.size(); 
   char myStr[1000],header[500]; 
   // prepare the header 
   sprintf(header,"#diff,diffCor,err;first-row-is-raw;second-is-ABA;third-is-opt;forth-is-bar-opt");  

   std::ofstream outfile;
   outfile.open(outpath.c_str());
   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      outfile << header << std::endl; 
      for(int i=0;i<N;i++){
	 sprintf(myStr,"%.3lf,%.3lf,%.3lf",x1[i],x2[i],dx[i]);
	 outfile << myStr << std::endl;
      }
      std::cout << "The data has been written to file: " << outpath << std::endl;
   }  
   return 0;
}
