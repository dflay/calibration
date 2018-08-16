// Process the PP data  
// produce an output file based on the input run list containing 
// the time stamp and raw frequencies  

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
#include "MovingAverage.h"

#include "./include/date.h"
#include "./include/Constants.h"
#include "./include/fixedProbeEvent.h"
#include "./include/trolleyAnaEvent.h"
#include "./include/sccEvent.h" 

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

int GetStatsForPP(perturbation_t ppPert,
                  std::vector<plungingProbeAnaEvent_t> data,std::vector<double> &time,
                  std::vector<double> &freq,std::vector<double> &freqErr,
                  std::vector<double> &freq_free,std::vector<double> &freqErr_free,
                  std::vector<double> &temp,std::vector<double> &tempErr,
                  double &min,double &max);

int PrintToFile(std::string outpath,std::vector<double> Time,
                std::vector<double> Freq,std::vector<double> dFreq,std::vector<double> Temp,std::vector<double> dTemp);

int Process_pp_prod(std::string configFile){

   int rc=0;
   int method = gm2fieldUtil::Constants::kPhaseDerivative;

   InputManager *inputMgr = new InputManager();
   inputMgr->Load(configFile);
   inputMgr->Print();

   std::string anaDate = inputMgr->GetAnalysisDate();
   bool isBlind        = inputMgr->IsBlind();
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

   char outPath[500]; 
   sprintf(outPath,"%s/pp-swap-data_pr-%02d_%s.csv",outDir,probeNumber,anaDate.c_str());
   std::string outpath_raw = outPath;  
   sprintf(outPath,"%s/pp-swap-data_free-prot_pr-%02d_%s.csv",outDir,probeNumber,anaDate.c_str());
   std::string outpath_free = outPath;  

   blind_t blind;
   ImportBlinding(blind);
   double blindValue = blind.value_pp;

   std::vector<int> run;
   std::vector<std::string> label;
   inputMgr->GetRunList(run);
   inputMgr->GetRunLabels(label);

   const int NRUNS = run.size();

   double yMin=0,yMax=0;

   std::vector<double> Freq,FreqErr,Temp,TempErr,Time;

   TString ppPlotPath;  

   // PP data 
   std::vector<plungingProbeAnaEvent_t> ppData;

   int NPP=0; 
   for(int i=0;i<NRUNS;i++){
      std::cout << "Getting PP data for run " << run[i] << "..." << std::endl;
      rc = GetPlungingProbeData(run[i],method,ppData);
      if(rc!=0){
	 std::cout << "No data!" << std::endl;
	 return 1;
      }
      // bring in your frequencies 
      NPP = ppData.size();
      for(int j=0;j<NPP;j++) rc = ModifyPlungingProbeData(kLeastSquaresPhase,ppData[j]);
      std::cout << "--> Done." << std::endl;
   }

   if(isBlind) ApplyBlindingPP(blindValue,ppData);

   // load in perturbation data 
   char inpath[200]; 
   perturbation_t ppPert;
   sprintf(inpath,"./input/perturbation/pp-pert.csv");
   LoadPerturbationData(inpath,ppPert);

   int M=0;

   std::vector<double> FreqFree,FreqFreeErr;

   // get stats for the PP 
   rc = GetStatsForPP(ppPert,ppData,Time,Freq,FreqErr,FreqFree,FreqFreeErr,Temp,TempErr,yMin,yMax); 

   // print to file 
   rc = PrintToFile(outpath_raw ,Time,Freq    ,FreqErr    ,Temp,TempErr); 
   rc = PrintToFile(outpath_free,Time,FreqFree,FreqFreeErr,Temp,TempErr); 

   // make plot of all frequencies for all shots 
   TGraph *gPP     = GetPPTGraph1("TimeStamp","freq",ppData); 
   gm2fieldUtil::Graph::SetGraphParameters(gPP,20,kBlack);
   gPP->SetMarkerSize(gMarkerSize); 

   TGraph *gPP_free = gm2fieldUtil::Graph::GetTGraphErrors(Time,FreqFree,FreqFreeErr);
   gm2fieldUtil::Graph::SetGraphParameters(gPP_free,20,kBlack);
   gPP_free->SetMarkerSize(gMarkerSize); 
   
   TCanvas *c1 = new TCanvas("c1","PP Data",1200,600);
   c1->Divide(1,2);

   c1->cd(1); 
   gPP->SetTitle( Form("Run %d",run[0]) );
   gPP->Draw("alp");  
   gPP->GetYaxis()->SetRangeUser(yMin,yMax);
   gm2fieldUtil::Graph::UseTimeDisplay(gPP); 
   gPP->Draw("alp");  
   c1->Update();

   c1->cd(2); 
   gPP_free->SetTitle( Form("Run %d",run[0]) );
   gPP_free->Draw("alp");  
   gm2fieldUtil::Graph::UseTimeDisplay(gPP_free); 
   gPP_free->Draw("alp");  
   c1->Update();

   // save the plot  
   c1->cd(); 
   ppPlotPath = Form("%s/pp-swap-data_pr-%02d.png",plotDir,probeNumber);
   c1->Print(ppPlotPath);

   return 0;
}
//______________________________________________________________________________
int GetStatsForPP(perturbation_t ppPert,
                  std::vector<plungingProbeAnaEvent_t> data,std::vector<double> &time,
                  std::vector<double> &freq,std::vector<double> &freqErr,
                  std::vector<double> &freq_free,std::vector<double> &freqErr_free,
                  std::vector<double> &temp,std::vector<double> &tempErr,
                  double &min,double &max){

   // gather stats for the PP; recall a single PP-DAQ run has N traces in it. 
   // here, a PP-DAQ run is a single swap!

   min = 50E+3; 
   max = 0; 
 
   int M=0,nmrDAQ_run=0;
   const int N = data.size();
 
   double arg=0;
   double mean_freq=0,stdev_freq=0;
   double mean_temp=0,stdev_temp=0;
   double freqFree=0,freqFreeErr=0;
   std::vector<double> x,xt;
 
   for(int i=0;i<N;i++){
      M = data[i].numTraces;
      nmrDAQ_run = data[i].run; 
      if(nmrDAQ_run==888 || nmrDAQ_run==834 || nmrDAQ_run==789){
	 std::cout << "WARNING: Known issue in run " << nmrDAQ_run << "!  Skipping..." << std::endl;
	 continue;
      } 
      for(int j=0;j<M;j++){
         // we want to add back in the LO since we're comparing against the TRLY eventually 
	 arg = data[i].freq[j] + data[i].freq_LO[j]; 
	 x.push_back(arg);
         xt.push_back(data[i].temp[j]); 
         if( min>data[i].freq[j] ) min = data[i].freq[j];  
         if( max<data[i].freq[j] ) max = data[i].freq[j];  
      }
      // get average on the run (swap) 
      mean_freq  = gm2fieldUtil::Math::GetMean<double>(x); 
      stdev_freq = gm2fieldUtil::Math::GetStandardDeviation<double>(x);
      mean_temp  = gm2fieldUtil::Math::GetMean<double>(xt); 
      stdev_temp = gm2fieldUtil::Math::GetStandardDeviation<double>(xt);
      // apply free proton corrections 
      GetOmegaP_free(ppPert,mean_freq,stdev_freq,mean_temp,stdev_temp,freqFree,freqFreeErr);
      // fill output vectors 
      time.push_back( data[i].time[0]/1E+9 ); 
      freq.push_back(mean_freq);
      freqErr.push_back(stdev_freq); 
      freq_free.push_back(freqFree);
      // WARNING: We'll fold in free-proton errors seprately.  we want shot errors here
      freqErr_free.push_back(stdev_freq);   
      temp.push_back(mean_temp); 
      tempErr.push_back(stdev_temp);
      // clean up for next run 
      x.clear();
      xt.clear();  
   }

   min -= 2; 
   max += 2;

   return 0;
}
//______________________________________________________________________________
int PrintToFile(std::string outpath,std::vector<double> Time,
                std::vector<double> Freq,std::vector<double> dFreq,std::vector<double> Temp,std::vector<double> dTemp){

   const int N = Time.size(); 
   char myStr[1000]; 

   std::ofstream outfile;
   outfile.open( outpath.c_str() );
   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      for(int i=0;i<N;i++){
	 sprintf(myStr,"%.0lf,%.3lf,%.3lf,%.3lf,%.3lf",Time[i],Freq[i],dFreq[i],Temp[i],dTemp[i]);
	 outfile << myStr << std::endl;
      } 
      std::cout << "The data has been written to file: " << outpath << std::endl;
      outfile.close();
   }  
   return 0;
}