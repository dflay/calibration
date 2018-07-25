// Process the trolley data; apply a drift correction during the run, 
// produce an output file based on the input run list containing 
// the raw frequencies and the drift-corrected frequencies for all probes  

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
#include "./src/Consolidate.C"
#include "./src/CustomUtilities.C"
#include "./src/CustomMath.C"
#include "./src/CustomGraph.C"
#include "./src/CustomImport.C"
#include "./src/CustomExport.C"
#include "./src/CustomAlgorithms.C"

double gMarkerSize = 0.8; 

int CheckTrolleyData(std::vector<fixedProbeEvent_t> fxprData,std::vector<trolleyAnaEvent_t> &trlyData);
int GetStatsForSingleProbe(int probe,std::vector<trolleyAnaEvent_t> data,double &mean,double &stdev,double &min,double &max);
// int GetStatsForAllProbes(std::vector<trolleyAnaEvent_t> data,std::vector<double> &MEAN,std::vector<double> &STDEV);

int PrintToFileSingle(std::string outpath,int run,double time,double t,double x,double dx,double y,double dy);

int PrintToFile(std::string outpath,int run,
                std::vector<double> x,std::vector<double> dx,
                std::vector<double> y,std::vector<double> dy);

int ProcessTRLY(std::string configFile){

   int rc=0;
   int method = gm2fieldUtil::Constants::kPhaseDerivative;

   InputManager *inputMgr = new InputManager();
   inputMgr->Load(configFile);
   inputMgr->Print();

   std::string anaDate = inputMgr->GetAnalysisDate();
   bool isBlind        = inputMgr->IsBlind();
   int fxpr_tag        = inputMgr->GetFixedProbeListTag(); 

   date_t theDate;
   GetDate(theDate);

   char plotDir[200];
   sprintf(plotDir,"./plots/%02d-%02d-%02d",theDate.month,theDate.day,theDate.year-2000);
   rc = MakeDirectory(plotDir);

   char outDir[200];
   sprintf(outDir,"./output");
   if(isBlind)  sprintf(outDir,"%s/blinded"  ,outDir);
   if(!isBlind) sprintf(outDir,"%s/unblinded",outDir);
   sprintf(outDir,"%s/%02d-%02d-%02d",outDir,theDate.month,theDate.day,theDate.year-2000);
   rc = MakeDirectory(outDir);

   char outPath[500];
   sprintf(outPath,"%s/trly-data_%s.csv",outDir,anaDate.c_str());
   std::string outpath = outPath;

   blind_t blind;
   ImportBlinding(blind);
   double blindValue = blind.value_tr;

   std::vector<int> run;
   std::vector<std::string> label;
   inputMgr->GetRunList(run);
   inputMgr->GetRunLabels(label);

   const int NRUNS = run.size();
   
   std::cout << "Getting fixed probe data..." << std::endl;
   char fxpr_path[500]; 
   sprintf(fxpr_path,"./input/probe-lists/fxpr-list_set-%d.csv",fxpr_tag); 
   std::string fxprPath = fxpr_path;
   std::vector<int> fxprList;
   gm2fieldUtil::Import::ImportData1<int>(fxprPath,"csv",fxprList);

   std::vector<gm2field::fixedProbeFrequency_t> fxprData;

   // details for start and stop indices for fxpr data  
   int NN=0;
   int NFP = fxprList.size();
   int startIndex = fxprList[0];
   int stopIndex  = fxprList[NFP-1];
   unsigned long long t0     = 0;
   unsigned long long tStart=0,tStop=0,tStep=1E+9;

   std::vector<fixedProbeEvent_t> fxprDataAvg;

   TGraph **gFPAVG = new TGraph*[NRUNS]; 
   TString fxprPlots = Form("%s/fxpr_runs-%04d-%04d.png",plotDir,run[0],run[NRUNS-1]);

   // for graphing 
   double xMin = tStart/1E+9;
   double xMax = tStop/1E+9;

   double time,freq,freqErr,freq_cor,freqErr_cor; 
   std::vector<double> FREQ,FREQERR,FREQ_cor,FREQERR_cor; 
   std::vector<trolleyAnaEvent_t> trlyData,trlyDataCor;     

   TCanvas *c1 = new TCanvas("c1","TRLY plots",1200,1000);
   c1->Divide(4,5);

   TCanvas *c2 = new TCanvas("c2","FXPR plots",1200,1000);
   c2->Divide((NRUNS+1)/2,2);

   TString trlyPlots; 
  
   double yMin=0,yMax=0; 

   // for gathering info on trly probe positions  
   std::vector<double> rr,yy,zz;
   std::vector< std::vector<double> > r,y,z;  // first index is run, second is probe

   const int NTRLY = 17; 
   // main analysis loop 
   for(int i=0;i<NRUNS;i++){
      std::cout << "------------------ RUN " << run[i] << " ------------------" << std::endl;
      std::cout << "Getting fixed probe data..." << std::endl;
      rc = gm2fieldUtil::RootHelper::GetFPFrequencies(run[i],fxprData);
      NN = fxprData.size();  
      tStart = fxprData[0].GpsTimeStamp[startIndex];
      tStop  = fxprData[NN-1].GpsTimeStamp[stopIndex];
      GetAverageFXPRVectorsNew(method,t0,tStart,tStop,tStep,fxprList,fxprData,fxprDataAvg);
      gFPAVG[i] = GetTGraphNew(fxprDataAvg); 
      gm2fieldUtil::Graph::SetGraphParameters(gFPAVG[i],20,kBlack);
      std::cout << "Getting trolley data..." << std::endl;
      rc = GetTrolleyData("",run[i],method,trlyData);
      if(rc!=0){
	 std::cout << "No data!" << std::endl;
	 return 1;
      }
      if(isBlind) ApplyBlindingTRLY(blindValue,trlyData);
      std::cout << "--> Done." << std::endl;
      rc = CheckTrolleyData(fxprDataAvg,trlyData);
      // apply P2P drift correction
      std::cout << "Applying drift correction..." << std::endl;
      rc = CorrectTRLYForDriftDuringMeasurement(fxprDataAvg,trlyData,trlyDataCor);
      std::cout << "--> Done." << std::endl;
      // get stats for all trolley probes
      TGraph **gTRLY     = new TGraph*[NTRLY]; 
      TGraph **gTRLY_cor = new TGraph*[NTRLY]; 
      for(int j=0;j<NTRLY;j++){ 
         // get stats
         time = trlyData[i].time[j]/1E+9; 
	 rc = GetStatsForSingleProbe(j,trlyData   ,freq    ,freqErr    ,yMin,yMax);
	 rc = GetStatsForSingleProbe(j,trlyDataCor,freq_cor,freqErr_cor,yMin,yMax);
         // print to file for a single probe
	 sprintf(outPath,"%s/trly-data_pr-%02d_%s.csv",outDir,j,anaDate.c_str());
	 outpath = outPath;
	 rc = PrintToFileSingle(outpath,run[i],time,0,freq,freqErr,freq_cor,freqErr_cor); 
         // get probe positions
         rr.push_back( trlyData[0].r[j] ); 
         yy.push_back( trlyData[0].y[j] ); 
         zz.push_back( trlyData[0].phi[j] ); 
         // make plots
         gTRLY[j]     = GetTRLYTGraph(j,"GpsTimeStamp","freq",trlyData);
         gTRLY_cor[j] = GetTRLYTGraph(j,"GpsTimeStamp","freq",trlyDataCor);
         gm2fieldUtil::Graph::SetGraphParameters(gTRLY[j]    ,21,kBlack);
         gm2fieldUtil::Graph::SetGraphParameters(gTRLY_cor[j],20,kRed);
         gTRLY[j]->SetMarkerSize(gMarkerSize); 
         gTRLY_cor[j]->SetMarkerSize(gMarkerSize);
         c1->cd(j+1); 
         gTRLY[j]->SetTitle( Form("Probe %02d",j) );
         gTRLY[j]->Draw("alp");  
         gTRLY[j]->GetYaxis()->SetRangeUser(yMin,yMax); 
         gm2fieldUtil::Graph::UseTimeDisplay(gTRLY[j]); 
         gTRLY[j]->Draw("alp");  
         gTRLY_cor[j]->Draw("same"); 
         c1->Update();  
      }
      // save probe positions to vector-of-vector 
      r.push_back(rr);  
      y.push_back(yy);  
      z.push_back(zz);  
      // save trly plots 
      trlyPlots = Form("%s/trly_run-%05d.png",plotDir,run[i]);
      c1->cd(); 
      c1->Print(trlyPlots);
      // fxpr plots 
      c2->cd(i+1);
      gFPAVG[i]->Draw("alp"); 
      gm2fieldUtil::Graph::SetGraphLabels(gFPAVG[i],Form("FXPR Data Run %d",run[i]),"","Frequency (Hz)"); 
      gm2fieldUtil::Graph::UseTimeDisplay(gFPAVG[i]); 
      gFPAVG[i]->Draw("alp"); 
      c2->Update(); 
      // set up for next run
      trlyData.clear();
      trlyDataCor.clear(); 
      fxprData.clear();
      fxprDataAvg.clear();
      rr.clear(); 
      yy.clear(); 
      zz.clear(); 
   }

   c2->cd(); 
   c2->Print(fxprPlots);

   // now get averages for trly probes 
   double mean_r,mean_y,mean_z;
   double stdev_r,stdev_y,stdev_z;
   std::vector<int> probe; 
   std::vector<double> R,DR,Y,DY,Z,DZ; 
   for(int i=0;i<NTRLY;i++){
      probe.push_back(i); 
      // for a given trolley probe, sum over values 
      for(int j=0;j<NRUNS;j++){
	 rr.push_back(r[j][i]);
	 yy.push_back(y[j][i]);
	 zz.push_back(z[j][i]);
      }
      // now get the stats
      mean_r  = gm2fieldUtil::Math::GetMean<double>(rr);  
      stdev_r = gm2fieldUtil::Math::GetStandardDeviation<double>(rr);  
      mean_y  = gm2fieldUtil::Math::GetMean<double>(yy);  
      stdev_y = gm2fieldUtil::Math::GetStandardDeviation<double>(yy);  
      mean_z  = gm2fieldUtil::Math::GetMean<double>(zz);  
      stdev_z = gm2fieldUtil::Math::GetStandardDeviation<double>(zz);  
      // save results 
      R.push_back(mean_r); 
      DR.push_back(stdev_r); 
      Y.push_back(mean_y); 
      DY.push_back(stdev_y); 
      Z.push_back(mean_z); 
      DZ.push_back(stdev_z);
      // clear vectors 
      rr.clear();
      yy.clear();
      zz.clear(); 
   }

   char outpath_trly_pos[500]; 
   sprintf(outpath_trly_pos,"%s/trly-pos.csv",outDir); 
   PrintTRLYPositions(outpath_trly_pos,probe,R,DR,Y,DY,Z,DZ); 

   return 0;
}
//______________________________________________________________________________
int GetStatsForSingleProbe(int probe,std::vector<trolleyAnaEvent_t> data,double &mean,double &stdev,double &min,double &max){
   min = 60E+3;
   max = 0;
   double fLO = 61.74E+6; 
   const int N = data.size();
   std::vector<double> x; 
   for(int i=0;i<N;i++){
      x.push_back( fLO + data[i].freq[probe] );
      if(min>data[i].freq[probe]) min = data[i].freq[probe]; 
      if(max<data[i].freq[probe]) max = data[i].freq[probe]; 
   } 
   mean  = gm2fieldUtil::Math::GetMean<double>(x); 
   stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(x); 
   min  -= 5; 
   max  += 5; 
   return 0;
}
//______________________________________________________________________________
int PrintToFileSingle(std::string outpath,int run,double time,double temp,double x,double dx,double y,double dy){
   // print a single probe to file 
   // x = uncorrected, y = drift corrected 
   char myStr[1000]; 
   // prepare the line
   sprintf(myStr,"%d,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf",run,time,temp,x,dx,y,dy);

   std::ofstream outfile;
   outfile.open(outpath,ios::app);
   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      outfile << myStr << std::endl;
      std::cout << "The data has been written to file: " << outpath << std::endl;
   }  
   return 0;
}
//______________________________________________________________________________
int CheckTrolleyData(std::vector<fixedProbeEvent_t> fxprData,std::vector<trolleyAnaEvent_t> &trlyData){
   // check end time of trolley data; if it exceeds the FXPR, chop it down to size

   const int NT    = trlyData.size();
   const int NF    = fxprData.size();
   const int NTRLY = 17; 

   int cntr=0; 
   double theTime=0;
   double max_time = fxprData[NF-1].time;

   // count how many events we need to cut off 
   bool isBad=false;
   for(int i=0;i<NT;i++){
      for(int j=0;j<NTRLY;j++){
         theTime = trlyData[i].time[j]/1E+9; 
         // std::cout << gm2fieldUtil::GetStringTimeStampFromUTC(theTime)  << " " 
         //           << gm2fieldUtil::GetStringTimeStampFromUTC(max_time) << std::endl;
	 if(theTime>max_time) isBad = true;   // if any trolley probe has a higher time than the FXPR, count it 
      }
      if(isBad) cntr++;
      // reset for next event 
      isBad = false; 
   }

   if(cntr>0) std::cout << "[CheckTrolleyData]: Will remove " << cntr << " events" << std::endl;

   // now cut the vector down if needed 
   for(int i=0;i<cntr;i++) trlyData.pop_back(); 

   return 0; 

}
//______________________________________________________________________________
int PrintToFile(std::string outpath,int run,
                std::vector<double> x,std::vector<double> dx,
                std::vector<double> y,std::vector<double> dy){
   // print a single probe to file 
   // x = uncorrected, y = drift corrected 
   const int N = x.size();
   char myStr[1000]; 
   // prepare the first line
   sprintf(myStr,"%d,%.3lf,%.3lf,%.3lf,%.3lf",run,x[0],dx[0],y[0],dy[0]);
   // put together all the lines 
   for(int i=1;i<N-1;i++) sprintf(myStr,"%s,%.3lf,%.3lf,%.3lf,%.3lf",myStr,x[i],dx[i],y[i],dy[i]);
   sprintf(myStr,"%s,%.3lf,%.3lf,%.3lf,%.3lf",myStr,x[N-1],dx[N-1],y[N-1],dy[N-1]);

   std::ofstream outfile;
   outfile.open(outpath,ios::app);
   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      outfile << myStr << std::endl;
      std::cout << "The data has been written to file: " << outpath << std::endl;
   }  
   return 0;
}
