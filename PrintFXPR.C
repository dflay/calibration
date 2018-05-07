// Print fixed probe data to a csv file    

#include <cstdlib>
#include <vector>

#include "TCanvas.h"

#include "RootTreeStructs.h"
#include "gm2fieldConstants.h"
#include "gm2fieldMath.h"
#include "gm2fieldGraph.h"
#include "gm2fieldRootHelper.h"
#include "gm2fieldImport.h"
#include "gm2fieldExport.h"

int GetStats(int probeNumber,int method,std::vector<gm2field::fixedProbeFrequency_t> fxprData,double &mean,double &stdev); 

int PrintFXPR(){

   int runNumber=0;
   std::cout << "Enter run number: "; 
   std::cin  >> runNumber;

   int method = gm2fieldUtil::Constants::kPhaseDerivative; 

   // get all frequencies 
   std::vector<gm2field::fixedProbeFrequency_t> fxprData; 
   int rc = gm2fieldUtil::RootHelper::GetFPFrequencies(runNumber,fxprData); 
   if (rc!=0) {
      std::cout << "No data.  Exiting..." << std::endl;
      return 1;
   }

   // get header info (yoke, layer, etc)   
   std::vector<gm2field::fixedProbeHeader_t> fxprHeader;
   rc = gm2fieldUtil::RootHelper::GetFPHeader(runNumber,fxprHeader);
   if (rc!=0) {
      std::cout << "No data.  Exiting..." << std::endl;
      return 1;
   }

   // get position info 
   std::vector<gm2field::fixedProbePosition_t> fxprPos; 
   rc = gm2fieldUtil::RootHelper::GetFPPositions(runNumber,fxprPos); 
   if (rc!=0) {
      std::cout << "No data.  Exiting..." << std::endl;
      return 1;
   }
 
   double mean=0,stdev=0;

   std::string fxprPath = "./input/probe-lists/fxpr-list.csv";
   std::vector<int> probe;
   gm2fieldUtil::Import::ImportData1<int>(fxprPath,"csv",probe); 

   const int N = probe.size(); 
   for(int i=0;i<N;i++){
      // rc = gm2fieldUtil::Export::PrintFixedProbeFrequencyToCSV(runNumber,probe[i],method,fxprData);
      GetStats(probe[i],method,fxprData,mean,stdev);
      std::cout << "probe = "   << Form("%03d",probe[i])  
	        << " azi ID = " << fxprHeader[0].AziId[probe[i]] 
                << " azi angle = " << fxprPos[0].Phi[probe[i]] 
	        << " yoke = "   << (char)fxprHeader[0].YokeId[probe[i]] 
	        << " layer = "  << (char)fxprHeader[0].LayerId[probe[i]] 
	        << " rad ID = " << (char)fxprHeader[0].RadId[probe[i]] 
	        << " mean = "   << Form("%.3lf",mean) << " stdev = " << Form("%.3lf",stdev) << std::endl; 
   }
  
   std::cout << "Other proes to consider: "; 
   for(int i=0;i<N;i++){
      for(int j=0;j<378;j++){
	 if( fxprPos[0].Phi[j]>180 && fxprPos[0].Phi[j]<195 ){
	    std::cout << Form("angle = %.3lf probe %d",fxprPos[0].Phi[j],j) << std::endl;
	 }
      }
   }
 
   // std::vector<double> FREQ;
   // const int NEvents = fxprData.size();
   // for(int i=0;i<100;i++){
   //    for(int j=0;j<N;j++) FREQ.push_back( fxprData[i].Frequency[j][method] );
   //    mean  = gm2fieldUtil::Math::GetMean<double>(FREQ); 
   //    stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(FREQ);
   //    std::cout << mean << " " << stdev << std::endl; 
   // } 

   return 0;
}
//______________________________________________________________________________
int GetStats(int probeNumber,int method,std::vector<gm2field::fixedProbeFrequency_t> fxprData,double &mean,double &stdev){
   // compute mean and standard deviation of data for a given probe number
   const int N = fxprData.size();
   std::vector<double> freq; 
   for(int i=0;i<N;i++) freq.push_back(fxprData[i].Frequency[probeNumber][method]);
   mean  = gm2fieldUtil::Math::GetMean<double>(freq); 
   stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(freq); 
   return 0; 
} 
