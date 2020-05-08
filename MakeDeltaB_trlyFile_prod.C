// Open all Delta-B result files for the analyzed trolley Delta-B(x,y)
// data and produce a single file with results    

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

#include "./src/InputManager.C"
#include "./src/FitFuncs.C"
#include "./src/TRLYFuncs.C"
#include "./src/CalibFuncs.C"
#include "./src/CustomUtilities.C"
#include "./src/CustomMath.C"
#include "./src/CustomImport.C"
#include "./src/CustomExport.C"
#include "./src/CustomAlgorithms.C"
#include "./src/CustomGraph.C"

int PrintToFile_trlyDB(const char *outpath,std::vector<deltab_t> x,std::vector<deltab_t> y); 

int MakeDeltaB_trlyFile_prod(std::string configFile){

   int rc=0;
   int method = gm2fieldUtil::Constants::kPhaseDerivative;

   InputManager *inputMgr = new InputManager();
   inputMgr->Load(configFile);
   inputMgr->Print();

   std::string anaDate    = inputMgr->GetAnalysisDate();
   std::string blindLabel = inputMgr->GetBlindLabel(); 
   bool isBlind           = inputMgr->IsBlind();
   int runPeriod          = inputMgr->GetRunPeriod(); 

   date_t theDate;  
   GetDate(theDate);  

   std::string outDir  = GetPath("output",isBlind,blindLabel,theDate.getDateString());

   char outpath[500],inpath[500]; 
   sprintf(outpath,"./input/delta-b/trly_xy_run-%d.csv",runPeriod);

   std::vector<deltab_t> dBx,dBy;     

   for(int i=0;i<17;i++){
      // dBx data  
      sprintf(inpath,"%s/dB-trly_final-location_rad-grad_pr-%02d.csv",outDir.c_str(),i+1); 
      LoadDeltaBData(inpath,dBx);
      // dBy data  
      sprintf(inpath,"%s/dB-trly_final-location_vert-grad_pr-%02d.csv",outDir.c_str(),i+1); 
      LoadDeltaBData(inpath,dBy);
   }

   std::cout << "TRLY dBx VALS" << std::endl;
   const int NN = dBx.size();
   for(int i=0;i<NN;i++){
      std::cout << Form("raw = %.3lf +/- %.3lf, ABA = %.3lf +/- %.3lf",
                        dBx[i].dB,dBx[i].dB_err,dBx[i].dB_fxpr,dBx[i].dB_fxpr_err) << std::endl;
   }
   std::cout << "TRLY dBy VALS" << std::endl;
   for(int i=0;i<NN;i++){
      std::cout << Form("raw = %.3lf +/- %.3lf, ABA = %.3lf +/- %.3lf",
                        dBy[i].dB,dBy[i].dB_err,dBy[i].dB_fxpr,dBy[i].dB_fxpr_err) << std::endl;
   }

   rc = PrintToFile_trlyDB(outpath,dBx,dBy);
   
   delete inputMgr;

   return 0;
}
//______________________________________________________________________________
int PrintToFile_trlyDB(const char *outpath,std::vector<deltab_t> x,std::vector<deltab_t> y){
   char outStr[500]; 

   const int N = x.size();

   std::ofstream outfile;
   outfile.open(outpath); 
   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      for(int i=0;i<N;i++){
	 sprintf(outStr,"%02d,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf"
                       ,i+1,x[i].dB,x[i].dB_err,x[i].dB_fxpr,x[i].dB_fxpr_err,
                            y[i].dB,y[i].dB_err,y[i].dB_fxpr,y[i].dB_fxpr_err);
	 outfile << outStr << std::endl;
      }
      sprintf(outStr,"Data printed to file: %s",outpath);
      PrintMessage("MakeDeltaB_trlyFile_prod",outStr);
   }

   return 0;
}
