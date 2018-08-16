// Determine the gradient in the field using the PP (scan) data    
// in the region close to the TRLY probe of interest  

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
#include "TFitResult.h"

#include "RootTreeStructs.h"
#include "gm2fieldMath.h"
#include "gm2fieldGraph.h"
#include "gm2fieldRootHelper.h"
#include "gm2fieldUnits.h"
#include "TemperatureSensor.h"

#include "./include/Constants.h"
#include "./include/drift.h"
#include "./include/plungingProbeAnaEvent.h"
#include "./include/trolleyAnaEvent.h"
#include "./include/fixedProbeEvent.h"

#include "./src/MyFits.C"
#include "./src/FitErr.C"
#include "./src/InputManager.C"
#include "./src/FXPRFuncs.C"
#include "./src/TRLYFuncs.C"
#include "./src/Consolidate.C"
#include "./src/CustomMath.C"
#include "./src/CustomGraph.C"
#include "./src/CustomImport.C"
#include "./src/CustomExport.C"
#include "./src/CustomAlgorithms.C"
#include "./src/CustomUtilities.C"
#include "./src/DeltaBFuncs.C"

int GetAziGrad(TGraph *g,double &grad,double &gradErr);

int LocalScanGrad_pp_prod(std::string configFile){

   std::cout << "------------------------------------" << std::endl;
   std::cout << "GET SHIMMED GRADIENT (USING PP SCAN)" << std::endl;

   int rc=0;
   int method = gm2fieldUtil::Constants::kPhaseDerivative;

   InputManager *inputMgr = new InputManager(); 
   inputMgr->UseAxis();         // need to grab the axis data in the JSON file 
   inputMgr->Load(configFile);
   inputMgr->Print(); 

   std::string date    = inputMgr->GetAnalysisDate(); 
   std::string fitFunc = inputMgr->GetValue("fit");
   bool isBlind        = inputMgr->IsBlind();
   int probeNumber     = inputMgr->GetTrolleyProbe(); 
   int axis            = inputMgr->GetAxis();
   int fxprSet         = inputMgr->GetFixedProbeListTag();  

   date_t theDate;
   GetDate(theDate);

   std::string Axis,gradType; 
   if(axis==0) Axis = "x"; 
   if(axis==1) Axis = "y"; 
   if(axis==2) Axis = "z"; 

   if( Axis.compare("x")==0 ) gradType = "rad"; 
   if( Axis.compare("y")==0 ) gradType = "vert"; 
   if( Axis.compare("z")==0 ) gradType = "azi";

   // make output directories 
   char outdir[200];
   if(isBlind)  sprintf(outdir,"./output/blinded/%s"  ,theDate.getDateString().c_str());
   if(!isBlind) sprintf(outdir,"./output/unblinded/%s",theDate.getDateString().c_str());
   rc = MakeDirectory(outdir);

   char plotDir[200];
   sprintf(plotDir,"./plots/%s",theDate.getDateString().c_str());
   MakeDirectory(plotDir);

   blind_t blind;
   ImportBlinding(blind);
   double blindValue = blind.value_pp;

   // outpaths  
   char outpath[200],tag[10];
   sprintf(tag    ,"%s-grad",gradType.c_str()); 
   sprintf(outpath,"%s/%s_pr-%02d.csv" ,outdir,tag,probeNumber);
   std::string strTag = tag;

   std::vector<int> run,subRun;
   std::vector<std::string> label;
   inputMgr->GetRunList(run); 
   inputMgr->GetRunLabels(label); 

   const int NRUN = run.size();
   if(NRUN==0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   // find the trolley run
   int midasRun=0; 
   for(int i=0;i<NRUN;i++){
      if(label[i].compare("midas-run")==0){
	 midasRun = run[i]; 
      }else{
	 subRun.push_back(run[i]);
      }
   }

   double ZERO[3] = {0,0,0}; 

   int NSR = subRun.size();
   for(int i=0;i<NSR;i++){
      if(subRun[i]==-1){
	 std::cout << "Invalid subrun number " << subRun[i] << "! There appears to be no scan data.  Exiting." << std::endl;
	 PrintToFile(outpath,strTag,ZERO,ZERO,ZERO,ZERO);
	 return 1;
      }
   }

   if(NSR<3){
      std::cout << "WARNING: Only two data points to fit!  Using a 1st-order polynomial instead of " << fitFunc << std::endl;
      fitFunc = "pol1";
   }

   // trolley probe coordinates 
   trolleyProbePosition_t trlyPos; 
   rc = GetTrolleyProbePositions(trlyPos);
   double pos[2] = {trlyPos.r[probeNumber-1],trlyPos.y[probeNumber-1]}; 

   // we do something different for the z axis (see below) 
   double trly_coord=0;
   if(axis!=2) trly_coord = pos[axis]; 

   // PP data 
   std::vector<plungingProbeAnaEvent_t> ppData,ppEvent; 
   std::cout << "Getting run " << midasRun << std::endl; 
   rc = GetPlungingProbeData(midasRun,method,ppData);
   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   const int NEV = ppData.size(); 
   std::cout << "PP events: " << NEV << std::endl;
   for(int i=0;i<NEV;i++) rc = ModifyPlungingProbeData(kLeastSquaresPhase,ppData[i]);
  
   if(isBlind) ApplyBlindingPP(blindValue,ppData);

   // now parse PP data using subrun list 
   rc = FilterPlungingProbeData(subRun,ppData,ppEvent); 
   std::cout << "Filtered PP data to use NMR-DAQ runs: " << std::endl;
   int M=0;
   const int NPP = ppEvent.size();
   double mean=0,stdev=0;
   std::vector<double> F;
   for(int i=0;i<NPP;i++){
      M = ppEvent[i].numTraces; 
      for(int j=0;j<M;j++) F.push_back(ppEvent[i].freq[j]); 
      mean  = gm2fieldUtil::Math::GetMean<double>(F); 
      stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(F); 
      std::cout << Form("%d: x = %.3lf mm, y = %.3lf mm, z = %.3lf mm, F = %.3lf +/- %.3lf Hz",
                        ppEvent[i].run,ppEvent[i].r[0],ppEvent[i].y[0],ppEvent[i].phi[0],
                        mean,stdev) << std::endl;
      // clean up for next event 
      F.clear();
   } 

   // fixed probe data
   char fxpr_path[200]; 
   sprintf(fxpr_path,"./input/probe-lists/fxpr-list_set-%d.csv",fxprSet); 
   std::string fxprPath = fxpr_path;
   std::vector<int> fxprList;
   gm2fieldUtil::Import::ImportData1<int>(fxprPath,"csv",fxprList);

   std::vector<gm2field::fixedProbeFrequency_t> fxprData;
   rc = gm2fieldUtil::RootHelper::GetFPFrequencies(midasRun,fxprData);
   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }

   const int NN = fxprData.size();
   const int NFP = fxprList.size();
   unsigned long long t0     = 0;
   unsigned long long tStart = fxprData[0].GpsTimeStamp[fxprList[0]];
   unsigned long long tStop  = fxprData[NN-1].GpsTimeStamp[fxprList[NFP-1]];
   unsigned long long tStep  = 1E+9;

   std::vector<fixedProbeEvent_t> fxprDataAvg;
   GetAverageFXPRVectorsNew(method,t0,tStart,tStop,tStep,fxprList,fxprData,fxprDataAvg);

   // Plunging probe plots 
   TGraphErrors *g = GetPPTGraphErrors2(Axis.c_str(),"freq",ppEvent);
   // TGraphErrors *g = GetPPScanGraph(Axis.c_str(),"freq",ppEvent,trly_coord);
   gm2fieldUtil::Graph::SetGraphParameters(g,20,kBlack);

   // Fixed probe plot 
   TGraph *gFXPR = GetTGraphNew(fxprDataAvg);
   gm2fieldUtil::Graph::SetGraphParameters(gFXPR,20,kBlack);  

   TString xAxisLabel = Form("%s (mm)",Axis.c_str());

   TCanvas *c1 = new TCanvas("c1","PP Data",1200,600);
  
   c1->cd(); 
   g->Draw("ap");
   gm2fieldUtil::Graph::SetGraphLabels(g,"Shimmed Field",xAxisLabel,"Frequency (Hz)"); 
   g->Draw("ap");
   TFitResultPtr fitResult = g->Fit(fitFunc.c_str(),"QS"); 
   c1->Update(); 

   TString plotPath = Form("%s/pp-shimmed-scan-%s_pr-%02d.png",plotDir,Axis.c_str(),probeNumber); 
   c1->cd();
   c1->Print(plotPath); 

   TCanvas *c2 = new TCanvas("c2","FXPR Data",1200,600);

   c2->cd(); 
   gFXPR->Draw("alp");
   gm2fieldUtil::Graph::SetGraphLabels(gFXPR,"Fixed Probe Average","","Frequency (Hz)"); 
   gm2fieldUtil::Graph::UseTimeDisplay(gFXPR);  
   gFXPR->Draw("alp");
   c2->Update(); 

   plotPath = Form("%s/pp-shimmed-scan_fxpr-avg_pr-%02d.png",plotDir,probeNumber); 
   c2->cd();
   c2->Print(plotPath);

   // get fit info 
   TF1 *myFit = g->GetFunction(fitFunc.c_str()); 

   const int NPAR = myFit->GetNpar();
   double par[NPAR],parErr[NPAR]; 
   for(int i=0;i<NPAR;i++){
      par[i]       = myFit->GetParameter(i); 
      parErr[i]    = myFit->GetParError(i); 
   }

   // get gradient based on fit type and which trolley probe we're looking at  
   double dBdr=0,dBdr_err=0,dBdr_cor=0,dBdr_cor_err=0;
   double PR[3],ER[3];
   int slopeIndex=0;

   double sum_sq=0;

   double X[3] = {trly_coord,0.,0.}; 
 
   if(axis!=2){
      std::cout << "The trolley coordinate is " << trly_coord << " mm" << std::endl;
      dBdr     = myFit->Derivative(trly_coord);  // evaluate at the trolley probe position of interest 
      dBdr_err = GetFitError(myFit,fitResult,MyPolyFitFuncDerivative,X);   
      // if( fitFunc.compare("pol1")==0 ){
      //    dBdr_err = parErr[1]; 
      // }else if( fitFunc.compare("pol2")==0 ){
      //    sum_sq   = TMath::Power(parErr[1],2.)  
      //             + TMath::Power(2.*parErr[2]*trly_coord,2.);   
      //    dBdr_err = TMath::Sqrt(sum_sq);  
      // }else if( fitFunc.compare("pol3")==0 ){
      //    sum_sq   = TMath::Power(parErr[1],2.)  
      //             + TMath::Power(2.*parErr[2]*trly_coord,2.)  
      //             + TMath::Power(9.*parErr[3]*trly_coord*trly_coord,2.); 
      //    dBdr_err = TMath::Sqrt(sum_sq);  
      // }
   }else{
      rc = GetAziGrad(g,dBdr,dBdr_err);
   }

   // WARNING: not using the drift-corrected result.  Those probes are too far away!  
   PR[0] = dBdr; 
   PR[1] = dBdr;  
   PR[2] = 0;

   ER[0] = dBdr_err; 
   ER[1] = dBdr_err; 
   ER[2] = 0;

   double drift[3]     = {0,0,0};
   double drift_err[3] = {0,0,0};

   // assume the difference as the drift error for now 
   drift_err[1] = TMath::Abs(PR[1]-PR[0]); 

   std::cout << "============================ RESULTS ============================" << std::endl;
   std::cout << Form("%s:        %.3lf +/- %.3lf Hz/mm",gradType.c_str(),PR[0],ER[0]) << std::endl;
   std::cout << Form("%s [fxpr]: %.3lf +/- %.3lf +/- %.3lf Hz/mm",gradType.c_str(),PR[1],ER[1],drift_err[1]) << std::endl;
   std::cout << Form("%s [trly]: %.3lf +/- %.3lf +/- %.3lf Hz/mm",gradType.c_str(),PR[2],ER[2],drift_err[2]) << std::endl;

   PrintToFile(outpath,strTag,PR,ER,drift,drift_err);

   delete inputMgr; 
   
   return 0;
}
// //______________________________________________________________________________
// double GetPPOrigin(int axis,std::vector<plungingProbeAnaEvent_t> data){
//    // find the origin of the scan 
//    std::vector<double> x; 
//    const int N = data.size();
//    for(int i=0;i<N;i++){
//       if(axis==0) x.push_back(data[i].r[0]); 
//       if(axis==1) x.push_back(data[i].y[0]); 
//       if(axis==2) x.push_back(data[i].phi[0]); 
//    }
//  
//    double x0=0; 
//    // sort in ascending order 
//    std::sort(x.begin(),x.end()); 
//    if(N>2){
//       x0 = x[1];
//    }else{
//       // how do we tell what is the origin here? 
//       x0 = x[0];
//       if(x[0]>x[1]) x0 = x[1];
//    }
//    return x0; 
// }
//______________________________________________________________________________
int GetAziGrad(TGraph *g,double &grad,double &gradErr){
   // will probably always be two data points
   // just take slope between points 
   const int N = g->GetN(); 
   double *x   = g->GetX();
   double *y   = g->GetY();
   double *ey  = g->GetEY();

   double df   = y[1]-y[0]; 
   double dx   = x[1]-x[0];

   // compute gradient 
   grad    = df/dx; 
   gradErr = TMath::Abs(grad);   // take absolute value 

   return 0; 
}
