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
#include "TLine.h"

#include "RootTreeStructs.h"
#include "gm2fieldMath.h"
#include "gm2fieldGraph.h"
#include "gm2fieldRootHelper.h"
#include "gm2fieldUnits.h"
#include "TemperatureSensor.h"
#include "MovingAverage.h"
#include "Blinder.h"

#include "./include/Constants.h"
#include "./include/drift.h"
#include "./include/plungingProbeAnaEvent.h"
#include "./include/trolleyAnaEvent.h"
#include "./include/fixedProbeEvent.h"

#include "./src/CustomUtilities.C"
#include "./src/CustomImport.C"
#include "./src/CustomExport.C"
#include "./src/CustomAlgorithms.C"
#include "./src/CustomMath.C"
#include "./src/CustomGraph.C"
#include "./src/OscFuncs.C" 
#include "./src/BlindFuncs.C"
#include "./src/MyFits.C"
#include "./src/FitErr.C"
#include "./src/InputManager.C"
#include "./src/TRLYFuncs.C"

int GetAziGrad(TGraph *g,double &grad,double &gradErr);
int GetCoordinates(std::vector<calibSwap_t> data,std::vector<double> &r); 
int RedefineOrigin(std::vector<plungingProbeAnaEvent_t> &data,std::vector<double> r0); 

int LocalScanGrad_pp_prod(std::string configFile){

   std::cout << "------------------------------------" << std::endl;
   std::cout << "GET SHIMMED GRADIENT (USING PP SCAN)" << std::endl;

   int rc=0;
   int prMethod = gm2fieldUtil::Constants::kPhaseDerivative;
   int ppMethod = plungingProbeAnalysis::kLeastSquaresPhase;

   InputManager *inputMgr = new InputManager(); 
   inputMgr->UseAxis();         // need to grab the axis data in the JSON file 
   inputMgr->Load(configFile);
   inputMgr->Print(); 

   std::string prodVersion   = inputMgr->GetProductionTag();
   std::string nmrAnaVersion = inputMgr->GetNMRANATag();
   std::string date          = inputMgr->GetAnalysisDate(); 
   std::string fitFunc       = inputMgr->GetValue("fit");
   std::string blindLabel    = inputMgr->GetBlindLabel();
   std::string cutFile       = inputMgr->GetCutFile();

   bool isBlind              = inputMgr->IsBlind();
   bool useOscCor            = inputMgr->GetOscCorStatus();
   int probeNumber           = inputMgr->GetTrolleyProbe(); 
   int axis                  = inputMgr->GetAxis();
   int runPeriod             = inputMgr->GetRunPeriod();  

   char cutPath[200];
   sprintf(cutPath,"./input/json/run-%d/%s",runPeriod,cutFile.c_str());
   std::string cutpath = cutPath;

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
   std::string plotDir = GetPath("plots" ,isBlind,blindLabel,theDate.getDateString());
   std::string outDir  = GetPath("output",isBlind,blindLabel,theDate.getDateString());

   int blindUnits  = inputMgr->GetBlindUnits(); 
   double blindMag = inputMgr->GetBlindScale(); 
   gm2fieldUtil::Blinder *myBlind = new gm2fieldUtil::Blinder(blindLabel,blindMag,blindUnits);
   double blindValue = myBlind->GetBlinding(1); // in Hz

   // outpaths  
   char outpath[200],tag[10];
   sprintf(tag    ,"%s-grad",gradType.c_str()); 
   sprintf(outpath,"%s/%s_pr-%02d.csv" ,outDir.c_str(),tag,probeNumber);
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
   std::vector<int> mRun; 
   for(int i=0;i<NRUN;i++){
      if(label[i].compare("midas-run")==0){
	 midasRun = run[i];
	 mRun.push_back(run[i]);  
      }else{
	 subRun.push_back(run[i]);
      }
   }

   int MR = mRun.size();
   std::cout << "Found " << MR << " MIDAS runs" << std::endl; 

   double ZERO[3] = {0,0,0}; 

   char msg[200];
   int NSR = subRun.size();
   for(int i=0;i<NSR;i++){
      if(subRun[i]==-1){
	 sprintf(msg,"[LocalScanGrad_pp_prod]: Invalid subrun number %d for probe %02d!",subRun[i],probeNumber);
	 Logger::PrintMessage(Logger::kERR,"default",msg,'a');
	 PrintToFile(outpath,strTag,ZERO,ZERO,ZERO,ZERO);
	 return 1;
      }
   }

   if(NSR<3){
      std::cout << "WARNING: Only two data points to fit!  Using a 1st-order polynomial instead of " << fitFunc << std::endl;
      fitFunc = "pol1";
   }

   // PP data 
   std::vector<plungingProbeAnaEvent_t> ppData,ppEvent,ppInput; 
   std::cout << "Getting PP data for run " << midasRun << "..." << std::endl; 
   for(int i=0;i<MR;i++) rc = GetPlungingProbeData(mRun[i],prMethod,ppMethod,ppInput,prodVersion,nmrAnaVersion,cutPath);
   if(rc!=0){
      std::cout << "No data!" << std::endl;
      return 1;
   }
  
   if(isBlind) ApplyBlindingPP(blindValue,ppInput); 

   // reference time for oscillation correction
   // need first time of first PP event  
   double t0 = ppInput[0].time[0]/1E+9;

   std::vector<int> fxprList;
   inputMgr->GetFXPRList(fxprList);

   bool subtractDrift = inputMgr->GetFXPRDriftStatus();  
   int period = inputMgr->GetNumEventsTimeWindow(); // 10;
   std::vector<averageFixedProbeEvent_t> fxprData;
   for(int i=0;i<MR;i++){
      rc = GetFixedProbeData_avg(mRun[i],prMethod,fxprList,fxprData,prodVersion,subtractDrift,period,0);
      if(rc!=0){
         std::cout << "No data!" << std::endl;
         return 1;
      }
   }

   // oscillation correction 
   if(useOscCor){
      rc = CorrectOscillation_pp(fxprData,ppInput,ppData,t0);
   }else{
      CopyPlungingProbe(ppInput,ppData);
   }

   // trolley probe coordinates 
   // std::vector<calibSwap_t> trData; 
   // trolleyProbePosition_t trlyPos; 
   // rc = GetTrolleyProbePositions(trlyPos);
   // double pos[2] = {trlyPos.r[probeNumber-1],trlyPos.y[probeNumber-1]}; 

   // now parse PP data using subrun list 
   rc = FilterPlungingProbeData(subRun,ppData,ppEvent); 
   std::cout << "Filtered PP data to use NMR-DAQ runs: " << std::endl;
   int M=0;
   const int NPP = ppEvent.size();
   double mean_x=0,mean_y=0,mean_z=0,mean=0,stdev=0,max=0,min=50E+3;
   std::vector<double> xx,yy,zz,F;
   for(int i=0;i<NPP;i++){
      M = ppEvent[i].numTraces; 
      for(int j=0;j<M;j++){
	 xx.push_back(ppEvent[i].r[j]);
	 yy.push_back(ppEvent[i].y[j]);
	 zz.push_back(ppEvent[i].phi[j]);
	 F.push_back(ppEvent[i].freq[j]);
      } 
      mean_x  = gm2fieldUtil::Math::GetMean<double>(xx); 
      mean_y  = gm2fieldUtil::Math::GetMean<double>(yy); 
      mean_z  = gm2fieldUtil::Math::GetMean<double>(zz); 
      mean    = gm2fieldUtil::Math::GetMean<double>(F); 
      stdev   = gm2fieldUtil::Math::GetStandardDeviation<double>(F);
      if(max<mean) max = mean;  
      if(min>mean) min = mean;  
      std::cout << Form("%d: x = %.3lf mm, y = %.3lf mm, z = %.3lf mm, F = %.3lf +/- %.3lf Hz",
                        ppEvent[i].run,mean_x,mean_y,mean_z,mean,stdev) << std::endl;
      // clean up for next event 
      xx.clear();
      yy.clear();
      zz.clear();
      F.clear();
   } 

   min -= 8; 
   max += 8; 

   // plunging probe swap data 
   std::vector<calibSwap_t> ppCalib_data;
   char inpath_pp[200]; 
   sprintf(inpath_pp,"%s/pp-swap-data_pr-%02d_%s.csv",outDir.c_str(),probeNumber,date.c_str());
   rc = LoadCalibSwapData(inpath_pp,ppCalib_data);  

   // trolley swap data 
   std::vector<calibSwap_t> trlyCalib_data;
   char inpath_tr[200]; 
   sprintf(inpath_tr,"%s/trly-swap-data_pr-%02d_%s.csv",outDir.c_str(),probeNumber,date.c_str());
   rc = LoadCalibSwapData(inpath_tr,trlyCalib_data);  

   // get the coordinates of the PP and trly probe during rapid swapping  
   std::vector<double> ppCoord,trlyCoord;
   rc = GetCoordinates(ppCalib_data  ,ppCoord  );  
   rc = GetCoordinates(trlyCalib_data,trlyCoord);  
   
   // redefine PP origin based on average location in rapid swapping.  We want
   // the gradient at this point   
   rc = RedefineOrigin(ppEvent,ppCoord);

   // Want to evaluate the gradient at the average location of the PP. 
   // not sure if this is right... 
   // we do something different for the z axis (see below) 
   // double X = (-1.)*(ppCoord[axis] - trlyCoord[axis]);  // WARNING: we impose a NEGATIVE gradient, so our convention flips here
   // if(axis!=2) X = 0;
   // std::cout << "[LocalScanGrad_pp_prod]: Evaluating gradient at z = " << X << std::endl;
   // this seems like the most obvious choice, since we redefine the origin according the PP swap location (average)  
   double X = 0;
   
   TLine *evalPoint = new TLine(X,min,X,max);
   evalPoint->SetLineColor(kBlue);
   evalPoint->SetLineStyle(2); 
   evalPoint->SetLineWidth(2); 

   // Plunging probe plots 
   TGraphErrors *g = GetPPTGraphErrors2(Axis.c_str(),"freq",ppEvent);
   // TGraphErrors *g = GetPPScanGraph(Axis.c_str(),"freq",ppEvent,X);
   gm2fieldUtil::Graph::SetGraphParameters(g,20,kBlack);

   // Fixed probe plot 
   TGraph *gFXPR = GetFXPRTGraph_avg("GpsTimeStamp","freq","NONE",fxprData);
   gm2fieldUtil::Graph::SetGraphParameters(gFXPR,20,kBlack);  

   TString xAxisLabel;
   if(ppCoord[axis]>=0){
      xAxisLabel = Form("%s - %.3lf (mm)",Axis.c_str(),ppCoord[axis]);
   }else{
      xAxisLabel = Form("%s + %.3lf (mm)",Axis.c_str(),TMath::Abs(ppCoord[axis]) );
   }

   TCanvas *c1 = new TCanvas("c1","PP Data",1200,600);
  
   c1->cd(); 
   g->Draw("ap");
   gm2fieldUtil::Graph::SetGraphLabels(g,"Shimmed Field (Gradient Evaluated at Dotted Line)",xAxisLabel,"Frequency (Hz)"); 
   g->Draw("ap");
   g->GetYaxis()->SetRangeUser(min,max); 
   TFitResultPtr fitResult = g->Fit(fitFunc.c_str(),"QS"); 
   evalPoint->Draw("same"); 
   c1->Update(); 

   TString plotPath = Form("%s/pp-shimmed-scan-%s_pr-%02d.png",plotDir.c_str(),Axis.c_str(),probeNumber); 
   c1->cd();
   c1->Print(plotPath); 

   TCanvas *c2 = new TCanvas("c2","FXPR Data",1200,600);

   c2->cd(); 
   gFXPR->Draw("alp");
   gm2fieldUtil::Graph::SetGraphLabels(gFXPR,"Fixed Probe Average","","Frequency (Hz)"); 
   gm2fieldUtil::Graph::UseTimeDisplay(gFXPR);  
   gFXPR->Draw("alp");
   c2->Update(); 

   plotPath = Form("%s/pp-shimmed-scan-%s_fxpr-avg_pr-%02d.png",plotDir.c_str(),Axis.c_str(),probeNumber); 
   c2->cd();
   c2->Print(plotPath);

   // get fit info 
   TF1 *myFit = g->GetFunction(fitFunc.c_str()); 

   std::cout << "Fit Parameters: " << std::endl;
   const int NPAR = myFit->GetNpar();
   double par[NPAR],parErr[NPAR]; 
   for(int i=0;i<NPAR;i++){
      par[i]       = myFit->GetParameter(i); 
      parErr[i]    = myFit->GetParError(i);
      std::cout << Form("p[%d] = %.3lf +/- %.3lf",i,par[i],parErr[i]) << std::endl; 
   }

   // get gradient based on fit type and which trolley probe we're looking at  
   double dBdr=0,dBdr_err=0,dBdr_cor=0,dBdr_cor_err=0;
   double PR[3],ER[3];
   int slopeIndex=0;

   double sum_sq=0;

   double XX[3] = {X,0.,0.}; 
 
   if(axis!=2){
      std::cout << "The trolley coordinate is " << X << " mm" << std::endl;
      dBdr     = myFit->Derivative(X);  // evaluate at the trolley probe position of interest 
      dBdr_err = GetFitError(myFit,fitResult,MyPolyFitFuncDerivative,XX);   
      // if( fitFunc.compare("pol1")==0 ){
      //    dBdr_err = parErr[1]; 
      // }else if( fitFunc.compare("pol2")==0 ){
      //    sum_sq   = TMath::Power(parErr[1],2.)  
      //             + TMath::Power(2.*parErr[2]*X,2.);   
      //    dBdr_err = TMath::Sqrt(sum_sq);  
      // }else if( fitFunc.compare("pol3")==0 ){
      //    sum_sq   = TMath::Power(parErr[1],2.)  
      //             + TMath::Power(2.*parErr[2]*X,2.)  
      //             + TMath::Power(9.*parErr[3]*X*X,2.); 
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
//______________________________________________________________________________
int RedefineOrigin(std::vector<plungingProbeAnaEvent_t> &data,std::vector<double> r0){
   int M=0;
   const int N = data.size();
   for(int i=0;i<N;i++){
      M = data[i].numTraces; 
      for(int j=0;j<M;j++){
	 data[i].r[j]   -= r0[0];
	 data[i].y[j]   -= r0[1];
	 data[i].phi[j] -= r0[2];
      }  
   }
   return 0;
}
//______________________________________________________________________________
int GetCoordinates(std::vector<calibSwap_t> data,std::vector<double> &r){
   const int N = data.size();
   std::vector<double> x,y,z; 
   for(int i=0;i<N;i++){
      x.push_back( data[i].r   ); 
      y.push_back( data[i].y   ); 
      z.push_back( data[i].phi ); 
   }

   double mean_x = gm2fieldUtil::Math::GetMean<double>(x); 
   double mean_y = gm2fieldUtil::Math::GetMean<double>(y); 
   double mean_z = gm2fieldUtil::Math::GetMean<double>(z);

   r.push_back(mean_x);  
   r.push_back(mean_y);  
   r.push_back(mean_z);  

   return 0;

}
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
