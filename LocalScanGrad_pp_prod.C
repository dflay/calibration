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
#include "./src/SystFuncs.C"

int GetAziGrad(TGraph *g,double &grad,double &gradErr);
int GetLinearGrad(TGraph *g,double &grad,double &gradErr); 
int GetCoordinates(std::vector<calibSwap_t> data,std::vector<double> &r); 
int RedefineOrigin(std::vector<plungingProbeAnaEvent_t> &data,std::vector<double> r0); 

int GetShimmedGradient_alt(int probe,int axis,std::string outDir,
                           std::vector<double> par,std::vector<double> parErr,
                           std::vector<double> par_syst,std::vector<double> parErr_syst,
                           double &grad,double &grad_err); 

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
   std::string date          = inputMgr->GetRunDate(); 
   std::string fitFunc       = inputMgr->GetValue("fit");
   std::string fitFunc_syst  = "pol1";    // use to estimate systematic uncertainty
   std::string blindLabel    = inputMgr->GetBlindLabel();
   std::string cutFile       = inputMgr->GetCutFile();

   bool isBlind              = inputMgr->IsBlind();
   bool useOscCor            = inputMgr->GetOscCorStatus();
   bool shimGradAlt          = inputMgr->GetShimGradAltStatus(); 
   int probeNumber           = inputMgr->GetTrolleyProbe(); 
   int axis                  = inputMgr->GetAxis();
   int runPeriod             = inputMgr->GetRunPeriod();
   // systematics 
   bool isSyst               = inputMgr->GetSystStatus();
   bool varyFit              = inputMgr->GetSystFitStatus("shim");
   int systDirNum            = inputMgr->GetSystDirNum();

   double tempCorValue = 0;
   bool useTempCor_pp  = inputMgr->GetTempCorStatus_pp();
   if(useTempCor_pp) tempCorValue = inputMgr->GetTempCor_pp(); 

   char cutPath[200];
   sprintf(cutPath,"./input/json/run-%d/%s",runPeriod,cutFile.c_str());
   std::string cutpath = cutPath;

   // date_t theDate;
   // GetDate(theDate);
   std::string theDate = inputMgr->GetAnaDate();

   std::string Axis,gradType; 
   if(axis==0) Axis = "x"; 
   if(axis==1) Axis = "y"; 
   if(axis==2) Axis = "z"; 

   if( Axis.compare("x")==0 ) gradType = "rad"; 
   if( Axis.compare("y")==0 ) gradType = "vert"; 
   if( Axis.compare("z")==0 ) gradType = "azi";

   // make output directories 
   std::string plotDir = GetPath("plots" ,isBlind,blindLabel,theDate,isSyst,systDirNum);
   std::string outDir  = GetPath("output",isBlind,blindLabel,theDate,isSyst,systDirNum);

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
   // std::cout << "Found " << MR << " MIDAS runs" << std::endl;
   // for(int i=0;i<MR;i++) std::cout << mRun[i] << std::endl; 

   double ZERO[3] = {0,0,0}; 

   char msg[200];
   int NSR = subRun.size();
   // std::cout << "Found " << NSR << " subruns" << std::endl;
   for(int i=0;i<NSR;i++){
      // std::cout << subRun[i] << std::endl;
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
   bool useNMRANA = true;
   std::vector<plungingProbeAnaEvent_t> ppEvent,ppInput,ppInput2,ppInput3; 
   for(int i=0;i<MR;i++){
      std::cout << "Getting PP data from MIDAS run " << mRun[i] << "..." << std::endl; 
      rc = GetPlungingProbeData(mRun[i],prMethod,ppMethod,ppInput,prodVersion,nmrAnaVersion,cutPath,useNMRANA,tempCorValue);
      if(rc!=0){
	 std::cout << "No data!" << std::endl;
	 return 1;
      }
   }

   if(isBlind) ApplyBlindingPP(blindValue,ppInput);

   // const int NI = ppInput.size();
   // for(int i=0;i<NI;i++){
   //    std::cout << Form("PP DAQ run %d, trace %d, freq = %.3lf Hz",ppInput[i].run,ppInput[i].traceNumber[0],ppInput[i].freq[0]) << std::endl;
   // }
 
   // now parse PP data using subrun list 
   rc = FilterPlungingProbeData(subRun,ppInput,ppInput2); 

   std::vector<int> fxprList;
   inputMgr->GetFXPRList(fxprList);

   // special time cut for the FXPR data (run 2 only)
   // time cut for the FXPR data
   char cut_path[200];
   sprintf(cut_path,"./input/json/run-%d/extra-cuts.json",runPeriod); 
   cutpath = cut_path;
   Cut *myCut = new Cut();
   unsigned long long tMin=0,tMax=-1;
   rc = GetFXPRCutTime(cutpath,probeNumber,axis,tMin,tMax);

   // filter PP data further if necessary
   if(probeNumber==8 &&runPeriod==2){
      rc = myCut->FilterPPData(runPeriod,probeNumber,"shim",Axis,ppInput2,ppInput3,cutpath);
   }else{
      CopyPlungingProbe(ppInput2,ppInput3);
   }

   delete myCut;

   bool subtractDrift = inputMgr->GetFXPRDriftStatus();  
   int period = inputMgr->GetNumEventsTimeWindow(); 
   std::vector<averageFixedProbeEvent_t> fxprData;
   for(int i=0;i<MR;i++){
      rc = GetFixedProbeData_avg(mRun[i],prMethod,fxprList,fxprData,prodVersion,subtractDrift,period,tMin,tMax);
      if(rc!=0){
         std::cout << "No data!" << std::endl;
         return 1;
      }
   }

   // reference time for oscillation correction
   // need first time of first PP event  
   double t0_osc = ppInput3[0].time[0]/1E+9;

   // oscillation correction
   if(useOscCor){
      rc = CorrectOscillation_pp(fxprData,ppInput3,ppEvent,t0_osc);
   }else{
      CopyPlungingProbe(ppInput3,ppEvent);
   }

   // trolley probe coordinates 
   // std::vector<calibSwap_t> trData; 
   // trolleyProbePosition_t trlyPos; 
   // rc = GetTrolleyProbePositions(trlyPos);
   // double pos[2] = {trlyPos.r[probeNumber-1],trlyPos.y[probeNumber-1]};
 
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
      std::cout << Form("PP DAQ Run %d: x = %.3lf mm, y = %.3lf mm, z = %.3lf mm, F = %.3lf +/- %.3lf Hz",
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
   sprintf(inpath_pp,"%s/pp-swap-data_pr-%02d.csv",outDir.c_str(),probeNumber);
   rc = LoadCalibSwapData(inpath_pp,ppCalib_data);  

   // trolley swap data 
   std::vector<calibSwap_t> trlyCalib_data;
   char inpath_tr[200]; 
   sprintf(inpath_tr,"%s/trly-swap-data_pr-%02d.csv",outDir.c_str(),probeNumber);
   rc = LoadCalibSwapData(inpath_tr,trlyCalib_data);  

   // get the AVG coordinates of the PP and trly probe during rapid swapping  
   std::vector<double> ppCoord,trlyCoord;
   rc = GetCoordinates(ppCalib_data  ,ppCoord  );  
   rc = GetCoordinates(trlyCalib_data,trlyCoord);  
 
   // redefine PP origin based on average location in rapid swapping.  We want
   // the gradient at this point   
   rc = RedefineOrigin(ppEvent,ppCoord);

   // std::cout << "[LocalScanGrad_pp_prod]: Coordinates during rapid swapping:" << std::endl;
   // std::cout << Form("                         PP: %.3lf mm, TRLY = %.3lf mm",ppCoord[axis],trlyCoord[axis]) << std::endl;

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
   
   TGraphErrors *g2 = GetPPTGraphErrors2(Axis.c_str(),"freq",ppEvent);
   gm2fieldUtil::Graph::SetGraphParameters(g2,20,kBlack);

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
   c1->Divide(1,2); 
  
   c1->cd(1); 
   g->Draw("ap");
   gm2fieldUtil::Graph::SetGraphLabels(g,"Shimmed Field (Gradient Evaluated at Dotted Line)",xAxisLabel,"Frequency (Hz)"); 
   g->Draw("ap");
   g->GetYaxis()->SetRangeUser(min,max); 
   TFitResultPtr fitResult = g->Fit(fitFunc.c_str(),"QS"); 
   evalPoint->Draw("same"); 
   c1->Update(); 

   c1->cd(2); 
   g2->Draw("ap");
   gm2fieldUtil::Graph::SetGraphLabels(g,"Shimmed Field (Linear Fit)",xAxisLabel,"Frequency (Hz)"); 
   g2->Draw("ap");
   g2->GetYaxis()->SetRangeUser(min,max); 
   TFitResultPtr fitResult_syst = g2->Fit(fitFunc_syst.c_str(),"QS"); 
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
   TF1 *myFit      = g->GetFunction(fitFunc.c_str()); 
   TF1 *myFit_syst = g2->GetFunction(fitFunc_syst.c_str()); 

   std::cout << "Standard fit (pol2) results: " << std::endl;
   const int NPAR = myFit->GetNpar();
   std::vector<std::string> parLabel; 
   std::vector<double> par,parErr; 
   for(int i=0;i<NPAR;i++){
      parLabel.push_back( "p"+std::to_string(i) ); 
      par.push_back( myFit->GetParameter(i) ); 
      parErr.push_back( myFit->GetParError(i) );
      std::cout << Form("%s = %.3lf +/- %.3lf",parLabel[i].c_str(),par[i],parErr[i]) << std::endl; 
   }

   std::cout << "Linear fit results: " << std::endl;
   const int NPAR_syst = myFit_syst->GetNpar();

   // systematic uncertainty estimate 
   std::vector<std::string> parLabel_syst; 
   std::vector<double> par_syst,parErr_syst; 

   for(int i=0;i<NPAR_syst;i++){
      parLabel_syst.push_back( "p"+std::to_string(i) ); 
      par_syst.push_back( myFit_syst->GetParameter(i) ); 
      parErr_syst.push_back( myFit_syst->GetParError(i) );
      std::cout << Form("%s = %.3lf +/- %.3lf",parLabel_syst[i].c_str(),par_syst[i],parErr_syst[i]) << std::endl; 
   }
   std::cout << "------------------" << std::endl;

   // save parameters to file
   char outpath_par[200]; 
   sprintf(outpath_par,"%s/shimmed-grad-%s_pars_pr-%02d.csv",outDir.c_str(),Axis.c_str(),probeNumber); 
   rc = PrintToFile(outpath_par,parLabel,par,parErr);  

   // save parameters to file
   sprintf(outpath_par,"%s/shimmed-grad-%s_pol1-pars_pr-%02d.csv",outDir.c_str(),Axis.c_str(),probeNumber); 
   rc = PrintToFile(outpath_par,parLabel_syst,par_syst,parErr_syst);  

   // get gradient based on fit type and which trolley probe we're looking at  
   double dBdr=0,dBdr_err=0,dBdr_cor=0,dBdr_cor_err=0;
   double PR[3],ER[3];
   int slopeIndex=0;

   double sum_sq=0;

   double XX[3] = {X,0.,0.}; 
   
   int NPTS = ppEvent.size(); // how many data points do we have? 

   // std::cout << "The trolley coordinate is " << X << " mm" << std::endl;
   if(NPTS>2){
      dBdr     = myFit->Derivative(X);  // evaluate at the trolley probe position of interest 
      dBdr_err = GetFitError(myFit,fitResult,MyPolyFitFuncDerivative,XX);   
    }else if(NPTS==2){
       // two or less points, use a linear fit 
       rc = GetLinearGrad(g,dBdr,dBdr_err);
    }else{
       std::cout << "[LocalScanGrad_pp_prod]: ERROR! Insufficient data points = " 
                 << NPTS << "!  Assuming 1000 Hz/mm" << std::endl;
       dBdr     = 1000.; 
       dBdr_err = 1000.; 
    }

   double dBdr_alt=0,dBdr_alt_err=0;
   if(shimGradAlt){
      std::cout << "[LocalScanGrad_pp_prod]: Using ALTERNATE evaluation!" << std::endl;
      rc = GetShimmedGradient_alt(probeNumber,axis,outDir,par,parErr,par_syst,parErr_syst,dBdr_alt,dBdr_alt_err);
      if(dBdr_alt==-1000){
	 // fall back to old method if this one fails
	 dBdr_alt     = dBdr; 
	 dBdr_alt_err = dBdr_err; 
       }
   }
 
   // if(axis!=2){
   //    // std::cout << "The trolley coordinate is " << X << " mm" << std::endl;
   //    dBdr     = myFit->Derivative(X);  // evaluate at the trolley probe position of interest 
   //    dBdr_err = GetFitError(myFit,fitResult,MyPolyFitFuncDerivative,XX);   
   //    // if( fitFunc.compare("pol1")==0 ){
   //    //    dBdr_err = parErr[1]; 
   //    // }else if( fitFunc.compare("pol2")==0 ){
   //    //    sum_sq   = TMath::Power(parErr[1],2.)  
   //    //             + TMath::Power(2.*parErr[2]*X,2.);   
   //    //    dBdr_err = TMath::Sqrt(sum_sq);  
   //    // }else if( fitFunc.compare("pol3")==0 ){
   //    //    sum_sq   = TMath::Power(parErr[1],2.)  
   //    //             + TMath::Power(2.*parErr[2]*X,2.)  
   //    //             + TMath::Power(9.*parErr[3]*X*X,2.); 
   //    //    dBdr_err = TMath::Sqrt(sum_sq);  
   //    // }
   // }else{
   //    rc = GetAziGrad(g,dBdr,dBdr_err);
   // }

   // WARNING: not using the drift-corrected result.  Those probes are too far away!  
   PR[0] = dBdr; 
   PR[1] = dBdr;  
   PR[2] = 0;

   ER[0] = dBdr_err; 
   ER[1] = dBdr_err; 
   ER[2] = 0;

   if(shimGradAlt){
      PR[0] = dBdr_alt; 
      PR[1] = dBdr_alt;  
      ER[0] = dBdr_alt_err; 
      ER[1] = dBdr_alt_err; 
   }

   if(isSyst && varyFit){
      std::cout << "[LocalScanGrad_pp_prod]: SYSTEMATIC VARIATION! Will vary fit result within (Gaussian) uncertainties" << std::endl;
      rc = systFunc::RandomizeFitValues(3,PR,ER); 
   }

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
int GetLinearGrad(TGraph *g,double &grad,double &gradErr){
   // for the case where we have two data points 
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
//______________________________________________________________________________
int GetShimmedGradient_alt(int probe,int axis,std::string outDir,
                           std::vector<double> par,std::vector<double> parErr,
                           std::vector<double> par_syst,std::vector<double> parErr_syst,
                           double &grad,double &grad_err){
   // alternative way of getting the shimmed gradient
   // evaluate shimmed gradient using shimmed field fit parameters of polynomial form 
   // f(q) = p0 + p1*q + p2*q^2
   // evaluate at q0 = delta dB/imp_grad  

   int NN = par.size();
   for(int i=0;i<NN;i++){
      std::cout << Form("%d, %.3lf, %.3lf",i,par[i],parErr[i]) << std::endl;
   }
   std::cout << "----" << std::endl;
   NN = par_syst.size();
   for(int i=0;i<NN;i++){
      std::cout << Form("%d, %.3lf, %.3lf",i,par_syst[i],parErr_syst[i]) << std::endl;
   }

   int rc=0;

   std::string gradName="none"; 
   if(axis==0) gradName = "rad"; 
   if(axis==1) gradName = "vert"; 
   if(axis==2) gradName = "azi"; 
     
   // load PP, TRLY delta-B data
   std::vector<deltab_t> pp,trly;
   char inpath[200];
   sprintf(inpath,"%s/dB-pp_final-location_%s-grad_pr-%02d.csv",outDir.c_str(),gradName.c_str(),probe); 
   LoadDeltaBData(inpath,pp); 

   sprintf(inpath,"%s/dB-trly_final-location_%s-grad_pr-%02d.csv",outDir.c_str(),gradName.c_str(),probe); 
   LoadDeltaBData(inpath,trly);

   // load imposed gradient data 
   // x and y (all probes) 
   // radial
   std::vector<double> igx,igxe;
   sprintf(inpath,"%s/imposed-grad-x_2D.csv",outDir.c_str());
   std::string path = inpath;
   rc = gm2fieldUtil::Import::ImportData2<double,double>(path,"csv",igx,igxe);
   // vertical 
   std::vector<double> igy,igye;
   sprintf(inpath,"%s/imposed-grad-y_2D.csv",outDir.c_str());
   path = inpath;
   rc = gm2fieldUtil::Import::ImportData2<double,double>(path,"csv",igy,igye); 
   // azi 
   double azi_grad=0,azi_grad_err=0;
   sprintf(inpath,"%s/imposed-grad_z_pr-%02d.csv",outDir.c_str(),probe);
   rc = LoadImposedAziGradData(inpath,azi_grad,azi_grad_err);

   double dB_tr=0,dB_pp=0,dB_tr_err=0,dB_pp_err=0;
   double dBdq_imp=0,dBdq_imp_err=0; 
 
   // delta-B (ABA is default)  
   dB_pp     = pp[axis].dB_fxpr;
   dB_pp_err = pp[axis].dB_fxpr_err;
   dB_tr     = trly[axis].dB_fxpr; 
   dB_tr_err = trly[axis].dB_fxpr_err; 
   // use raw if necessary
   if(dB_tr==0){
      dB_tr     = trly[axis].dB;  
      dB_tr_err = trly[axis].dB_err;
   }  
   if(dB_pp==0){
      dB_pp     = pp[axis].dB; 
      dB_pp_err = pp[axis].dB_err; 
   } 

   // imposed gradient 
   if(axis==0){
      dBdq_imp     = igx[probe-1]; 
      dBdq_imp_err = igxe[probe-1]; 
   }else if(axis==1){
      dBdq_imp     = igy[probe-1]; 
      dBdq_imp_err = igye[probe-1]; 
   }else if(axis==2){
      dBdq_imp     = azi_grad; 
      dBdq_imp_err = azi_grad_err; 
   }

   std::cout << Form("dB(TR) = %.3lf +/- %.3lf Hz",dB_tr,dB_tr_err)       << std::endl; 
   std::cout << Form("dB(PP) = %.3lf +/- %.3lf Hz",dB_pp,dB_pp_err)       << std::endl; 
   std::cout << Form("dBi/dq = %.3lf +/- %.3lf Hz",dBdq_imp,dBdq_imp_err) << std::endl; 

   double q0     = (dB_tr-dB_pp)/dBdq_imp;  // NOTE the Delta dB order! We evaluate where the TRLY is!  
   double q0_err = 0.;                      // FIXME: probably need an uncertainty due to q0... 
                                            // error on Delta-B is supressed by dB/dq, and dB/dq is usually very good.  

   double fitErr=0,fitDiff=0;

   const int NP = par.size();
   if(NP<2){
      std::cout << "[LocalScanGrad_pp_prod::GetShimmedGradient_alt]: ERROR!  Number of fit parameters < 2!" << std::endl;
      grad     = -1000; 
      grad_err =  1000; 
   }else if(NP==2){
      // linear fit 
      grad     = par[1];
      grad_err = parErr[1];
   }else if(NP==3){
      // polynomial fit, order 2
      grad     = par[1] + 2.*par[2]*q0;
      fitErr   = TMath::Sqrt( parErr[1]*parErr[1] + 2.*2.*parErr[2]*parErr[2] );
      fitDiff  = TMath::Abs(grad-par_syst[1]);  // consider difference in the gradients of higher order fit and linear fit as an error
      grad_err = TMath::Sqrt( fitErr*fitErr + fitDiff*fitDiff );  // quadrature sum 
      std::cout << Form("[LocalScanGrad_pp_prod::GetShimmedGradient_alt]: grad = %.3lf +/- %.3lf",grad,grad_err) << std::endl;
   }
   return 0;
}
