// Compare PP and TRLY Delta-B measurements, imposed and shimmed gradients to obtain 
// misalignment uncertainties   
// TODO 
// - Combine Delta-B, shimmed field, imposed gradients to find error due to misalignment 

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

int GetShimmedGradients(const char *prefix,std::string date,int probe,std::vector<int> run,
                        std::vector<std::string> gradName,std::vector<grad_meas_t> &shim_grad_avg);

int Misalignment(std::string configFile){

   int rc=0;
   int method = gm2fieldUtil::Constants::kPhaseDerivative;

   InputManager *inputMgr = new InputManager();
   inputMgr->Load(configFile);
   inputMgr->Print();

   std::string fitFunc = "pol2";  // FIXME: should add this to the input manager... 
   std::string anaDate = inputMgr->GetAnalysisDate();
   bool isBlind        = inputMgr->IsBlind();
   int probe           = inputMgr->GetTrolleyProbe(); 
 
   std::vector<int> run;
   inputMgr->GetRunList(run); 
   const int NRUNS = run.size();

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

   char outPath[500],outPath_fxpr[500]; 
   sprintf(outPath     ,"%s/misalignment_results_%s.csv"     ,outDir,anaDate.c_str());
   sprintf(outPath_fxpr,"%s/misalignment_results_fxpr_%s.csv",outDir,anaDate.c_str());
   std::string outpath = outPath;  

   // load PP,TRLY Delta-B values
   std::vector<std::string> gradName;
   gradName.push_back("rad");
   gradName.push_back("vert");
   // gradName.push_back("azi");

   std::vector<std::string> axis; 
   axis.push_back("r"); 
   axis.push_back("y"); 
   // axis.push_back("z"); 

   const int NG = gradName.size(); 

   char inpath[500];  
   std::vector<deltab_t> pp,trly;
   for(int i=0;i<NG;i++){
      sprintf(inpath,"%s/dB-pp_final-location_%s-grad_%s.csv",outDir,gradName[i].c_str(),anaDate.c_str());
      LoadDeltaBData(inpath,pp);
      sprintf(inpath,"%s/dB-trly_final-location_%s-grad_pr-%02d_%s.csv",outDir,gradName[i].c_str(),probe,anaDate.c_str());
      LoadDeltaBData(inpath,trly);
   }

   // gather trolley probe coordinates
   char inpath_pos[500]; 
   sprintf(inpath_pos,"%s/trly-pos.csv",outDir);
   trolleyProbePosition_t trlyProbePos; 
   rc = LoadTrolleyPositionData(inpath_pos,trlyProbePos);

   double rad_coord  = trlyProbePos.r[probe]; 
   double vert_coord = trlyProbePos.y[probe]; 
   double azi_coord  = trlyProbePos.phi[probe]; 

   // load shimmed gradients.  We took NRUNS, so we'll average over them 
   std::vector<grad_meas_t> shim_grad;
   rc = GetShimmedGradients(outDir,anaDate,probe,run,gradName,shim_grad); 

   // load imposed gradients 
   char inpath_igrad[500];
   sprintf(inpath_igrad,"./input/gradients/rad-grad_vs_height.csv");
   std::vector<imposed_gradient_t> iRadGrad_vs_H;
   rc = LoadImposedGradientData(inpath_igrad,iRadGrad_vs_H); 
   sprintf(inpath_igrad,"./input/gradients/vert-grad_vs_radius.csv");
   std::vector<imposed_gradient_t> iVertGrad_vs_R;
   rc = LoadImposedGradientData(inpath_igrad,iVertGrad_vs_R); 

   TGraphErrors *gRadGrad  = GetTGraphErrors(iRadGrad_vs_H); 
   TGraphErrors *gVertGrad = GetTGraphErrors(iVertGrad_vs_R);

   gm2fieldUtil::Graph::SetGraphParameters(gRadGrad ,20,kBlack);
   gm2fieldUtil::Graph::SetGraphParameters(gVertGrad,20,kBlack);

   // fit these plots 
   TCanvas *c1 = new TCanvas("c1","Imposed Gradients",1200,800);
   c1->Divide(1,2);

   c1->cd(1); 
   gRadGrad->Draw("ap"); 
   gm2fieldUtil::Graph::SetGraphLabels(gRadGrad,"Radial Gradient vs Height Above Midplane","Height (mm)","Radial Gradient (Hz/mm)");
   gRadGrad->Draw("ap");
   gRadGrad->Fit(fitFunc.c_str()); 
   c1->Update(); 

   c1->cd(2); 
   gVertGrad->Draw("ap"); 
   gm2fieldUtil::Graph::SetGraphLabels(gVertGrad,"Vertical Gradient vs Radius","Radius (mm)","Vertical Gradient (Hz/mm)");
   gVertGrad->Draw("ap");
   gVertGrad->Fit(fitFunc.c_str()); 
   c1->Update();

   TString plotPath = Form("%s/imposed-gradients.png",plotDir);
   c1->cd();
   c1->Print(plotPath); 

   TF1 *fitRadGrad  = gRadGrad->GetFunction(fitFunc.c_str()); 
   TF1 *fitVertGrad = gVertGrad->GetFunction(fitFunc.c_str());

   // evaluate fit at probe coordinate to get gradients (Hz/mm)  
   double rad_grad  = fitRadGrad->Eval(vert_coord);
   double vert_grad = fitVertGrad->Eval(rad_coord);

   std::vector<double> imposedGrad; 
   imposedGrad.push_back(rad_grad);
   imposedGrad.push_back(vert_grad);

   // compute misalignment error in mm and in Hz 
   double arg=0,ddB=0;
   const int NPOS = gradName.size();  
   std::vector<double> dq,dq_fxpr;   // misalignment in mm 
   std::vector<double> ERR,ERR_fxpr; // error in Hz  
   for(int i=0;i<NPOS;i++){
      // raw 
      ddB = pp[i].dB - trly[i].dB;   // difference in Delta-B
      arg = ddB/imposedGrad[i];      // dividing by imposedGrad gives error in mm 
      dq.push_back(arg);
      // ABA/fxpr drift corrected 
      ddB = pp[i].dB_fxpr - trly[i].dB_fxpr; 
      arg = ddB/imposedGrad[i]; 
      dq_fxpr.push_back(arg);
      // now compute the error in Hz using the shimmed field gradients 
      arg = shim_grad[i].grad*dq[i]; 
      ERR.push_back(arg);  
      arg = shim_grad[i].grad_fxpr*dq_fxpr[i]; 
      ERR_fxpr.push_back(arg);  
   } 

   std::cout << "------------------------- RESULTS -------------------------" << std::endl;
   std::cout << Form("Misalignment error for trly probe %02d: ",probe) << std::endl;
   for(int i=0;i<NPOS;i++){
      std::cout << Form("axis: %s",axis[i].c_str()) << std::endl;
      std::cout << Form("PP Delta-B") << std::endl;
      std::cout << Form("raw = %.3lf +/- %.3lf Hz, ABA = %.3lf +/- %.3lf Hz",
                        pp[i].dB,pp[i].dB_err,pp[i].dB_fxpr,pp[i].dB_fxpr_err) << std::endl;
      std::cout << Form("TRLY Delta-B") << std::endl;
      std::cout << Form("raw = %.3lf +/- %.3lf Hz, ABA = %.3lf +/- %.3lf Hz",
                        trly[i].dB,trly[i].dB_err,trly[i].dB_fxpr,trly[i].dB_fxpr_err) << std::endl;
      std::cout << Form("imposed gradient = %.3lf Hz/mm, shimmed gradient = %.3lf +/- %.3lf Hz/mm",imposedGrad[i],shim_grad[i].grad,shim_grad[i].grad_err) << std::endl; 
      std::cout << Form("misalignment: %.3lf mm (%.3lf Hz), drift-cor: %.3lf mm (%.3lf Hz)",
                        dq[i],ERR[i],dq_fxpr[i],ERR_fxpr[i]) << std::endl;
      std::cout << "-------------------------------------------------------------" << std::endl;
   }

   rc = PrintToFile(outPath,axis,dq,ERR,dq_fxpr,ERR_fxpr);

   return 0;
}
//______________________________________________________________________________
int GetShimmedGradients(const char *prefix,std::string date,int probe,std::vector<int> run,
                        std::vector<std::string> gradName,std::vector<grad_meas_t> &shim_grad_avg){
   // load in shimmed gradients
   // we do this for N runs!
 
   const int NAXES = gradName.size(); 
   const int NRUNS = run.size();

   std::vector<double> grad,grad_fxpr; 
   double mean=0,stdev=0;
   double mean_fxpr=0,stdev_fxpr=0;
 
   char inpath[500]; 
   grad_meas_t meas; 
   std::vector<grad_meas_t> shim_grad; 
   for(int i=0;i<NAXES;i++){
      for(int j=0;j<NRUNS;j++){
	 sprintf(inpath,"%s/%s-grad_final-location_pr-%02d_run-%05d_%s.csv",prefix,gradName[i].c_str(),probe,run[j],date.c_str());
	 LoadGradientData(inpath,shim_grad);
         grad.push_back( shim_grad[i].grad ); 
         grad_fxpr.push_back( shim_grad[i].grad_fxpr ); 
         // set up for next run 
         shim_grad.clear(); 
      }
      // get stats
      mean  = gm2fieldUtil::Math::GetMean<double>(grad);
      stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(grad);    
      // FXPR drift corrected  
      mean_fxpr  = gm2fieldUtil::Math::GetMean<double>(grad_fxpr);
      stdev_fxpr = gm2fieldUtil::Math::GetStandardDeviation<double>(grad_fxpr);   
      // populate shim_grad_avg object 
      meas.name          = gradName[i]; 
      meas.grad          = mean; 
      meas.grad_err      = stdev; 
      meas.grad_fxpr     = mean_fxpr; 
      meas.grad_fxpr_err = stdev_fxpr;
      // push back into vector 
      shim_grad_avg.push_back(meas); 
      // clean up for next axis 
      grad.clear();
      grad_fxpr.clear();
   }
   return 0;
}

