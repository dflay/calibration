// Calculate the imposed gradient   

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

#include "gm2fieldMath.h"
#include "gm2fieldGraph.h"
#include "gm2fieldRootHelper.h"
#include "gm2fieldUnits.h"
#include "TemperatureSensor.h"
#include "MovingAverage.h"

#include "./include/Constants.h"
#include "./include/fixedProbeEvent.h"
#include "./include/trolleyAnaEvent.h"
#include "./include/sccEvent.h" 
#include "./include/date.h"

#include "./src/CustomUtilities.C"
#include "./src/CustomMath.C"
#include "./src/CustomImport.C"
#include "./src/CustomExport.C"
#include "./src/CustomAlgorithms.C"
#include "./src/CustomGraph.C"
#include "./src/OscFuncs.C"
#include "./src/FitFuncs.C"
#include "./src/TRLYFuncs.C"

int GetData(char axis,int i,std::vector< std::vector<double> > dB,std::vector< std::vector<double> > dB_err, 
            std::vector<double> &v,std::vector<double> &w,std::vector<double> &ew); 

int GetSlope(char axis,std::vector<int> probe,std::vector<double> dB,std::vector<double> dB_err,double &slope,double &err); 
int GetStats(int probe,int nev,std::vector<double> time,std::vector<trolleyAnaEvent_t> Data,
             std::vector<double> &MEAN,std::vector<double> &STDEV);

int CalculateDiff(std::vector<double> mu1,std::vector<double> sig1,
                  std::vector<double> mu2,std::vector<double> sig2,
                  std::vector<double> &DIFF,std::vector<double> &ERR);

int ImposedGrad_xy_prod(){
   
   int M=0;  

   int runPeriod=0;
   std::cout << "Enter run period: ";
   std::cin  >> runPeriod;

   json jConfig;
   char inpath[200];
   sprintf(inpath,"./input/json/run-%d/config.json",runPeriod); 
   std::string inpath_str = inpath; 
   int rc = gm2fieldUtil::Import::ImportJSON(inpath_str,jConfig); 

   std::string blindLabel  = jConfig["blinding"]["label"]; 
   std::string prodVersion = jConfig["prod-tag"]; 
   int nev                 = (int)jConfig["num-events-to-avg"]; 
   int fxprSet             = (int)jConfig["fxpr-set"]; 
   bool useOscCor          = (bool)( (int)jConfig["use-osc-cor"] );
   bool isBlind            = (bool)( (int)jConfig["blinding"]["enable"] ); 

   if(useOscCor) std::cout << "[ImposedGrad]: Using oscillation correction" << std::endl; 

   date_t theDate; 
   GetDate(theDate); 

   std::string outDir  = GetPath("output",isBlind,blindLabel,theDate.getDateString());

   std::vector<double> db,db_err;
   std::vector< std::vector<double> > dB,dB_err;  

   std::vector<std::string> gradName;
   gradName.push_back("rad"); 
   gradName.push_back("vert"); 

   const int NG = gradName.size(); 

   std::vector<deltab_t> trly;
   for(int i=0;i<17;i++){
      for(int j=0;j<NG;j++){
	 sprintf(inpath,"%s/dB-trly_final-location_%s-grad_pr-%02d.csv",outDir.c_str(),gradName[j].c_str(),i+1);
	 LoadDeltaBData(inpath,trly);
	 db.push_back( trly[j].dB_fxpr ); 
	 db_err.push_back( trly[j].dB_fxpr_err ); 
      }
      dB.push_back(db); 
      dB_err.push_back(db_err);
      std::cout << Form("probe %02d: dBx = %.3lf +/- %.3lf, dBy = %.3lf +/- %.3lf",
                        i+1,dB[i][0],dB_err[i][0],dB[i][1],dB_err[i][1]) << std::endl;
      // clean up
      db.clear(); 
      db_err.clear(); 
      trly.clear();
   }

   // use all probes in a specific height or radius
   int MM=0;
   std::vector<double> X,Y,EY; 
   double intercept=0,r=0,sum_sq=0,SLOPE=0;
   double slope[5] = {0,0,0,0,0}; 
   double err[5]   = {0,0,0,0,0};
   double pos[5] = {30.3,17.5,0,-17.5,-30.3};
 
   std::vector<double> VER,RG,ERG; 
   std::cout << "RADIAL GRADIENT vs HEIGHT" << std::endl;
   for(int i=0;i<5;i++){
      GetData('x',i+1,dB,dB_err,X,Y,EY); 
      rc = gm2fieldUtil::Math::LeastSquaresFitting(X,Y,intercept,SLOPE,r);
      slope[i] = SLOPE; 
      MM = X.size(); 
      for(int j=0;j<MM;j++) sum_sq += EY[j]*EY[j]; 
      err[i] = sqrt(sum_sq)/fabs(X[0]-X[MM-1]); // error estimate 
      std::cout << Form("pos = %.3lf mm, grad = %.3lf Hz/mm, err = %.3lf Hz/mm",pos[i],slope[i],err[i]) << std::endl; 
      VER.push_back(pos[i]);  
      RG.push_back(slope[i]);  
      ERG.push_back(err[i]);  
      // clean up
      X.clear();
      Y.clear();
      EY.clear();
   } 

   std::vector<double> RAD,VG,EVG; 
   std::cout << "VERTICAL GRADIENT vs RADIUS" << std::endl;
   for(int i=0;i<5;i++){
      GetData('y',i+1,dB,dB_err,X,Y,EY); 
      rc = gm2fieldUtil::Math::LeastSquaresFitting(X,Y,intercept,SLOPE,r);
      slope[i] = SLOPE; 
      MM = X.size(); 
      for(int j=0;j<MM;j++) sum_sq += EY[j]*EY[j]; 
      err[i] = sqrt(sum_sq)/fabs(X[0]-X[MM-1]); // error estimate 
      std::cout << Form("pos = %.3lf mm, grad = %.3lf Hz/mm, err = %.3lf Hz/mm",pos[i],slope[i],err[i]) << std::endl;
      RAD.push_back(pos[i]);  
      VG.push_back(slope[i]);  
      EVG.push_back(err[i]);  
      // clean up
      X.clear();
      Y.clear();
      EY.clear();
   } 

   TGraphErrors *gRadGrad_vs_Height  = gm2fieldUtil::Graph::GetTGraphErrors(VER,RG,ERG);
   gm2fieldUtil::Graph::SetGraphParameters(gRadGrad_vs_Height,20,kBlack); 

   TGraphErrors *gVertGrad_vs_Radius = gm2fieldUtil::Graph::GetTGraphErrors(RAD,VG,EVG);
   gm2fieldUtil::Graph::SetGraphParameters(gVertGrad_vs_Radius,20,kBlack); 

   TString TitleRvH = Form("Radial Gradient vs Height");
   TString TitleVvR = Form("Vertical vs Radius");

   TCanvas *c1 = new TCanvas("c1","Transverse Gradients",1200,600);
   c1->Divide(1,2);
 
   c1->cd(1); 
   gRadGrad_vs_Height->Draw("alp");
   gm2fieldUtil::Graph::SetGraphLabels(gRadGrad_vs_Height,TitleRvH,"Height Above Midplane (mm)","Gradient (Hz/mm)"); 
   gm2fieldUtil::Graph::SetGraphLabelSizes(gRadGrad_vs_Height,0.05,0.06); 
   gRadGrad_vs_Height->Draw("alp");
   c1->Update(); 

   c1->cd(2); 
   gVertGrad_vs_Radius->Draw("alp");
   gm2fieldUtil::Graph::SetGraphLabels(gVertGrad_vs_Radius,TitleVvR,"Radius (mm)","Gradient (Hz/mm)"); 
   gm2fieldUtil::Graph::SetGraphLabelSizes(gVertGrad_vs_Radius,0.05,0.06); 
   gVertGrad_vs_Radius->Draw("alp");
   c1->Update(); 

   // print to file 

   return 0;
}
//______________________________________________________________________________
int GetData(char axis,int i,std::vector< std::vector<double> > dB,std::vector< std::vector<double> > dB_err, 
            std::vector<double> &v,std::vector<double> &w,std::vector<double> &ew){

   // axis = x, y 
   // i    = coordinate index

   std::vector<int> PR; 
   if(axis=='x'){
      if(i==1){
	 // highest vertical spot
	 PR.push_back(13);
	 PR.push_back(11);
      }else if(i==2){ 
	 // second highest spot 
	 PR.push_back(14);
	 PR.push_back(4);
	 PR.push_back(10);
      }else if(i==3){
	 // center line 
	 PR.push_back(15);
	 PR.push_back(5);
	 PR.push_back(1);
	 PR.push_back(3);
	 PR.push_back(9);
      }else if(i==4){
	 // second lowest 
	 PR.push_back(16);
	 PR.push_back(2);
	 PR.push_back(8);
      }else if(i==5){
	 // lowest spot 
	 PR.push_back(17);
	 PR.push_back(7);
      }
   }else if(axis=='y'){
      if(i==1){
	 // largest radius
	 PR.push_back(8); 
	 PR.push_back(10);
      }else if(i==2){ 
	 // second highest  
	 PR.push_back(7);
	 PR.push_back(3);
	 PR.push_back(11);
      }else if(i==3){
	 // center line 
	 PR.push_back(6);
	 PR.push_back(2);
	 PR.push_back(1);
	 PR.push_back(4);
	 PR.push_back(12);
      }else if(i==4){
	 // second lowest 
	 PR.push_back(17);
	 PR.push_back(5);
	 PR.push_back(13);
      }else if(i==5){
	 // lowest radius 
	 PR.push_back(16);
	 PR.push_back(14);
      }
   }

   trolleyProbePosition_t trlyPos;
   int rc = GetTrolleyProbePositions(trlyPos);
 
   int pr=0; 
   const int NN = PR.size();
   for(int i=0;i<NN;i++){
      pr = PR[i]-1; 
      if(axis=='x'){
	 v.push_back( trlyPos.r[pr] ); 
	 w.push_back( dB[pr][0] ); 
	 ew.push_back( dB_err[pr][0] ); 
      }else if(axis=='y'){
	 v.push_back( trlyPos.y[pr] ); 
	 w.push_back( dB[pr][1] ); 
	 ew.push_back( dB_err[pr][1] ); 
      }
   } 

   return 0; 

}
//______________________________________________________________________________
int GetSlope(char axis,std::vector<int> probe,std::vector<double> dB,std::vector<double> dB_err,double &slope,double &err){
   // get the slope across the dB values 
   // choose the probes according to the input list 
   const int N = dB.size();
   const int M = probe.size(); 

   trolleyProbePosition_t trlyPos; 
   int rc = GetTrolleyProbePositions(trlyPos); 

   double thePos=0;
   std::vector<double> pos,v,ev;
   for(int i=0;i<M;i++){
      for(int j=0;j<N;j++){
	 if(probe[i]==j){ // found the probe of interest 
	    if(axis=='x') thePos = trlyPos.r[probe[i]];   
	    if(axis=='y') thePos = trlyPos.y[probe[i]];  
            pos.push_back(thePos);
	    v.push_back(dB[j]); 
	    ev.push_back(dB_err[j]);
	    std::cout << Form("MATCH probe %d, pos = %.3lf mm, dB = %.3lf +/- %.3lf Hz",probe[i],thePos,dB[j],dB_err[j]) << std::endl; 
	 }
      }
   }

   int MM = ev.size();
   if(MM!=M){
      std::cout << "Data set doesn't match expectation!  M = " << M << " MM = " << MM << std::endl;
      return 1;
   }

   double intercept=0,r=0;
   rc = gm2fieldUtil::Math::LeastSquaresFitting(pos,v,intercept,slope,r); 
  
   double sum_sq=0;
   for(int i=0;i<MM;i++) sum_sq += ev[i]*ev[i]; 
   err = sqrt(sum_sq)/fabs(pos[0]-pos[MM-1]); // error estimate 

   return 0; 

}
//______________________________________________________________________________
int CalculateDiff(std::vector<double> mu1,std::vector<double> sig1,
                  std::vector<double> mu2,std::vector<double> sig2,
                  std::vector<double> &DIFF,std::vector<double> &ERR){
   // compute shift in field
   // mu1,sig1 = mean and stdev for SCC on 
   // mu2,sig2 = mean and stdev for SCC off  
   const int N1 = mu1.size();
   const int N2 = mu2.size();
   if(N1!=N2){
      std::cout << "[CalculateDiff]: Error!  Vector sizes don't match! " << std::endl;
      return 1;
   }
   double diff=0,err=0;
   for(int i=0;i<N1;i++){
      diff = mu1[i] - mu2[i];
      err  = TMath::Sqrt( sig1[i]*sig1[i] + sig2[i]*sig2[i] );
      DIFF.push_back(diff); 
      ERR.push_back(err); 
   }
   return 0;
} 
//______________________________________________________________________________
int GetStats(int probe,int nev,std::vector<double> time,std::vector<trolleyAnaEvent_t> Data,
             std::vector<double> &MEAN,std::vector<double> &STDEV){

   const int N = time.size();
   int M=0,rc=0;
   double mean=0,stdev=0;
   std::vector<double> freq;
   std::vector<trolleyAnaEvent_t> Event; 
   for(int i=0;i<N;i++){
      // std::cout << "Looking for time " << gm2fieldUtil::GetStringTimeStampFromUTC(time[i]) << std::endl;
      // find events 
      rc = FilterSingle("freq",probe,nev,time[i],Data,freq); 
      // now get mean of events 
      mean  = gm2fieldUtil::Math::GetMean<double>(freq); 
      stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(freq);  
      // store result
      MEAN.push_back(mean); 
      STDEV.push_back(stdev); 
      // set up for next time 
      freq.clear(); 
   }
  
   // now restrict to less than 20 results if necessary 
   int cntr = MEAN.size();
   while(cntr>20){
      MEAN.pop_back();
      STDEV.pop_back();
      cntr = MEAN.size();
   }

   return 0;
}

