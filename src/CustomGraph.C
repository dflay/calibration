#include "../include/CustomGraph.h"
//______________________________________________________________________________
TGraphErrors *GetTGraphErrors(std::vector<imposed_gradient_t> data){
   std::vector<double> x,y,ey;
   const int N = data.size();
   for(int i=0;i<N;i++){
      x.push_back( data[i].pos ); 
      y.push_back( data[i].grad ); 
      ey.push_back( data[i].grad_err ); 
   }
   TGraphErrors *g = gm2fieldUtil::Graph::GetTGraphErrors(x,y,ey); 
   return g;
}
//______________________________________________________________________________
TGraph *GetPPTGraph1(TString xAxis,TString yAxis,std::vector<plungingProbeAnaEvent_t> data){
   // a custom plotter for the plunging probe analysis event
   // treat all sub-events as a regular event 
   std::vector<double> x,y; 
   FillPPVector1(xAxis,data,x); 
   FillPPVector1(yAxis,data,y); 
   TGraph *g = gm2fieldUtil::Graph::GetTGraph(x,y); 
   return g;
}
//______________________________________________________________________________
TGraph *GetPPTGraph2(TString xAxis,TString yAxis,std::vector<plungingProbeAnaEvent_t> data){
   // a custom plotter for the plunging probe analysis event
   // average over all sub-events  
   std::vector<double> x,y; 
   FillPPVector2(xAxis,data,x); 
   FillPPVector2(yAxis,data,y); 
   TGraph *g = gm2fieldUtil::Graph::GetTGraph(x,y); 
   return g;
}
//______________________________________________________________________________
TGraph *GetPPTGraph3(TString xAxis,TString yAxis,plungingProbeAnaEvent_t data){
   // a custom plotter for the plunging probe analysis event
   // these data have a single NMR-DAQ run consisting of an arbitrary number of sub-events  
   std::vector<double> x,y; 
   FillPPVector3(xAxis,data,x); 
   FillPPVector3(yAxis,data,y); 
   TGraph *g = gm2fieldUtil::Graph::GetTGraph(x,y); 
   return g;
}
//______________________________________________________________________________
TGraphErrors *GetPPTGraphErrors2(TString xAxis,TString yAxis,std::vector<plungingProbeAnaEvent_t> data){
   // a custom plotter for the plunging probe analysis event
   TString yAxis_err = Form("%s_err",yAxis.Data()); 
   std::vector<double> x,y,ey; 
   FillPPVector2(xAxis,data,x); 
   FillPPVector2(yAxis,data,y); 
   FillPPVector2(yAxis_err,data,ey);
   // const int N = x.size(); 
   // std::cout << N << " data points" << std::endl;
   // for(int i=0;i<N;i++){
   //    std::cout << Form("***** F = %.3lf +/- %.3lf Hz",y[i],ey[i]) << std::endl;
   // }
   // std::cout << "---------------------------" << std::endl; 
   TGraphErrors *g = gm2fieldUtil::Graph::GetTGraphErrors(x,y,ey); 
   return g;
}
//______________________________________________________________________________
TGraphErrors *GetPPScanGraph(TString xAxis,TString yAxis,std::vector<plungingProbeAnaEvent_t> data,double x0){
   // a custom plotter for the plunging probe analysis event
   TString yAxis_err = Form("%s_err",yAxis.Data()); 
   std::vector<double> x,y,ey; 
   FillPPVector2(xAxis,data,x); 
   FillPPVector2(yAxis,data,y); 
   FillPPVector2(yAxis_err,data,ey);

   // subtract off the offset 
   std::vector<double> X; 
   const int N = x.size();
   for(int i=0;i<N;i++) X.push_back(x[i]-x0);

   TGraphErrors *g = gm2fieldUtil::Graph::GetTGraphErrors(X,y,ey); 

   return g;
}
// //______________________________________________________________________________
// TGraph2D *GetAzimuthalProjection(std::vector<plungingProbeAnaEvent_t> data,int units){
//    double REF = 50E+3; // reference mixdown frequency (right way to do it...)  
//    double arg=0;
//    std::vector<double> x,y,z;
//    const int N = data.size();
//    for(int i=0;i<N;i++){
//       // get coordinates of probe 
//       x.push_back(data[i].r);  
//       y.push_back(data[i].y); 
//       arg = data[i].freq;
//       if(units==gm2fieldUtil::Units::kHz) arg /= 1E+3; 
//       if(units==gm2fieldUtil::Units::ppm) arg /= 61.79; 
//       if(units==gm2fieldUtil::Units::ppb) arg /= 0.06179; 
//       if( gm2fieldUtil::Math::IsInfOrNaN<double>(arg) ) std::cout << "ERROR: " << arg << std::endl;
//       z.push_back(arg);
//       // set up for next probe  
//    }
//    TGraph2D *g = new TGraph2D(N,&x[0],&y[0],&z[0]);
//    return g;
// }
//______________________________________________________________________________
int FillPPVector1(TString axis,std::vector<plungingProbeAnaEvent_t> data,std::vector<double> &x){
   // get all sub-events  
   int M=0;
   double mean=0;
   std::vector<double> v; 
   const int N = data.size();
   for(int i=0;i<N;i++){
      M = data[i].numTraces;
      if(axis=="run")        for(int j=0;j<M;j++) x.push_back( (double)data[i].run );   
      if(axis=="TimeStamp")  for(int j=0;j<M;j++) x.push_back(data[i].time[j]/1E+9 );   
      if(axis=="x")          for(int j=0;j<M;j++) x.push_back(data[i].r[j]         ); 
      if(axis=="y")          for(int j=0;j<M;j++) x.push_back(data[i].y[j]         );  
      if(axis=="z")          for(int j=0;j<M;j++) x.push_back(data[i].phi[j]       );  
      if(axis=="temp")       for(int j=0;j<M;j++) x.push_back(data[i].temp[j]      );  
      if(axis=="temp_err")   for(int j=0;j<M;j++) x.push_back(data[i].temp_err[j]  );  
      if(axis=="freq")       for(int j=0;j<M;j++) x.push_back(data[i].freq[j]      );  
      if(axis=="freq_err")   for(int j=0;j<M;j++) x.push_back(data[i].freq_err[j]  );  
      if(axis=="freq_LO")    for(int j=0;j<M;j++) x.push_back(data[i].freq_LO[j]   );  
      if(axis=="freq_RF")    for(int j=0;j<M;j++) x.push_back(data[i].freq_RF[j]   ); 
   }            
   return 0; 
}
//______________________________________________________________________________
int FillPPVector2(TString axis,std::vector<plungingProbeAnaEvent_t> data,std::vector<double> &x){
   // average over all traces in an event
   int M=0; 
   double mean=0,stdev=0;
   std::vector<double> w; 
   const int N = data.size();
   if(axis=="run"){
      for(int i=0;i<N;i++) x.push_back( (double)data[i].run );   
   }else if(axis=="TimeStamp"){
      for(int i=0;i<N;i++) x.push_back( data[i].time[0]/1E+9  );  // choose first time value here 
   }else{
      for(int i=0;i<N;i++){
	 M = data[i].numTraces;
	 if(axis=="x")        for(int j=0;j<M;j++) w.push_back(data[i].r[j]      ); 
	 if(axis=="y")        for(int j=0;j<M;j++) w.push_back(data[i].y[j]      );  
	 if(axis=="z")        for(int j=0;j<M;j++) w.push_back(data[i].phi[j]    );  
	 if(axis=="temp")     for(int j=0;j<M;j++) w.push_back(data[i].temp[j]   );  
	 if(axis=="temp_err") for(int j=0;j<M;j++) w.push_back(data[i].temp[j]   );  
	 if(axis=="freq")     for(int j=0;j<M;j++) w.push_back(data[i].freq[j]   );  
	 if(axis=="freq_err") for(int j=0;j<M;j++) w.push_back(data[i].freq[j]   );  
	 if(axis=="freq_LO")  for(int j=0;j<M;j++) w.push_back(data[i].freq_LO[j]);  
	 if(axis=="freq_RF")  for(int j=0;j<M;j++) w.push_back(data[i].freq_RF[j]); 
	 mean  = gm2fieldUtil::Math::GetMean<double>(w); 
	 stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(w); 
	 if(axis=="temp_err" || axis=="freq_err"){
	    x.push_back(stdev);
         }else{
	    x.push_back(mean);
         }
	 // clean up for next event 
	 w.clear();
      }
   }            
   return 0; 
}
//______________________________________________________________________________
int FillPPVector3(TString axis,plungingProbeAnaEvent_t data,std::vector<double> &x){
   // grab all events 
   int M=0; 
   double mean=0,stdev=0;
   std::vector<double> w; 
   const int N = data.numTraces;
   for(int i=0;i<N;i++){
      if(axis=="run")       x.push_back( (double)data.run );  
      if(axis=="TimeStamp") x.push_back(data.time[i]/1E+9); 
      if(axis=="x")         x.push_back(data.r[i]      ); 
      if(axis=="y")         x.push_back(data.y[i]      );  
      if(axis=="z")         x.push_back(data.phi[i]    );  
      if(axis=="temp")      x.push_back(data.temp[i]   );  
      if(axis=="temp_err")  x.push_back(data.temp[i]   );  
      if(axis=="freq")      x.push_back(data.freq[i]   );  
      if(axis=="freq_err")  x.push_back(data.freq[i]   );  
      if(axis=="freq_LO")   x.push_back(data.freq_LO[i]);  
      if(axis=="freq_RF")   x.push_back(data.freq_RF[i]); 
   }  

   if(axis=="freq"){
      for(int i=0;i<N;i++) if(x[i]>20E+3) std::cout << x[i] << std::endl;
   }
          
   return 0; 
}
//______________________________________________________________________________
TGraphErrors *GetTGraphErrors(int method,unsigned long long timeStart,unsigned long long timeStop,unsigned long long timeStep,
                              std::vector<int> fxprList,std::vector<gm2field::fixedProbeFrequency_t> fxprData){
   // a plot between timeStart and timeStop of the average field from probes in the list fxprList
   // the user specifies the step size in time 
   double err=30.*0.06179;  // assume a 30 ppb error 
   unsigned long long t0 = 0; 
   std::vector<unsigned long long> time;
   std::vector<double> TIME,FREQ,ERR; 
   GetAverageFXPRVectors(method,t0,timeStart,timeStop,timeStep,fxprList,fxprData,time,FREQ); 
   int NPTS = time.size();
   for(int i=0;i<NPTS;i++){
      TIME.push_back(time[i]/1E+9);
      ERR.push_back(err);
   }
   TGraphErrors *g = gm2fieldUtil::Graph::GetTGraphErrors(TIME,FREQ,ERR);
   return g; 
}
//______________________________________________________________________________
TGraph *GetTGraphNew(std::vector<fixedProbeEvent_t> fxprData){
   std::vector<double> TIME,FREQ; 
   int NPTS = fxprData.size();
   for(int i=0;i<NPTS;i++){
      TIME.push_back(fxprData[i].time);
      FREQ.push_back(fxprData[i].freq);
   }
   TGraph *g = gm2fieldUtil::Graph::GetTGraph(TIME,FREQ);
   return g; 
}
//______________________________________________________________________________
TGraph *GetTGraph(int method,unsigned long long timeStart,unsigned long long timeStop,unsigned long long timeStep,
                  std::vector<int> fxprList,std::vector<gm2field::fixedProbeFrequency_t> fxprData){
   // a plot between timeStart and timeStop of the average field from probes in the list fxprList
   // the user specifies the step size in time
   unsigned long long t0 = 0; 
   std::vector<unsigned long long> time;
   std::vector<double> TIME,FREQ; 
   std::cout << "Getting vectors..." << std::endl; 
   GetAverageFXPRVectors(method,t0,timeStart,timeStop,timeStep,fxprList,fxprData,time,FREQ);
   std::cout << "--> Done." << std::endl; 
   int NPTS = time.size();
   std::cout << "Processing " << NPTS << " events..." << std::endl;
   for(int i=0;i<NPTS;i++) TIME.push_back(time[i]/1E+9);
   std::cout << "--> Done." << std::endl; 
   std::cout << "Getting graph..." << std::endl;
   TGraph *g = gm2fieldUtil::Graph::GetTGraph(TIME,FREQ);
   std::cout << "--> Done." << std::endl; 
   return g; 
}
//______________________________________________________________________________
TGraph *GetTGraph(int method,unsigned long long t0,unsigned long long timeStart,unsigned long long timeStop,unsigned long long timeStep,
                  std::vector<int> fxprList,std::vector<gm2field::fixedProbeFrequency_t> fxprData){
   // a plot between timeStart and timeStop of the average field from probes in the list fxprList
   // the user specifies the step size in time 
   std::vector<unsigned long long> time;
   std::vector<double> TIME,FREQ; 
   GetAverageFXPRVectors(method,t0,timeStart,timeStop,timeStep,fxprList,fxprData,time,FREQ); 
   int NPTS = time.size();
   for(int i=0;i<NPTS;i++) TIME.push_back(time[i]/1E+9);
   TGraph *g = gm2fieldUtil::Graph::GetTGraph(TIME,FREQ);
   return g; 
}
//______________________________________________________________________________
TGraph *GetTGraph2Runs(int method,std::vector<int> probe,
                       unsigned long long tStart,unsigned long long tStop,unsigned long long tStep,
                       std::vector<gm2field::fixedProbeFrequency_t> fxprData){

   // determine average FXPR frequency across a list of probes,  
   // using a range defined by tStart < t < tStop with a step size tStep 
   // the times tStart, tStop are necessarily end/start points of runs that 
   // are far in time (more than a second or two)   

   // populate all the data we have in fxprData into vectors 
   std::vector<unsigned long long> time;
   std::vector<double> TIME,FREQ; 
   GetAverageFXPRVectorsAlt(method,tStart,tStop,tStep,probe,fxprData,time,FREQ); 
   int NPTS = time.size();
   for(int i=0;i<NPTS;i++) TIME.push_back(time[i]/1E+9);

   // now get a graph 
   TGraph *g = gm2fieldUtil::Graph::GetTGraph(TIME,FREQ);
   return g;
}
//______________________________________________________________________________
TGraph *GetInterpolatedTGraph(int method,std::vector<int> probe,
                              unsigned long long tStart,unsigned long long tStop,unsigned long long tStep,
                              std::vector<gm2field::fixedProbeFrequency_t> fxprData,
                              std::vector<double> &stats){

   // determine average FXPR frequency across a list of probes,  
   // using a range defined by tStart < t < tStop with a step size tStep 
   // the times tStart, tStop are necessarily end/start points of runs that 
   // are far in time (more than a second or two)   

   // populate all the data we have in fxprData into vectors 
   std::vector<unsigned long long> time;
   std::vector<double> TIME,FREQ; 
   GetAverageFXPRVectorsAlt(method,tStart,tStop,tStep,probe,fxprData,time,FREQ); 
   int NPTS = time.size();
   for(int i=0;i<NPTS;i++) TIME.push_back(time[i]/1E+9);

   // do a very crude linear fit
   // take average over previous ten seconds to get the starting field value 
   std::vector<double> tt,ff; 
   double aFreq,stdev;
   unsigned long long aTime = 0;
   for(int i=10;i>=1;i--){
      aTime = tStart - ( (double)i )*1E+9;  
      GetAverageFXPR(method,aTime,probe,fxprData,aFreq,stdev);
      ff.push_back(aFreq); 
   }
   double fStart     = gm2fieldUtil::Math::GetMean<double>(ff); 
   double fStart_err = gm2fieldUtil::Math::GetStandardDeviation<double>(ff);  
   ff.clear(); 
   // take average over first ten seconds of last run to get the stop field value 
   for(int i=1;i<=10;i++){
      aTime = tStop + ( (double)i )*1E+9;  
      GetAverageFXPR(method,aTime,probe,fxprData,aFreq,stdev);
      ff.push_back(aFreq); 
   }
   double fStop     = gm2fieldUtil::Math::GetMean<double>(ff); 
   double fStop_err = gm2fieldUtil::Math::GetStandardDeviation<double>(ff); 
   // store these data 
   stats.push_back(fStart); 
   stats.push_back(fStart_err); 
   stats.push_back(fStop); 
   stats.push_back(fStop_err); 
   // now fill in the points
   double t0 = tStart/1E+9;
   double t1 = tStop/1E+9;
   double f0 = fStart;
   double f1 = fStop;
   NPTS = (tStop-tStart)/(tStep);
   double theTime=0,theFreq=0;
   double dt = tStep/1E+9;
   for(int i=0;i<NPTS;i++){
      theTime = t0 + ( (double)i )*dt;
      theFreq = gm2fieldUtil::Math::LinearInterpolation(theTime,t0,f0,t1,f1);
      // std::cout << gm2fieldUtil::GetStringTimeStampFromUTC(theTime) << " " << Form("%.3lf",theFreq) << std::endl;
      TIME.push_back(theTime);
      FREQ.push_back(theFreq);
   }
   // now get a graph 
   TGraph *g = gm2fieldUtil::Graph::GetTGraph(TIME,FREQ);
   return g;
}
//______________________________________________________________________________
TGraph *GetTRLYTGraph(int probe,TString xAxis,TString yAxis,std::vector<trolleyAnaEvent_t> data,double sf){
   // a custom plotter for the trolley analysis event 
   std::vector<double> x,y; 
   FillTRVector(probe,xAxis,data,sf,x); 
   FillTRVector(probe,yAxis,data,sf,y); 
   TGraph *g = gm2fieldUtil::Graph::GetTGraph(x,y); 
   return g;
}
//______________________________________________________________________________
TGraphErrors *GetSlicePlot(char axis,std::vector<trolleyAnaEvent_t> trlyData){
   // get plot from trolley probes to estimate vertical and transverse gradients 
   const int NP = 5;
   int radProbe[NP]  = {14,4,0,2,8 };
   int vertProbe[NP] = { 5,1,0,3,11};
   int probe[NP];
   if(axis=='r') for(int i=0;i<NP;i++) probe[i] = radProbe[i];
   if(axis=='v') for(int i=0;i<NP;i++) probe[i] = vertProbe[i];
   int k=0;
   double mean_freq=0,stdev_freq=0;
   std::vector<double> f,X,Y,EY;
   const int N = trlyData.size();
   for(int j=0;j<NP;j++){
      k = probe[j];
      for(int i=0;i<N;i++) f.push_back( trlyData[i].freq[k] );
      // for a given probe, get the mean frequency over all events  
      mean_freq  = gm2fieldUtil::Math::GetMean<double>(f);
      stdev_freq = gm2fieldUtil::Math::GetStandardDeviation<double>(f);
      Y.push_back(mean_freq);
      EY.push_back(stdev_freq);
      // std::cout << Form("PROBE %02d DETAILS: r = %.3lf, y = %.3lf, f = %.3lf +/- %.3lf",
      //                   k,trlyData[0].r[k],trlyData[0].y[k],mean_freq,stdev_freq) << std::endl;
      if(axis=='r') X.push_back( trlyData[0].r[k] );
      if(axis=='v') X.push_back( trlyData[0].y[k] );
      // clean up for next probe 
      f.clear();
   }
   TGraphErrors *g = gm2fieldUtil::Graph::GetTGraphErrors(X,Y,EY);
   return g;
}
//______________________________________________________________________________
TGraphErrors *GetSlicePlot(char axis,std::vector<trolleyAnaEvent_t> trlyData,
                           std::vector<double> &X,std::vector<double> &Y,std::vector<double> &EY){
   // get plot from trolley probes to estimate vertical and transverse gradients 
   const int NP = 5;
   int radProbe[NP]  = {14,4,0,2,8 };
   int vertProbe[NP] = { 5,1,0,3,11};
   int probe[NP];
   if(axis=='r') for(int i=0;i<NP;i++) probe[i] = radProbe[i];
   if(axis=='v') for(int i=0;i<NP;i++) probe[i] = vertProbe[i];
   int k=0;
   double mean_freq=0,stdev_freq=0;
   std::vector<double> f;
   const int N = trlyData.size();
   for(int j=0;j<NP;j++){
      k = probe[j];
      for(int i=0;i<N;i++) f.push_back( trlyData[i].freq[k] );
      // for a given orobe, get the mean frequency over all events  
      mean_freq  = gm2fieldUtil::Math::GetMean<double>(f);
      stdev_freq = gm2fieldUtil::Math::GetStandardDeviation<double>(f);
      Y.push_back(mean_freq);
      EY.push_back(stdev_freq);
      if(axis=='r') X.push_back( trlyData[0].r[k] );
      if(axis=='v') X.push_back( trlyData[0].y[k] );
      // clean up for next probe 
      f.clear();
   }
   TGraphErrors *g = gm2fieldUtil::Graph::GetTGraphErrors(X,Y,EY);
   return g;
}
//______________________________________________________________________________
TGraph2D *GetAzimuthalProjection(std::vector<trolleyAnaEvent_t> data,int units){
   double REF = 50E+3; // reference mixdown frequency (right way to do it...)  
   double arg=0,mean=0;
   std::vector<double> x,y,z,w;
   const int M = data.size();
   for(int i=0;i<NUM_TRLY;i++){
      // get coordinates of trolley probes 
      x.push_back(data[0].r[i]);  
      y.push_back(data[0].y[i]); 
      // average over all events for each probe  
      for(int j=0;j<M;j++){
         arg = data[j].freq[i]; 
         if(units==gm2fieldUtil::Constants::kHz) arg /= 1E+3; 
         if(units==gm2fieldUtil::Constants::ppm) arg /= 61.79; 
         if(units==gm2fieldUtil::Constants::ppb) arg /= 0.06179; 
         if( gm2fieldUtil::Math::IsInfOrNaN<double>(arg) ) std::cout << "ERROR: " << arg << std::endl;
	 w.push_back(arg);
      }
      mean = gm2fieldUtil::Math::GetMean<double>(w);
      z.push_back(mean); 
      // std::cout << x[i] << " " << y[i] << " " << z[i] << std::endl;
      // set up for next probe  
      w.clear();
   }
   TGraph2D *g = new TGraph2D(NUM_TRLY,&x[0],&y[0],&z[0]);
   return g;
}
//______________________________________________________________________________
TGraphErrors *GetSCCTestGraphTRLY(int probe,TString xAxis,TString yAxis,std::vector< std::vector<sccTrlyEvent_t> > data){
   std::vector<double> x,y,ey;
   const int N = data[probe].size();
   for(int i=0;i<N;i++){
      if(xAxis=="trial") x.push_back( (double)i+1 );  
      if(xAxis=="bare")  x.push_back( (double)data[probe][i].freq_bare );  
      if(xAxis=="scc")   x.push_back( (double)data[probe][i].freq_scc  );  
      if(xAxis=="diff")  x.push_back( (double)data[probe][i].freq_diff );  
      if(yAxis=="trial") y.push_back( (double)i+1 );  
      if(yAxis=="bare")  y.push_back( (double)data[probe][i].freq_bare );  
      if(yAxis=="scc")   y.push_back( (double)data[probe][i].freq_scc  );  
      if(yAxis=="diff")  y.push_back( (double)data[probe][i].freq_diff );  
      if(yAxis=="bare")  ey.push_back( (double)data[probe][i].freq_bare_err );  
      if(yAxis=="scc")   ey.push_back( (double)data[probe][i].freq_scc_err  );  
      if(yAxis=="diff")  ey.push_back( (double)data[probe][i].freq_diff_err );  
   }
   TGraphErrors *g = gm2fieldUtil::Graph::GetTGraphErrors(x,y,ey); 
   return g; 
}
//______________________________________________________________________________
TGraph *GetTRLYPositionsGraph(){
   // get a plot of trolley positions 
   trolleyProbePosition_t trlyPos;
   int rc = GetTrolleyProbePositions(trlyPos);
   if(rc!=0) return NULL; 

   std::vector<double> x,y;
   for(int i=0;i<NUM_TRLY;i++){
      x.push_back( trlyPos.r[i]*0.99 );
      y.push_back( trlyPos.y[i]*0.99 );
   }

   TGraph *g = gm2fieldUtil::Graph::GetTGraph(x,y);
   g->SetMarkerStyle(20);
   g->SetMarkerSize(1.0);

   return g;
}
//______________________________________________________________________________
TGraph *GetSCCPlot(int type,std::vector<gm2field::surfaceCoils_t> data){
   // construct the sum of the top (type=1), bottom (type=0) or azi (type=-1) coil currents 
   // to see what the SCC config is
   double sum=0;
   std::vector<double> x,y;
   int M=4; // for azi coils
   if(type==0||type==1) M = 100;
   if( TMath::Abs(type)>1 ) return NULL;
   int NEV = data.size();
   for(int i=0;i<NEV;i++){
      for(int j=0;j<M;j++){
         if(type==0)  sum += data[i].BotCurrents[j];
         if(type==1)  sum += data[i].TopCurrents[j];
         if(type==-1) sum  = data[i].AzCurrents[0];
      }
      if(type==0)  x.push_back( data[i].BotTime[0]/1E+9 );
      if(type==1)  x.push_back( data[i].TopTime[0]/1E+9 );
      if(type==-1) x.push_back( data[i].TopTime[0]/1E+9 ); // don't have an associated Az time...  
      y.push_back(sum);
      // set up for next event 
      sum = 0;
   }
   TGraph *g = gm2fieldUtil::Graph::GetTGraph(x,y);
   return g;
}
//______________________________________________________________________________
int FillTRVector(int probe,TString axis,std::vector<trolleyAnaEvent_t> data,double sf,std::vector<double> &x){
   const int N = data.size();
   if(axis=="GpsTimeStamp") for(int i=0;i<N;i++) x.push_back( data[i].time[probe]/1E+9 );  
   if(axis=="r")            for(int i=0;i<N;i++) x.push_back( sf*data[i].r[probe]      );  
   if(axis=="y")            for(int i=0;i<N;i++) x.push_back( sf*data[i].y[probe]      );  
   if(axis=="phi")          for(int i=0;i<N;i++) x.push_back( sf*data[i].phi[probe]    );  
   if(axis=="temp")         for(int i=0;i<N;i++) x.push_back( sf*data[i].temp[probe]   );  
   if(axis=="freq")         for(int i=0;i<N;i++) x.push_back( sf*data[i].freq[probe]   );  
   return 0; 
}
//______________________________________________________________________________
TGraphErrors *GetTRLYTGraph_aziScan(int probe,double thr,TString yAxis,std::vector<trolleyAnaEvent_t> data){
   double v=0;
   double v_prev = data[0].phi[probe]; 
   double angle=0,mean=0,stdev=0;
   std::vector<double> X,Y,x,y,ey; 
   const int N = data.size();
   for(int i=1;i<N;i++){
      v = data[i].phi[probe]; 
      if( TMath::Abs(v-v_prev)<thr ){
	 // very close in angle, add to vector 
	 angle = v;
	 if(yAxis=="r")    Y.push_back( data[i].r[probe]         );  
	 if(yAxis=="y")    Y.push_back( data[i].y[probe]         );  
	 if(yAxis=="temp") Y.push_back( data[i].temp[probe]      );  
	 if(yAxis=="freq") Y.push_back( data[i].freq[probe]      );  
      }else{
	 // outside threshold, average and move on 
	 mean  = gm2fieldUtil::Math::GetMean<double>(Y); 
	 stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(Y);
         x.push_back(angle);
	 y.push_back(mean);
	 ey.push_back(stdev);
	 // clear vectors
	 Y.clear();
      }
      v_prev = v; 
   }
   // make the plot 
   TGraphErrors *g = gm2fieldUtil::Graph::GetTGraphErrors(x,y,ey);
   return g;
}

