#include "../include/CustomMath.h" 
//______________________________________________________________________________
int GetStats_vec(int lo,int hi,std::vector<double> x,double &mean,double &stdev){
   const int N = x.size();
   int start=0,stop=hi;
   if(lo>0) start = lo; 
   if(hi>N) stop  = N;  
   std::vector<double> v;
   for(int i=start;i<stop;i++) v.push_back(x[i]); 
   mean  = gm2fieldUtil::Math::GetMean<double>(v); 
   stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(v); 
   return 0;  
}
//______________________________________________________________________________
int GetStats(std::string varName,plungingProbeAnaEvent_t data,double &mean,double &stdev){
   const int N = data.numTraces;
   if(N==0) return -1;

   std::vector<double> x;
   for(int i=0;i<N;i++){
      if( varName.compare("r")==0    ) x.push_back(data.r[i]   ); 
      if( varName.compare("y")==0    ) x.push_back(data.y[i]   ); 
      if( varName.compare("phi")==0  ) x.push_back(data.phi[i] ); 
      if( varName.compare("freq")==0 ) x.push_back(data.freq[i]); 
      if( varName.compare("temp")==0 ) x.push_back(data.temp[i]); 
   }

   mean  = gm2fieldUtil::Math::GetMean<double>(x); 
   stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(x);

   return 0;   
}
//______________________________________________________________________________
int GetStats(std::string varName,std::vector<plungingProbeAnaEvent_t> data,double &mean,double &stdev){
   const int N = data.size();
   if(N==0) return -1;

   int M=0;
   std::vector<double> x;
   for(int i=0;i<N;i++){
      M = data[i].numTraces; 
      for(int j=0;j<M;j++){
	 if( varName.compare("r")==0    ) x.push_back(data[i].r[j]   ); 
	 if( varName.compare("y")==0    ) x.push_back(data[i].y[j]   ); 
	 if( varName.compare("phi")==0  ) x.push_back(data[i].phi[j] ); 
	 if( varName.compare("freq")==0 ) x.push_back(data[i].freq[j]); 
	 if( varName.compare("temp")==0 ) x.push_back(data[i].temp[j]); 
      }
   }

   mean  = gm2fieldUtil::Math::GetMean<double>(x); 
   stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(x);

   return 0;   
}
//______________________________________________________________________________
int GetStats(std::string varName,std::vector<calibSwap_t> data,double &mean,double &stdev){
   const int N = data.size();
   if(N==0) return -1;

   std::vector<double> x;
   for(int i=0;i<N;i++){
      if( varName.compare("freq")==0 ) x.push_back(data[i].freq); 
      if( varName.compare("temp")==0 ) x.push_back(data[i].temp); 
   }

   mean  = gm2fieldUtil::Math::GetMean<double>(x); 
   stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(x);

   return 0;   
}
//______________________________________________________________________________
int GetStats(std::string varName,int probe,std::vector<trolleyAnaEvent_t> data,double &mean,double &stdev){
   // get stats for the trolley data 
   const int N = data.size();
   if(N==0) return -1; 

   std::vector<double> x;
   for(int i=0;i<N;i++){
      if( varName.compare("freq")==0 ) x.push_back(data[i].freq[probe]); 
      if( varName.compare("temp")==0 ) x.push_back(data[i].temp[probe]); 
   }

   mean  = gm2fieldUtil::Math::GetMean<double>(x); 
   stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(x);

   return 0; 
}
//______________________________________________________________________________
int GetStats(std::vector<double> x,double &mean,double &stdev){
   // compute mean and standard deviation of data 
   const int N = x.size();
   mean  = gm2fieldUtil::Math::GetMean<double>(x);
   stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(x);
   // std::cout << N << " events" << std::endl; 
   // for(int i=0;i<N;i++) std::cout << x[i] << std::endl; 
   // std::cout << "MEAN = " << mean << " " << "STDEV = " << stdev << std::endl; 
   // std::cout << "------------------" << std::endl;
   return 0;
}
//______________________________________________________________________________
int GetStats(int probeNumber,int method,std::vector<gm2field::fixedProbeFrequency_t> Data,double &mean,double &stdev){
   // compute mean and standard deviation of data for a given probe number
   const int N = Data.size();
   std::vector<double> freq;
   for(int i=0;i<N;i++){
      freq.push_back(Data[i].Frequency[probeNumber][method]);
   }
   mean  = gm2fieldUtil::Math::GetMean<double>(freq);
   stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(freq);
   return 0;
}
//______________________________________________________________________________
int GetStats(int probeNumber,int method,std::vector<gm2field::trolleyProbeFrequency_t> Data,double &mean,double &stdev){
   // compute mean and standard deviation of data for a given probe number
   const int N = Data.size();
   std::vector<double> freq;
   for(int i=0;i<N;i++){
      freq.push_back(Data[i].Frequency[probeNumber][method]);
   }
   mean  = gm2fieldUtil::Math::GetMean<double>(freq);
   stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(freq);
   return 0;
}
//______________________________________________________________________________
int GetStats(TString axis,std::vector<gm2field::psFeedback_t> Data,double &mean,double &stdev){
   // compute mean and standard deviation of data for a given probe number
   const int N = Data.size();
   std::vector<double> x;
   for(int i=0;i<N;i++){
      if(axis=="filtered_mean_freq") x.push_back(Data[i].filtered_mean_freq);
      if(axis=="field_setpoint")     x.push_back(Data[i].field_setpoint);
      if(axis=="current")            x.push_back(Data[i].current);
   }
   mean  = gm2fieldUtil::Math::GetMean<double>(x);
   stdev = gm2fieldUtil::Math::GetStandardDeviation<double>(x);
   return 0;
}
//______________________________________________________________________________
double MultipoleFitFunc(double *x,double *par){
   // multipole expansion of the magnetic field 
   // i = 0: num pars   
   // i = 1: dipole
   // i = 2: quadrupole 
   // i = 3: sextupole
   // i = 4: octupole 
   // i = 5: decupole 
   const int N = par[0]; 
   double R0   = 4.5;                    // in cm 
   double X    = x[0];                   // in cm 
   double Y    = x[1];                   // in cm 
   double r    = TMath::Sqrt(X*X + Y*Y); // in cm  
   double thr  = TMath::ATan(Y/X);       // in rad

   double j=0;
   double f = par[1];  
   for(int i=2;i<N;i+=2){
      j  = (double)i/2.;
      f += TMath::Power(r/R0,j)*( par[i]*TMath::Cos(j*thr) + par[i+1]*TMath::Sin(j*thr) );  
   }
   return f; 
}
