#include "../include/CustomImport.h"
//______________________________________________________________________________
int LoadCalibSwapData(const char *inpath,std::vector<calibSwap_t> &data){

   int i=0;
   std::string stime,sf,sfe,st,ste;

   calibSwap_t dataPt;

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,stime,',');
         std::getline(infile,sf   ,',');
         std::getline(infile,sfe  ,',');
         std::getline(infile,st   ,',');
         std::getline(infile,ste);
	 dataPt.time    = std::atof( stime.c_str() );
	 dataPt.freq    = std::atof( sf.c_str()    );
	 dataPt.freqErr = std::atof( sfe.c_str()   );
	 data.push_back(dataPt); 
      }
      infile.close();
      data.pop_back();
   }

   return 0;
}
//______________________________________________________________________________
int LoadImposedGradientData(const char *inpath,imposed_gradient_t &data){

   int i=0;
   std::string sp,sgrad,sgrad_err;

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,sp,',');
         std::getline(infile,sgrad,',');
         std::getline(infile,sgrad_err);
	 data.pos      = std::atof( sp.c_str()        );
	 data.grad     = std::atof( sgrad.c_str()     );
	 data.grad_err = std::atof( sgrad_err.c_str() );
      }
      infile.close();
   }

   return 0;
}
//______________________________________________________________________________
int LoadImposedGradientData(const char *inpath,std::vector<imposed_gradient_t> &data){

   int i=0;
   std::string sp,sgrad,sgrad_err;

   imposed_gradient_t dataPt;

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,sp,',');
         std::getline(infile,sgrad,',');
         std::getline(infile,sgrad_err);
	 dataPt.pos      = std::atof( sp.c_str()        );
	 dataPt.grad     = std::atof( sgrad.c_str()     );
	 dataPt.grad_err = std::atof( sgrad_err.c_str() );
	 data.push_back(dataPt); 
      }
      infile.close();
      data.pop_back();
   }

   return 0;
}
//______________________________________________________________________________
int LoadImposedAziGradData(const char *inpath,int probe,double &dBdz){

   int ipr=0;
   std::string stp,sbg,sig_82,sg,sg_per_A,ssc;  // bare grad, grad at 0.82 A, grad, grad/Amp,shim current (A) 

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,stp     ,',');
         std::getline(infile,sbg     ,',');
         std::getline(infile,sig_82  ,',');
         std::getline(infile,sg      ,',');
         std::getline(infile,sg_per_A,',');
         std::getline(infile,ssc);
         ipr = std::atoi( stp.c_str() ); 
         if(ipr==probe) dBdz = std::atof( sg.c_str() ); 
      }
      infile.close();
   }

   return 0;
}
//______________________________________________________________________________
int LoadTrolleyDeltaBData(const char *inpath,trolleyDeltaB_t &data){

   int i=0;
   std::string sp,snq,snqE,ssq,ssqE;

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,sp  ,',');
         std::getline(infile,snq ,',');
         std::getline(infile,snqE,',');
         std::getline(infile,ssq ,',');
         std::getline(infile,ssqE);
         data.probeID[i]      = i+1; 
	 data.normQuad[i]     = std::atof( snq.c_str()  );
	 data.normQuad_err[i] = std::atof( snqE.c_str() );
	 data.skewQuad[i]     = std::atof( ssq.c_str()  );
	 data.skewQuad_err[i] = std::atof( ssqE.c_str() );
	 i++;
      }
      infile.close();
   }

   return 0;
}
//______________________________________________________________________________
int LoadTrolleyPositionData(const char *inpath,trolleyProbePosition_t &data){

   int i=0;
   std::string sp,sr,sdr,sy,sdy,sphi,sdphi;

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,sp,',');
         std::getline(infile,sr,',');
         std::getline(infile,sdr,',');
         std::getline(infile,sy,',');
         std::getline(infile,sdy,',');
         std::getline(infile,sphi,',');
         std::getline(infile,sdphi);
	 data.r[i]    = std::atof( sr.c_str()    );
	 data.dr[i]   = std::atof( sdr.c_str()   );
	 data.y[i]    = std::atof( sy.c_str()    );
	 data.dy[i]   = std::atof( sdy.c_str()   );
	 data.phi[i]  = std::atof( sphi.c_str()  );
	 data.dphi[i] = std::atof( sdphi.c_str() );
	 i++;
      }
      infile.close();
   }

   return 0;
}
//______________________________________________________________________________
int ImportNMRANAData(const char *inpath,std::vector<nmrAnaEvent_t> &Data){
   // load data from the NMR-ANA framework 

   nmrAnaEvent_t inData;

   int N=0,k=0;
   int irun,ipulse,izc,NUM_PULSES=0;
   double inoise,inc,ifa,ifb,ifc,ifa_ph,ifb_ph,ifc_ph,ivmax;

   // const int MAX = 1000;
   // char buf[MAX];

   ifstream infile;
   infile.open(inpath,ios::in);

   if(infile.fail()){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      // cout << "Opening the file: " << inpath << endl;
      // for(int i=0;i<1;i++) infile.getline(buf,MAX);
      while( !infile.eof() ){
         infile >> irun >> ipulse >> ivmax >> inoise >> izc >> inc >> ifa >> ifb >> ifc >> ifa_ph >> ifb_ph >> ifc_ph;
         inData.run     = irun;
         inData.ampl    = ivmax;
         inData.noise   = inoise;
         inData.pulse   = ipulse; 
         inData.zc      = izc; 
         inData.nc      = inc; 
         inData.freq[0] = ifa;
         inData.freq[1] = ifb;
         inData.freq[2] = ifc;
         inData.freq[3] = ifa_ph;
         inData.freq[4] = ifb_ph;
         inData.freq[5] = ifc_ph;
         Data.push_back(inData);
      }
      infile.close();
      Data.pop_back();
   }

   return 0;
}
//______________________________________________________________________________
int ImportPPEncoderCalib_csv(const char *inpath,std::vector<double> &y){

   int N=0,ix1=0,ix2=0,ix3=0;
   string sx1,sx2,sx3;

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,sx1,',');
         std::getline(infile,sx2,',');
         std::getline(infile,sx3);
         ix1 = std::atof( sx1.c_str() );
         ix2 = std::atof( sx2.c_str() );
         ix3 = std::atof( sx3.c_str() );
         y.push_back(ix1);
         y.push_back(ix2);
         y.push_back(ix3);
      }
      infile.close();
      N = y.size();
      if(N>3) y.pop_back();
   }

   return 0;
}
//______________________________________________________________________________
int ImportDeltaBFileList_csv(const char *inpath,
                    std::vector<int> &x1,std::vector<std::string> &x2,
                    std::vector<double> &x3){

   int ix1=0;
   double ix3;
   string sx1,sx2,sx3;

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,sx1,',');
         std::getline(infile,sx2,',');
         std::getline(infile,sx3);
         ix1 = atoi( sx1.c_str() );
         ix3 = atof( sx3.c_str() );
         x1.push_back(ix1);
         x2.push_back(sx2);
         x3.push_back(ix3);
      }
      infile.close();
      x1.pop_back();
      x2.pop_back();
      x3.pop_back();
   }

   return 0;
}
//______________________________________________________________________________
int LoadPerturbationData(const char *inpath,perturbation_t &pert){

   int cntr=0;
   string s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12;

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,s1,',');
         std::getline(infile,s2,',');
         std::getline(infile,s3,',');
         std::getline(infile,s4,',');
         std::getline(infile,s5,',');
         std::getline(infile,s6,',');
         std::getline(infile,s7,',');
         std::getline(infile,s8,',');
         std::getline(infile,s9,',');
         std::getline(infile,s10,',');
         std::getline(infile,s11,',');
         std::getline(infile,s12);
	 cntr++; 
	 if(cntr==1){
	    pert.sigma         = std::atof( s1.c_str() );
	    pert.sigma_err     = std::atof( s2.c_str() );
	    pert.chi           = std::atof( s3.c_str() );
	    pert.chi_err       = std::atof( s4.c_str() );
	    pert.eps           = std::atof( s5.c_str() );
	    pert.eps_err       = std::atof( s6.c_str() );
	    pert.delta_m       = std::atof( s7.c_str() );
	    pert.delta_m_err   = std::atof( s8.c_str() );
	    pert.delta_eps     = std::atof( s9.c_str() );
	    pert.delta_eps_err = std::atof( s10.c_str() );
	    pert.delta_mag     = std::atof( s11.c_str() );
	    pert.delta_mag_err = std::atof( s12.c_str() );
	 }
      }
      infile.close();
   }
  
   return 0;
}

//______________________________________________________________________________
int LoadFieldData(const char *inpath,nmr_meas_t &data){

   int cntr=0;
   std::string sL,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12;

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,sL,',');
         std::getline(infile,s1,',');
         std::getline(infile,s2,',');
         std::getline(infile,s3,',');
         std::getline(infile,s4,',');
         std::getline(infile,s5,',');
         std::getline(infile,s6,',');
         std::getline(infile,s7,',');
	 std::getline(infile,s8,',');
	 std::getline(infile,s9,',');
	 std::getline(infile,s10,',');
	 std::getline(infile,s11,',');
	 std::getline(infile,s12);
	 cntr++;
	 if(cntr==1){
	    data.name          = sL;
	    data.freq          = std::atof( s1.c_str() );
	    data.freq_err      = std::atof( s2.c_str() );
	    data.freq_fxpr     = std::atof( s3.c_str() );
	    data.freq_fxpr_err = std::atof( s4.c_str() );
	    data.freq_trly     = std::atof( s5.c_str() );
	    data.freq_trly_err = std::atof( s6.c_str() );
            data.p2p_err       = std::atof( s8.c_str() ); 
            data.r2r_err       = 0.;  // no run-to-run uncertainty here  
	 }
      }
      infile.close();
   }

   std::cout << data.name << " " 
             << Form("raw     = %.3lf +/- %.3lf",data.freq,data.freq_err)           << " " 
             << Form("fxpr    = %.3lf +/- %.3lf",data.freq_fxpr,data.freq_fxpr_err) << " "  
             << Form("trly    = %.3lf +/- %.3lf",data.freq_trly,data.freq_trly_err) << " " 
             << Form("p2p_err = %.3lf",data.p2p_err)                                << std::endl; 

   return 0;
}
//______________________________________________________________________________
int LoadGradientData(const char *inpath,grad_meas_t &x){

   grad_meas_t data;

   std::string sL,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12;

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,sL,',');
         std::getline(infile,s1,',');
         std::getline(infile,s2,',');
         std::getline(infile,s3,',');
         std::getline(infile,s4,',');
         std::getline(infile,s5,',');
         std::getline(infile,s6,',');
         std::getline(infile,s7,',');
         std::getline(infile,s8,',');
         std::getline(infile,s9,',');
         std::getline(infile,s10,',');
         std::getline(infile,s11,',');
         std::getline(infile,s12);
	 x.name           = sL;
         x.grad           = std::atof( s1.c_str() );
         x.grad_err       = std::atof( s2.c_str() );
         x.grad_fxpr      = std::atof( s3.c_str() );
         x.grad_fxpr_err  = std::atof( s4.c_str() );
         x.grad_trly      = std::atof( s5.c_str() );
	 x.grad_trly_err  = std::atof( s6.c_str() );
         x.drift_fxpr     = std::atof( s9.c_str() );
         x.drift_fxpr_err = std::atof( s10.c_str() );
         x.drift_trly     = std::atof( s11.c_str() );
	 x.drift_trly_err = std::atof( s12.c_str() );
      }
      infile.close();
   }

   // const int N = x.size();
   // for(int i=0;i<N;i++){
   //    std::cout << x[i].name << " " 
   //              << x[i].grad      << " +/- " << x[i].grad_err      << " " 
   //              << x[i].grad_fxpr << " +/- " << x[i].grad_fxpr_err << " " 
   //              << x[i].grad_trly << " +/- " << x[i].grad_trly_err << std::endl;
   // }
   // std::cout << "--------------" << std::endl; 

   return 0;
}
//______________________________________________________________________________
int LoadGradientData(const char *inpath,std::vector<grad_meas_t> &x){

   grad_meas_t data;

   std::string sL,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12;

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,sL,',');
         std::getline(infile,s1,',');
         std::getline(infile,s2,',');
         std::getline(infile,s3,',');
         std::getline(infile,s4,',');
         std::getline(infile,s5,',');
         std::getline(infile,s6,',');
         std::getline(infile,s7,',');
         std::getline(infile,s8,',');
         std::getline(infile,s9,',');
         std::getline(infile,s10,',');
         std::getline(infile,s11,',');
         std::getline(infile,s12);
	 data.name           = sL;
         data.grad           = std::atof( s1.c_str() );
         data.grad_err       = std::atof( s2.c_str() );
         data.grad_fxpr      = std::atof( s3.c_str() );
         data.grad_fxpr_err  = std::atof( s4.c_str() );
         data.grad_trly      = std::atof( s5.c_str() );
	 data.grad_trly_err  = std::atof( s6.c_str() );
         data.drift_fxpr     = std::atof( s9.c_str() );
         data.drift_fxpr_err = std::atof( s10.c_str() );
         data.drift_trly     = std::atof( s11.c_str() );
	 data.drift_trly_err = std::atof( s12.c_str() );
	 x.push_back(data);
      }
      infile.close();
      x.pop_back();
   }

   // const int N = x.size();
   // for(int i=0;i<N;i++){
   //    std::cout << x[i].name << " " 
   //              << x[i].grad      << " +/- " << x[i].grad_err      << " " 
   //              << x[i].grad_fxpr << " +/- " << x[i].grad_fxpr_err << " " 
   //              << x[i].grad_trly << " +/- " << x[i].grad_trly_err << std::endl;
   // }
   // std::cout << "--------------" << std::endl; 

   return 0;
}
//______________________________________________________________________________
int LoadDeltaBData_trlyXY(const char *inpath,int probe,std::vector<deltab_t> &x){

   int ipr=0;
   deltab_t data;
   std::string stp,s1,s2,s3,s4,s5,s6,s7,s8,s9;

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,stp,',');
         std::getline(infile,s1 ,',');
         std::getline(infile,s2 ,',');
         std::getline(infile,s3 ,',');
         std::getline(infile,s4 ,',');
         std::getline(infile,s5 ,',');
         std::getline(infile,s6 ,',');
         std::getline(infile,s7 ,',');
         std::getline(infile,s8);
         ipr = std::atoi( stp.c_str() );
	 // std::cout << "PROBE " << ipr << std::endl;
         if(ipr==probe){
	    // std::cout << "--> MATCH" << std::endl;
	    // radial 
	    data.name = "rad-grad";
	    data.dB             = std::atof( s1.c_str()  );
	    data.dB_err         = std::atof( s2.c_str()  );
	    data.dB_fxpr        = std::atof( s3.c_str()  );
	    data.dB_fxpr_err    = std::atof( s4.c_str()  );
	    data.dB_trly        = 0;
	    data.dB_trly_err    = 0;
	    x.push_back(data);
	    // vertical 
	    data.name = "vert-grad";
	    data.dB             = std::atof( s5.c_str()  );
	    data.dB_err         = std::atof( s6.c_str()  );
	    data.dB_fxpr        = std::atof( s7.c_str()  );
	    data.dB_fxpr_err    = std::atof( s8.c_str()  );
	    data.dB_trly        = 0;
	    data.dB_trly_err    = 0;
	    x.push_back(data);
	 }
      }
      infile.close();
      // x.pop_back();
   }

   // const int N = x.size();
   // std::cout << N << " entries found" << std::endl;
   // for(int i=0;i<N;i++){
   //    std::cout << x[i].name << " " 
   //              << x[i].dB         << " +/- " << x[i].dB_err         << " " 
   //              << x[i].dB_fxpr    << " +/- " << x[i].dB_fxpr_err    << " " 
   //              << x[i].dB_trly    << " +/- " << x[i].dB_trly_err    << " " 
   //              << x[i].drift_fxpr << " +/- " << x[i].drift_fxpr_err << " " 
   //              << x[i].drift_trly << " +/- " << x[i].drift_trly_err << std::endl;
   // }
   // std::cout << "--------------" << std::endl; 

   return 0;
}
//______________________________________________________________________________
int LoadDeltaBData(const char *inpath,std::vector<deltab_t> &x){

   deltab_t data;
   std::string sL,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12;

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,sL,',');
         std::getline(infile,s1,',');
         std::getline(infile,s2,',');
         std::getline(infile,s3,',');
         std::getline(infile,s4,',');
         std::getline(infile,s5,',');
         std::getline(infile,s6,','); 
         std::getline(infile,s7,',');
         std::getline(infile,s8,',');
         std::getline(infile,s9,',');
         std::getline(infile,s10,',');
         std::getline(infile,s11,',');
         std::getline(infile,s12);
	 data.name           = sL;
         data.dB             = std::atof( s1.c_str()  );
         data.dB_err         = std::atof( s2.c_str()  );
         data.dB_fxpr        = std::atof( s3.c_str()  );
         data.dB_fxpr_err    = std::atof( s4.c_str()  );
         data.dB_trly        = std::atof( s5.c_str()  );
	 data.dB_trly_err    = std::atof( s6.c_str()  );
         data.drift_fxpr     = std::atof( s9.c_str()  );   // WARNING: note that we're starting at the 9th entry! 
         data.drift_fxpr_err = std::atof( s10.c_str() );
         data.drift_trly     = std::atof( s11.c_str() );
	 data.drift_trly_err = std::atof( s12.c_str() );
	 x.push_back(data);
      }
      infile.close();
      x.pop_back();
   }

   // const int N = x.size();
   // for(int i=0;i<N;i++){
   //    std::cout << x[i].name << " " 
   //              << x[i].dB         << " +/- " << x[i].dB_err         << " " 
   //              << x[i].dB_fxpr    << " +/- " << x[i].dB_fxpr_err    << " " 
   //              << x[i].dB_trly    << " +/- " << x[i].dB_trly_err    << " " 
   //              << x[i].drift_fxpr << " +/- " << x[i].drift_fxpr_err << " " 
   //              << x[i].drift_trly << " +/- " << x[i].drift_trly_err << std::endl;
   // }
   // std::cout << "--------------" << std::endl; 

   return 0;
}
//______________________________________________________________________________
int FindStartIndexTRLY(std::string date,int runNumber){
   // find the true start index to keep events for the trolley data because of the 
   // issue with interactive vs continuous mode
   int index=22,theRun=0; 
   std::string sr,sL;

   char inpath[500]; 
   sprintf(inpath,"./input/runlists/%s/trly-mode.csv",date.c_str() );

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,sr,',');
         std::getline(infile,sL);
         theRun = std::atoi(sr.c_str()); 
         if(runNumber==theRun){
	    std::cout << "Found run " << runNumber << std::endl;
	    if( sL.compare("interactive")==0 ) index = 0;
	    if( sL.compare("continuous")==0  ) index = 22;
         }
      }
      infile.close();
   }

   std::cout << "--> For run " << runNumber << ", the start index is " << index << std::endl; 

   return index;
  
}
//______________________________________________________________________________
int ImportBlinding(blind_t &data){
   // load the blinding data 
   int cntr=0;
   std::string sv,sv2,serr_x,serr_y,serr_z,serr_f1,serr_f2;

   std::string inpath = "./misc/blind/values.csv"; 

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,sv,',');
         std::getline(infile,sv2,',');
         std::getline(infile,serr_x,',');
         std::getline(infile,serr_y,',');
         std::getline(infile,serr_z,',');
         std::getline(infile,serr_f1,',');
         std::getline(infile,serr_f2);
         cntr++;
	 if(cntr==1){ // for some reason we need this...
	    data.value_pp = std::atof( sv.c_str() ); 
	    data.value_tr = std::atof( sv2.c_str() ); 
	    data.err_x    = std::atof( serr_x.c_str() ); 
	    data.err_y    = std::atof( serr_y.c_str() ); 
	    data.err_z    = std::atof( serr_z.c_str() ); 
	    data.err_pp   = std::atof( serr_f1.c_str() ); 
	    data.err_tr   = std::atof( serr_f2.c_str() );
	 } 
      }
      infile.close();
   }

   return 0;
}
//______________________________________________________________________________
int SortRuns(std::vector<std::string> label,std::vector<int> allRuns,
             std::vector<int> &run,std::vector<int> &driftRun,std::vector<int> &index){

   const int NRUN = allRuns.size(); 
   int bareRun=0,bareIndex=0;
   int gradRun=0,gradIndex=0;
   int bare2Run=0,bare2Index=0;
   for(int i=0;i<NRUN;i++){
      if(label[i].compare("bare")==0){
         bareRun   = allRuns[i];
         bareIndex = i;
         run.push_back(allRuns[i]);
      }else if(label[i].compare("bare-2")==0){
         bare2Run   = allRuns[i];
         bare2Index = i;
      }else{
         gradRun   = allRuns[i];
         gradIndex = i;
         run.push_back(allRuns[i]);
      }
   }

   // drift runs 
   driftRun.push_back( allRuns[bareIndex]  ); 
   driftRun.push_back( allRuns[bare2Index] );

   // save indices 
   index.push_back( bareIndex  );  
   index.push_back( gradIndex  );  
   index.push_back( bare2Index );  

   return 0;

}
//______________________________________________________________________________
int SortRunsAlt(std::vector<std::string> label,std::vector<int> allRuns,
             std::vector<int> &run,std::vector<int> &index){

   const int NRUN = allRuns.size(); 
   int bareRun=0,bareIndex=0;
   int gradRun=0,gradIndex=0;
   int bare2Run=0,bare2Index=0;
   for(int i=0;i<NRUN;i++){
      if(label[i].compare("bare")==0){
         bareIndex = i;
      }else if(label[i].compare("bare-2")==0){
         bare2Index = i;
      }else{
         gradIndex = i;
      }
      run.push_back(allRuns[i]);
   }

   // save indices 
   index.push_back( bareIndex  );  
   index.push_back( gradIndex  );  
   index.push_back( bare2Index );  

   return 0;

}
//______________________________________________________________________________
int ImportResults(std::string inpath,result_t &data){
   // load the results data 
   int cntr=0;
   std::string sl,s1,s1err,s2,s2err,s3,s3err;

   ifstream infile;
   infile.open(inpath);
   if( infile.fail() ){
      cout << "Cannot open the file: " << inpath << endl;
      return 1;
   }else{
      while( !infile.eof() ){
         std::getline(infile,sl,',');
         std::getline(infile,s1,',');
         std::getline(infile,s1err,',');
         std::getline(infile,s2,',');
         std::getline(infile,s2err,',');
         std::getline(infile,s3,',');
         std::getline(infile,s3err);
         if( sl.compare("pp-raw")==0){
	    data.ppRaw          = std::atof( s1.c_str() ); 
	    data.ppRaw_err      = std::atof( s1err.c_str() );
	    data.ppRaw_fxpr     = std::atof( s2.c_str() ); 
	    data.ppRaw_fxpr_err = std::atof( s2err.c_str() ); 
	    data.ppRaw_trly     = std::atof( s3.c_str() ); 
	    data.ppRaw_trly_err = std::atof( s3err.c_str() ); 
	 }else if(sl.compare("pp-free")==0){
	    data.ppFree          = std::atof( s1.c_str() ); 
	    data.ppFree_err      = std::atof( s1err.c_str() );
	    data.ppFree_fxpr     = std::atof( s2.c_str() ); 
	    data.ppFree_fxpr_err = std::atof( s2err.c_str() ); 
	    data.ppFree_trly     = std::atof( s3.c_str() ); 
	    data.ppFree_trly_err = std::atof( s3err.c_str() ); 
         }else if(sl.compare("trly")==0){
	    data.trly          = std::atof( s1.c_str() ); 
	    data.trly_err      = std::atof( s1err.c_str() );
	    data.trly_fxpr     = std::atof( s2.c_str() ); 
	    data.trly_fxpr_err = std::atof( s2err.c_str() ); 
	    data.trly_trly     = std::atof( s3.c_str() ); 
	    data.trly_trly_err = std::atof( s3err.c_str() ); 
         }else if(sl.compare("drift-shim")==0){
	    data.driftShim          = std::atof( s1.c_str() ); 
	    data.driftShim_err      = std::atof( s1err.c_str() );
	    data.driftShim_fxpr     = std::atof( s2.c_str() ); 
	    data.driftShim_fxpr_err = std::atof( s2err.c_str() ); 
	    data.driftShim_trly     = std::atof( s3.c_str() ); 
	    data.driftShim_trly_err = std::atof( s3err.c_str() ); 
         }else if(sl.compare("diff")==0){
	    data.diff          = std::atof( s1.c_str() ); 
	    data.diff_err      = std::atof( s1err.c_str() );
	    data.diff_fxpr     = std::atof( s2.c_str() ); 
	    data.diff_fxpr_err = std::atof( s2err.c_str() ); 
	    data.diff_trly     = std::atof( s3.c_str() ); 
	    data.diff_trly_err = std::atof( s3err.c_str() ); 
         }
      }
      infile.close();
   }

   return 0;
}
