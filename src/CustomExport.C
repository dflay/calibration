#include "../include/CustomExport.h"
//______________________________________________________________________________
int PrintToFile_sccTimes(const char *outpath,std::vector<double> off,std::vector<double> on){
   // print SCC time stamps to file (as a string)  
   const int NOFF = off.size();
   const int NON  = on.size();
   if(NOFF!=NON){
      std::cout << "Unequal number of transitions! " << std::endl;
      std::cout << "off = " << NOFF << " on = " << NON << std::endl;
      return 1;
   }

   unsigned long theTime=0;
   char outStr[200];

   std::ofstream outfile;
   outfile.open(outpath);

   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      for(int i=0;i<NOFF;i++){
         theTime = (unsigned long)off[i];
         sprintf(outStr,"%s",gm2fieldUtil::GetStringTimeStampFromUTC(theTime).c_str());
         outfile << outStr << std::endl;
         theTime = (unsigned long)on[i];
         sprintf(outStr,"%s",gm2fieldUtil::GetStringTimeStampFromUTC(theTime).c_str());
         outfile << outStr << std::endl;
      }
      outfile.close();
   }
   return 0;
}
//______________________________________________________________________________
int PrintTRLYFrequencies(const char *outpath,int probe,std::vector<trolleyAnaEvent_t> data){
   const int N = data.size(); 
   char outStr[500]; 
   std::ofstream outfile;
   outfile.open(outpath);
   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      for(int i=0;i<N;i++){
	 sprintf(outStr,"%.0lf,%.3lf",data[i].time[probe]/1E+9,data[i].freq[probe]);
	 outfile << outStr << std::endl;
      }
      outfile.close();
      std::cout << "The data has been printed to the file: " << outpath << std::endl;
   }
   return 0;
} 
//______________________________________________________________________________
int PrintTRLYPositions(const char *outpath,std::vector<int> probe,
                       std::vector<double> r,std::vector<double> dr,
                       std::vector<double> y,std::vector<double> dy,
                       std::vector<double> z,std::vector<double> dz){
   const int N = probe.size();
   char outStr[500];
   std::ofstream outfile; 
   outfile.open(outpath);
   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      for(int i=0;i<N;i++){
	 sprintf(outStr,"%02d,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf",
                 probe[i],r[i],dr[i],y[i],dy[i],z[i],dz[i]);
         outfile << outStr << std::endl;
      }
      std::cout << "The data has been written to the file: " << outpath << std::endl;
      outfile.close();
   }
   return 0;
}
//______________________________________________________________________________
int PrintToFile(const char *outpath,std::vector<std::string> label,std::vector<double> x){

   char outStr[200]; 
   const int N = label.size(); 

   std::ofstream outfile;
   outfile.open(outpath);
   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      for(int i=0;i<N;i++){
	 sprintf(outStr,"%s,%.3lf",label[i].c_str(),x[i]);
	 outfile << outStr << std::endl;
      }
      outfile.close();
      std::cout << "The data has been written to the file: " << outpath << std::endl;
   }
   return 0;
}
//______________________________________________________________________________
int PrintToFile(const char *outpath,std::vector<std::string> label,
                std::vector<double> x1,std::vector<double> x2,std::vector<double> x3){

   char outStr[200]; 
   const int N = label.size(); 

   std::ofstream outfile;
   outfile.open(outpath);
   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      for(int i=0;i<N;i++){
	 sprintf(outStr,"%s,%.3lf,%.3lf,%.3lf",label[i].c_str(),x1[i],x2[i],x3[i]);
	 outfile << outStr << std::endl;
      }
      outfile.close();
      std::cout << "The data has been written to the file: " << outpath << std::endl;
   }
   return 0;
}
//______________________________________________________________________________
int PrintToFile(const char *outpath,std::vector<std::string> label,
                std::vector<double> x1,std::vector<double> x2,
                std::vector<double> x3,std::vector<double> x4){

   char outStr[200]; 
   const int N = label.size(); 

   std::ofstream outfile;
   outfile.open(outpath);
   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      for(int i=0;i<N;i++){
	 sprintf(outStr,"%s,%.3lf,%.3lf,%.3lf,%.3lf",label[i].c_str(),x1[i],x2[i],x3[i],x4[i]);
	 outfile << outStr << std::endl;
      }
      outfile.close();
      std::cout << "The data has been written to the file: " << outpath << std::endl;
   }
   return 0;
}
//______________________________________________________________________________
int PrintToFile(const char *outpath,std::vector<std::string> label,
                std::vector<double> x1,std::vector<double> x2,
                std::vector<double> x3,std::vector<double> x4,
                std::vector<double> x5,std::vector<double> x6){

   char outStr[200]; 
   const int N = label.size(); 

   std::ofstream outfile;
   outfile.open(outpath);
   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      for(int i=0;i<N;i++){
	 sprintf(outStr,"%s,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf",label[i].c_str(),x1[i],x2[i],x3[i],x4[i],x5[i],x6[i]);
	 outfile << outStr << std::endl;
      }
      outfile.close();
      std::cout << "The data has been written to the file: " << outpath << std::endl;
   }
   return 0;
}
//______________________________________________________________________________
int PrintToFile(const char *outpath,std::string label,double *x,double *x_err,double *y,double *y_err){
   // print results to file
   char outStr[200];
   sprintf(outStr,"%s,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf",
           label.c_str(),
           x[0],x_err[0],x[1],x_err[1],x[2],x_err[2],
           y[0],y_err[0],y[1],y_err[1],y[2],y_err[2]);

   std::ofstream outfile;
   outfile.open(outpath);
   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      outfile << outStr << std::endl;
      outfile.close();
      std::cout << "The data has been written to the file: " << outpath << std::endl;
   }
   return 0;
}
//______________________________________________________________________________
int PrintToFile(const char *outpath,std::string label,double *x,double *x_err){
   // print results to file
   char outStr[200];
   sprintf(outStr,"%s,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf",label.c_str(),x[0],x_err[0],x[1],x_err[1],x[2],x_err[2]);

   std::ofstream outfile;
   outfile.open(outpath);
   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      outfile << outStr << std::endl;
      outfile.close();
      std::cout << "The data has been written to the file: " << outpath << std::endl;
   }
   return 0;
}
//______________________________________________________________________________
int PrintToFile(const char *outpath,double *x){
   // print results to file
   char outStr[200];
   sprintf(outStr,"%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf",x[0],x[1],x[2],x[3],x[4],x[5],x[6]);

   std::ofstream outfile;
   outfile.open(outpath);
   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      outfile << outStr << std::endl;
      outfile.close();
      std::cout << "The data has been written to the file: " << outpath << std::endl;
   }
   return 0;
}
//______________________________________________________________________________
int PrintToFile(const char *outpath,std::string label,const int N,double *x){
   // print results to file
   char outStr[200];
   sprintf(outStr,"%s",label.c_str());
   for(int i=0;i<N;i++){
      sprintf(outStr,"%s,%.3lf",outStr,x[i]);
   }

   std::ofstream outfile;
   outfile.open(outpath,std::ios::app);
   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      outfile << outStr << std::endl;
      outfile.close();
      std::cout << "The data has been written to the file: " << outpath << std::endl;
   }
   return 0;
}
//______________________________________________________________________________
int PrintToFile(const char *outpath,std::string label,const int N,double *x,double *x_err){
   // print results to file
   char outStr[200];
   sprintf(outStr,"%s",label.c_str());
   for(int i=0;i<N;i++){
      sprintf(outStr,"%s,%.3lf,%.3lf",outStr,x[i],x_err[i]);
   }

   std::ofstream outfile;
   outfile.open(outpath,std::ios::app);
   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      outfile << outStr << std::endl;
      outfile.close();
      std::cout << "The data has been written to the file: " << outpath << std::endl;
   }
   return 0;
}
//______________________________________________________________________________
int PrintToFile(const char *outpath,std::vector<double> x){
   char outStr[200];
   const int N = x.size();

   std::ofstream outfile;
   outfile.open(outpath);
   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      for(int i=0;i<N;i++){
         sprintf(outStr,"%.3lf",x[i]);
         outfile << outStr << std::endl;
      }
   }
   return 0;
}
//______________________________________________________________________________
int PrintToFile(const char *outpath,std::vector<double> x1,std::vector<double> x2){
   char outStr[200];
   const int N = x1.size();

   std::ofstream outfile;
   outfile.open(outpath);
   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      for(int i=0;i<N;i++){
         sprintf(outStr,"%.3lf,%.3lf",x1[i],x2[i]);
         outfile << outStr << std::endl;
      }
   }
   return 0;
}
//______________________________________________________________________________
int PrintToFile(const char *outpath,std::vector<std::string> label,
                std::vector<double> x1,std::vector<double> x2){
   char outStr[200];
   const int N = x1.size();

   std::ofstream outfile;
   outfile.open(outpath);
   if( outfile.fail() ){
      std::cout << "Cannot open the file: " << outpath << std::endl;
      return 1;
   }else{
      for(int i=0;i<N;i++){
         sprintf(outStr,"%s,%.3lf,%.3lf",label[i].c_str(),x1[i],x2[i]);
         outfile << outStr << std::endl;
      }
   }
   return 0;
}
