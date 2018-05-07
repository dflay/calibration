#include "../include/CustomExport.h"
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

