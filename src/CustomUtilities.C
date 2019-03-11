#include "../include/CustomUtilities.h"
//______________________________________________________________________________
int GetDate(date_t &aDate){

   std::time_t t = std::time(0);   // get time now
   std::tm* now = std::localtime(&t);

   aDate.month = now->tm_mon + 1; 
   aDate.day   = now->tm_mday;
   aDate.year  = now->tm_year + 1900;

   return 0;
}
//______________________________________________________________________________
int MakeDirectory(const char *path){
   int rc=0;
   struct stat SB;

   const int SIZE = 200;
   char *dir_path = new char[SIZE];

   sprintf(dir_path,"%s",path);

   // check to see if the directory exists
   if( (stat(dir_path,&SB)==0)&&(S_ISDIR(SB.st_mode)) ){
      // the directory exists!  Do nothing.  
      // You can also check for other file types using other S_IS* macros:
      // S_ISDIR — directory
      // S_ISREG — regular file
      // S_ISCHR — character device
      // S_ISBLK — block device
      // S_ISFIFO — FIFO
      // S_ISLNK — symbolic link
      // S_ISSOCK — socket
   }else{
      rc = mkdir(dir_path,0700);
   }

   if(rc!=0) std::cout << "[MakeDirectory]: Cannot make the directory: " << path << std::endl;

   delete[] dir_path;

   return rc;
}
//______________________________________________________________________________
std::string GetPath(std::string base,bool isBlind,std::string blindLabel,std::string date){
   // create the path based upon the top-level directory, blinding, and the date provided
   int rc=0;
   char theDir[200];
   if(isBlind){
      sprintf(theDir,"./%s/blinded/%s/%s",base.c_str(),blindLabel.c_str(),date.c_str());
   }else{
      sprintf(theDir,"./%s/unblinded/%s",base.c_str(),date.c_str());
   }

   std::string THE_DIR = theDir; 
   // std::cout << "[GetPath]: Using directory '" << theDir << "'" << std::endl; 

   return THE_DIR;
}
