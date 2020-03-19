#include "../include/CustomUtilities.h"
//______________________________________________________________________________
int PrintMessage(std::string label,std::string msg,int lineNumber){
   if(lineNumber>=0){
      std::cout << "[" << label << "]: line " << lineNumber << ": " << msg << std::endl;
   }else{
      std::cout << "[" << label << "]: " << msg << std::endl;
   }
   return 0; 
}
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
std::string GetPath(std::string base,bool isBlind,std::string blindLabel,std::string date,
                    bool isSyst,int systDirNum){
   // create the path based upon the top-level directory, blinding, and the date provided
   int rc=0;
   char theDir[200];
   if(isBlind){
      sprintf(theDir,"./%s/blinded/%s/%s",base.c_str(),blindLabel.c_str(),date.c_str());
   }else{
      sprintf(theDir,"./%s/unblinded/%s",base.c_str(),date.c_str());
   }

   if(isSyst){
      sprintf(theDir,"%s/syst-%03d",theDir,systDirNum); 
   }

   std::string THE_DIR = theDir; 
   // std::cout << "[GetPath]: Using directory '" << theDir << "'" << std::endl; 

   return THE_DIR;
}
//______________________________________________________________________________
TLine **GetLines(int color,double min,double max,std::vector<double> x){
   const int N = x.size();
   TLine **L = new TLine*[N];
   for(int i=0;i<N;i++){
      L[i] = new TLine(x[i],min,x[i],max);
      L[i]->SetLineWidth(2);
      L[i]->SetLineStyle(2);
      L[i]->SetLineColor(color);
   }
   return L;
}
//______________________________________________________________________________
unsigned long long GetFXPRCutTime(std::string inpath,int probe,int index,bool &isMax){

   json jData;
   int rc = gm2fieldUtil::Import::ImportJSON(inpath,jData);
   // create probe key 
   char pr[20];
   sprintf(pr,"probe-%02d",probe);
   std::string probeStr = pr;
   
   auto it_key = jData.find(probeStr);  // this is an iterator 
   if (it_key!=jData.end() ){
      // not at the end of jData -- found the key 
   }else{
      std::cout << Form("[CustomUtilities::GetFXPRCutTime]: No key named: %s.  Returning 0.",probeStr.c_str()) << std::endl;
      isMax = false;
      return 0;
   }
   
   // get values
   std::string timeStr   = jData[probeStr]["fxpr"]["time"][index];
   std::string state_str = jData[probeStr]["fxpr"]["state"][index];

   // check for empty placeholders (that is, no cut) 
   if( timeStr.compare("NONE")==0 || state_str.compare("NONE")==0 ){
      isMax = false;
      return 0;
   }

   // is this time an upper bound? 
   isMax = false;
   if( state_str.compare("before")==0 ) isMax = true;  

  

   // determine if STD or DST  
   bool isDST = false;
   int M = timeStr.size();
   std::string substr = timeStr.substr(M-3,M-1);
   if( substr.compare("CST")==0 ){
      isDST = false; // standard time
   }else{
      isDST = true;  // daylight savings time 
   }
   unsigned long long TIME = 1E+9*gm2fieldUtil::GetUTCTimeStampFromString(timeStr,isDST);

   return TIME;
}
