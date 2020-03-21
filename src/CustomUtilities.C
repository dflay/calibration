#include "../include/CustomUtilities.h"
//______________________________________________________________________________
int SplitString(std::string delim,std::string myStr,std::vector<std::string> &out){
   // split a string by a delimiter
   stringstream ss(myStr);
   std::vector<string> result;

   while( ss.good() ){
      std::string substr;
      if(delim.compare(",")==0) std::getline(ss,substr,',');
      if(delim.compare(" ")==0) std::getline(ss,substr,' ');
      out.push_back(substr);
   } 
   return 0;
}
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
int GetFXPRCutTime(std::string inpath,int probe,int index,
                   unsigned long long &tMin,unsigned long long &tMax){

   // cut type: 0 = lower bound time stamp, 1 = upper bound time stamp, 2 = range
   // output is a lower and upper bound of cuts
    
   // first set default values 
   tMin = 0;
   tMax = -1;    

   json jData;
   int rc = gm2fieldUtil::Import::ImportJSON(inpath,jData);
   // create probe key 
   char pr[20];
   sprintf(pr,"probe-%02d",probe);
   std::string probeStr = pr;
  
   int cutType = 0;
 
   auto it_key = jData.find(probeStr);  // this is an iterator 
   if (it_key!=jData.end() ){
      // not at the end of jData -- found the key 
   }else{
      std::cout << Form("[CustomUtilities::GetFXPRCutTime]: No key named: %s.  Returning 0.",probeStr.c_str()) << std::endl;
      cutType = 0;
      return 0;
   }
   
   // get values
   std::string timeStr   = jData[probeStr]["fxpr"]["time"][index];
   std::string state_str = jData[probeStr]["fxpr"]["state"][index];

   // check for empty placeholders (that is, no cut) 
   if( timeStr.compare("NONE")==0 || state_str.compare("NONE")==0 ){
      std::cout << Form("[CustomUtilities::GetFXPRCutTime]: No cut detected!  Returning 0.") << std::endl;
      cutType = 0;
      return 0;
   }

   // check if the time string is actually more than one timestamp (i.e., defines a cut range)
   bool isCutRange=false;
   std::vector<std::string> tsv; 
   std::string::size_type n; 
   n = timeStr.find(',');  // search for a comma character from beginning of string  
   if(n==std::string::npos){
      // no comma, this is a single timestamp
      tsv.push_back(timeStr); 
   }else{
      // comma found, is a cut range 
      isCutRange = true;
      rc = SplitString(",",timeStr,tsv); // split to a vector 
   } 

   // verify cut type  
   if( state_str.compare("lower")==0 ){
      cutType = kLowerBound;
   }else if( state_str.compare("upper")==0){
      cutType = kUpperBound; 
   }else if( state_str.compare("range")==0){
      cutType = kRange; 
   } 

   // consistency check
   if(isCutRange==true && cutType!=2){
      std::cout << "[CustomUtilities::GetFXPRCutTime]: ERROR! Time vector size is > 1, state is not range!" << std::endl;
      std::cout << "                                   Assume lower bound of time = 0 and returning." << std::endl;
      cutType = 0;
      return 0;
   }

   // convert to UTC
   std::vector<unsigned long long> timeCut; 
   int M=0; 
   bool isDST=false;
   unsigned long long TIME=0;
   const int NT = tsv.size();
   std::string substr; 
   for(int i=0;i<NT;i++){
      // determine if STD or DST  
      M = tsv[i].size();
      substr = tsv[i].substr(M-3,M-1);
      if( substr.compare("CST")==0 ){
	 isDST = false; // standard time
      }else{
	 isDST = true;  // daylight savings time 
      }
      TIME = 1E+9*gm2fieldUtil::GetUTCTimeStampFromString(tsv[i],isDST);
      std::cout << Form("[CustomUtilities::GetFXPRCutTime]: time cut = %s (%.0lf)",tsv[i].c_str(),TIME/1E+9) << std::endl;
      timeCut.push_back(TIME); 
   }

    if(cutType==kLowerBound){
      // lower bound
      tMin = timeCut[0];
   }else if(cutType==kUpperBound){
      // upper bound 
      tMax = timeCut[0];
   }else if(cutType==kRange){
      // cut range
      tMin = timeCut[0];
      tMax = timeCut[1];
   }

   return 0;
}
