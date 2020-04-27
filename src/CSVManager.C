#include "../include/CSVManager.h"
//______________________________________________________________________________
CSVManager::CSVManager(){
   fHeaderExists = false;
}
//______________________________________________________________________________
CSVManager::~CSVManager(){

}
//______________________________________________________________________________
int CSVManager::ClearData(){
   fHeader.clear();
   int NROW = fData.size();
   int NCOL = fData[0].size();
   for(int i=0;i<NROW;i++){
      fData[i].clear();
   }
   fData.clear();
   return 0;
}
//______________________________________________________________________________
int CSVManager::SplitString(const char delim,const std::string inStr,std::vector<std::string> &out){
   // split a string by a delimiter
   std::stringstream ss(inStr);
   while( ss.good() ){
      std::string substr;
      std::getline(ss,substr,delim);
      out.push_back(substr);
   }
   return 0;
}
//______________________________________________________________________________
int CSVManager::ReadFile(const char *inpath,bool headerExists){

   // update variable 
   fHeaderExists = headerExists;
 
   std::string aLine;
   std::vector<std::string> line;

   std::ifstream infile;
   infile.open(inpath);

   if( infile.fail() ){
      std::cout << "[CSVManager::ReadFile]: Cannot open the file: " << inpath << std::endl;
      return 1;
   }else{
      while( !infile.eof() ){
	 std::getline(infile,aLine);
         line.push_back(aLine); 
      }
      line.pop_back();
      infile.close();
   }
   
   int NROW = line.size();
   std::vector<std::string> row;
   // fData.resize(NROW); 

   // now parse the data
   int rc=0,k=0,NCOL=0; 
   for(int i=0;i<NROW;i++){
      // split the line into a vector.  This is a single row 
      rc   = SplitString(',',line[i],row);
      // get number of columns 
      NCOL = row.size();
      // set up the size of the columns
      // fData[i].resize(NCOL);
      if(i==0 && fHeaderExists){
	 // this is the header
	 for(int j=0;j<NCOL;j++) fHeader.push_back(row[j]);  
      }else{
	 // not the header, fill the data 
	 fData.push_back(row);
      } 
      // clean up
      row.clear(); 
   }

   NROW = fData.size();

   char msg[200]; 
   if(fHeaderExists)  sprintf(msg,"[CSVManager::ReadFile]: Found header, %d rows, %d columns",NROW,NCOL);
   if(!fHeaderExists) sprintf(msg,"[CSVManager::ReadFile]: Found NO header, %d rows, %d columns",NROW,NCOL);
   std::cout << msg << std::endl;

   return 0; 
}
//______________________________________________________________________________
std::string CSVManager::GetElement_str(int rowIndex,int colIndex){
   int NROW = fData.size();
   int NCOL = fData[0].size();
   std::string data="NONE";
   if(NROW>0 && NCOL>0){
      data = fData[rowIndex][colIndex];
   }else{
      std::cout << "[CSVManager::GetElement]: NO data for row "
	 << rowIndex << ", col " << colIndex << std::endl;
   }
   return data;
}
//______________________________________________________________________________
int CSVManager::GetColumn_byIndex_str(int colIndex,std::vector<std::string> &data){
   // find the data by col index 
   int NROW = fData.size();
   std::string elem;
   for(int i=0;i<NROW;i++){
      elem = GetElement_str(i,colIndex);
      data.push_back(elem);
   }
   return 0;
}
//______________________________________________________________________________
int CSVManager::GetColumn_byName_str(std::string colName,std::vector<std::string> &data){
   // find the column index by name 
   int NCOL=0,NROW=0,k=-1;
   if(fHeaderExists){
      NCOL = fHeader.size();
      for(int i=0;i<NCOL;i++) if(fHeader[i].compare(colName)==0) k = i;
      if(k>=0){
	 GetColumn_byIndex_str(k,data);
      }else{
	 std::cout << "[CSVManager::GetColumn_byName]: Cannot find the key '" 
                   << colName << "' in header!" << std::endl;
	 return 1;
      }
   }else{
      std::cout << "[CSVManager::GetColumn_byName]: No header to search!";
      std::cout << "  Try CSVManager::GetColumn_byName" << std::endl;
      return 1;
   }
   return 0;
}
//______________________________________________________________________________
int CSVManager::GetHeader(std::vector<std::string> &header){
   int N=0;
   if(fHeaderExists){
      N = fHeader.size();
      for(int i=0;i<N;i++) header.push_back(fHeader[i]);
   }else{
      std::cout << "[CSVManager::GetHeader]: No header!" <<  std::endl;
      return 1;
   }
   return 0;
}
//______________________________________________________________________________
int CSVManager::Print(){

   int NROW = fData.size();
   int NCOL = fData[0].size();

   char myStr[200];
   if(fHeaderExists){ 
      sprintf(myStr,"%s",fHeader[0].c_str() );
      for(int i=1;i<NCOL;i++) sprintf(myStr,"%s,%s",myStr,fHeader[i].c_str() );
      std::cout << myStr << std::endl;
   }
  
   for(int i=0;i<NROW;i++){
      sprintf(myStr,"%s",fData[i][0].c_str()); 
      for(int j=1;j<NCOL;j++) sprintf(myStr,"%s,%s",myStr,fData[i][j].c_str() );
      std::cout << myStr << std::endl;
   }
 
   return 0;
}
