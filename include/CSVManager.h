#ifndef CSV_MANAGER_H
#define CSV_MANAGER_H

// a generic CSV file reader/writer 

#include <cstdlib> 
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

class CSVManager { 

   private:
      bool fHeaderExists;
      std::vector<std::string> fHeader;
      std::vector< std::vector<std::string> > fData;

      int SplitString(const char delim,const std::string inStr,std::vector<std::string> &out);

      template <typename T>
	 int CheckType(T data){
            int rc=0; 
	    if( std::is_arithmetic<T>::value ){
	       if( std::is_integral<T>::value ){
		  // it's an integer or boolean 
	       }else{
		  // it's a double or float 
	       }
	    }else{
	       rc = 1; // it's a string
            }
	    return rc;
	 } 

   public:
      CSVManager();
      ~CSVManager();

     int Print();
     int ClearData();
     int ReadFile(const char *inpath,bool header=false); 
     int WriteFile(const char *outpath);

     // getter methods
     std::string GetElement_str(int rowIndex,int colIndex);  
     int GetHeader(std::vector<std::string> &header);  
     int GetColumn_byIndex_str(int colIndex,std::vector<std::string> &data); 
     int GetColumn_byName_str(std::string colName,std::vector<std::string> &data);

     template <typename T>
	int GetColumn_byIndex(int colIndex,std::vector<T> &data){
	   // find the data by col index 
	   int NROW = fData.size();
	   T elem;
	   for(int i=0;i<NROW;i++){
	      elem = GetElement<T>(i,colIndex);
	      data.push_back(elem);
	   }
	   return 0;
	}

     template <typename T> 
     int GetColumn_byName(std::string colName,std::vector<T> &data){
	// find the column index by name 
	int NCOL=0,NROW=0,k=-1;
	if(fHeaderExists){
	   NCOL = fHeader.size();
	   for(int i=0;i<NCOL;i++) if(fHeader[i].compare(colName)==0) k = i;
	   if(k>=0){
	      GetColumn_byIndex<T>(k,data);
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

     template <typename T>
	T GetElement(int rowIndex,int colIndex){
	   int NROW = fData.size();
	   int NCOL = fData[0].size();
	   std::string data="NONE";
	   if(NROW>0 && NCOL>0){
	      data = fData[rowIndex][colIndex];
	   }else{
	      std::cout << "[CSVManager::GetElement]: NO data for row " 
                        << rowIndex << ", col " << colIndex << std::endl;
	   }
	   // now parse the string 
	   T x;
	   if( std::is_arithmetic<T>::value ){
	      if( std::is_integral<T>::value ){
		 // it's an integer or boolean 
		 x = std::atoi( data.c_str() );
	      }else{
		 // it's a double or float 
		 x = std::atof( data.c_str() );
	      }
	   }
	   return x;
	}

}; 

#endif 
