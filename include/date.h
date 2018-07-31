#ifndef DATE_H
#define DATE_H

typedef struct date{
   int month;
   int day;
   int year;
   std::string getDateString(){
      char dstr[200];
      sprintf(dstr,"%02d-%02d-%02d",month,day,year-2000); 
      std::string theStr = dstr;
      return theStr;
   } 
} date_t; 

#endif 
