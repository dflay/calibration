#include "../include/TRLYFuncs.h"
//______________________________________________________________________________
int GetTrolleyProbePosition(int index,double *pos){
   // get the trolley probe position relative to the center probe
   // array pos is indexed as (0,1) = (x,y) 
   double spacing = 17.5;
   double spacing_alt=0; 
   if(index==0){
      return 0.;
   }else if(index==3){
      pos[0] =  0;  
      pos[1] =  0;  
   }else if(index==9){
      pos[0] = -2.*spacing;
      pos[1] =  0;
   }else if(index==5){
      pos[0] =  1.*spacing;  
      pos[1] =  0;  
   }else if(index==15){
      pos[0] =  2.*spacing;
      pos[1] =  0;
   }else if(index==6){
      pos[0] =  0;  
      pos[1] = -2.*spacing;  
   }else if(index==2){
      pos[0] =  0;
      pos[1] =  -1.*spacing;
   }else if(index==4){
      pos[0] =  0;  
      pos[1] =  1.*spacing;  
   }else if(index==12){
      pos[0] =  0;
      pos[1] =  2.*spacing;
   }else if(index==7){
      pos[0] = -1.*spacing;  
      pos[1] = -1.*spacing_alt;  
   }else if(index==8){
      pos[0] =  1.*spacing_alt;
      pos[1] = -1.*spacing;
   }else if(index==10){
      pos[0] = -1.*spacing_alt;  
      pos[1] =  1.*spacing;  
   }else if(index==11){
      pos[0] = -1.*spacing;
      pos[1] =  1.*spacing_alt;
   }else if(index==13){
      pos[0] =  1.*spacing;  
      pos[1] =  1.*spacing_alt;  
   }else if(index==14){
      pos[0] =  2.*spacing_alt;
      pos[1] =  1.*spacing;
   }else if(index==16){
      pos[0] =  2.*spacing_alt;  
      pos[1] = -1.*spacing;  
   }else if(index==17){
      pos[0] =  1.*spacing_alt;
      pos[1] = -2.*spacing_alt;
   }
   return 0;
}
