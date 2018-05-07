#ifndef BLIND_H
#define BLIND_H

// a struct for blinding data 

typedef struct blind {
   double value_pp;
   double err_pp;
   double value_tr;
   double err_tr;
   double err_x;
   double err_y;
   double err_z;
} blind_t; 

#endif 
