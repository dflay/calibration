#ifndef IMAGE_PARAMETER_H 
#define IMAGE_PARAMETER_H

// PP image measurement input struct 

typedef struct imageParameter{
   double height;              // height of FXPR 
   int trial;                  // trial number 
   int midasRun;               // MIDAS run of measurement 
   int nmrDAQRun;              // NMR-DAQ run of measurement 
} imageParameter_t;

#endif  
