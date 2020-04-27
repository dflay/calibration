#ifndef MISALIGNMENT_H 
#define MISALIGNMENT_H

#define NUM_AXES_MISCOR 3

// a data structure to hold misalignment data 
typedef struct misalignment{
   double dx;             // in mm (raw) 
   double dy;             // in mm (raw)   
   double dz;             // in mm (raw)    

   double dx_err;         // in mm (raw) 
   double dy_err;         // in mm (raw)   
   double dz_err;         // in mm (raw)    

   double dx_aba;         // in mm (ABA)   
   double dy_aba;         // in mm (ABA)  
   double dz_aba;         // in mm (ABA)   

   double dx_aba_err;     // in mm (ABA)   
   double dy_aba_err;     // in mm (ABA)  
   double dz_aba_err;     // in mm (ABA)   
 
   double dx_opt;         // in mm (opt)     
   double dy_opt;         // in mm (opt)    
   double dz_opt;         // in mm (opt)    

   double dx_opt_err;     // in mm (opt)     
   double dy_opt_err;     // in mm (opt)    
   double dz_opt_err;     // in mm (opt)    

   double dz_bar_opt;
   double dz_bar_opt_err;
 
   // double dz_bar_opt;     // in mm (opt) uses barcode 
   // double dx_opt_err;     // in mm (opt)     
   // double dy_opt_err;     // in mm (opt)    
   // double dz_opt_err;     // in mm (opt)     
   // double dz_bar_opt_err; // in mm (opt) uses barcode 
   // double dB_x;           // in Hz (raw)     
   // double dB_y;           // in Hz (raw)  
   // double dB_z;           // in Hz (raw)   
   // double dB_x_aba;       // in Hz (ABA)       
   // double dB_y_aba;       // in Hz (ABA)    
   // double dB_z_aba;       // in Hz (ABA)   
   // double dB_x_opt;       // in Hz (opt)     
   // double dB_y_opt;       // in Hz (opt)      
   // double dB_z_opt;       // in Hz (opt)  
   // double dB_z_bar_opt; // in Hz (opt) uses barcode 
   // double dB_x_opt_err;  // in Hz (opt)     
   // double dB_y_opt_err;  // in Hz (opt)      
   // double dB_z_opt_err;  // in Hz (opt)  
   // double dB_z_bar_opt_Err; // in Hz (opt) uses barcode 
} misalignment_t;

// for the misalignment *correction* in Hz  
typedef struct misalignCor {
   std::string label;  
   double val;
   double err;
   // barcode value
   double val_bar;
   double err_bar;
} misalignCor_t;  

#endif 
