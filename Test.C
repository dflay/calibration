// test the input manager 

#include <cstdlib> 
#include <iostream>
#include <vector>
#include <string> 

#include "./src/InputManager.C"

int Test(){

   std::string inpath = "./test.json";
   InputManager *inputMgr = new InputManager(); 
   inputMgr->Load(inpath);

   inputMgr->Print();

   delete inputMgr; 

   return 0;
} 
