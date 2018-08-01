// test the input manager 

#include <cstdlib> 
#include <iostream>
#include <vector>
#include <string> 

#include "./src/InputManager.C"

int Test(){

   std::string inpath = "test.json";
   InputManager *inputMgr = new InputManager(); 
   inputMgr->Load(inpath);
   // inputMgr->Print();

   std::string testKey = "my-key"; 
   bool isThere = inputMgr->DoesKeyExist("my-key"); 
   std::cout << Form("The key '%s' was found: ",testKey.c_str()) << isThere << std::endl;
   if(isThere){
      std::cout << Form("key value: %s",inputMgr->GetValue(testKey).c_str()) << std::endl;
   }

   delete inputMgr; 

   return 0;
} 
