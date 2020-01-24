#ifndef ABA_H
#define ABA_H

// a class that calculates the ABA difference between a series of measurements 
#include <cstdlib>
#include <iostream>
#include <vector>

class ABA {

   private:
      int fVerbosity;
      bool fUseTimeWeight;

      void Init();

   public:
      ABA();
      ~ABA();

      void UseTimeWeight(bool v=true) { fUseTimeWeight = v; }
      void SetVerbosity(int v)        { fVerbosity     = v; }

      int GetVerbosity()          const { return fVerbosity;     }
      bool GetTimeWeightStatus()  const { return fUseTimeWeight; }

      int GetDifference(std::vector<double> A_time   ,std::vector<double> A,std::vector<double> A_err,
                        std::vector<double> B_time   ,std::vector<double> B,std::vector<double> B_err,
                        std::vector<double> &diff_aba,std::vector<double> &diff_aba_err);

};

#endif
