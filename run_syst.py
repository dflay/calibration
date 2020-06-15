# !/usr/bin/python 

# A script to run the calibration analysis in systematic mode

import os 
import sys
import json
import csv
import datetime

from myfuncs import utcNow 

isDebug = False

Nargs = len(sys.argv) 
if(Nargs!=4): 
  print("[run_syst]: Invalid input!")
  print("[run_syst]: Usage: python run_syst.py num-iter run-period probe-list") 
  print("[run_syst]: Example: python run_syst.py 10 1 1-17") 
  sys.exit(0) 

N         = int(sys.argv[1]) 
runPeriod = sys.argv[2] 
probeList = sys.argv[3] 
   
print("[run_syst]: ************* STARTING SYSTEMATIC ANALYSIS *************")
print("[run_syst]: Run period  = {0}".format(runPeriod)) 
print("[run_syst]: Probe list  = {0}".format(probeList)) 
print("[run_syst]: Num of iter = {0:01d}".format(N))
print("[run_syst]: ********************************************************")

startTime = utcNow() 
for i in xrange(0,N): 
   cmd = "python run_analysis_batch.py {0} {1}".format(runPeriod,probeList)
   if(isDebug): 
      print(cmd)
   else:  
      os.system(cmd) 
   print("[run_syst]: ********************************************************") 
   print("[run_syst]: ***************** END OF ITERATION {0:03d} *****************".format(i+1)) 
   print("[run_syst]: ********************************************************") 

endTime = utcNow() 

elapsedTime = float(endTime-startTime) 
print("[run_syst]: Elapsed time = {0:.3f} sec".format(elapsedTime))

