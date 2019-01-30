# !/usr/bin/python 

# A script to run systematic studies 

import os 
import sys
import json
import csv
import datetime

from ana_funcs import cleanup

debug = False 

Nargs = len(sys.argv) 
if(Nargs<2): 
  print("[run_analysis]: Invalid input!")
  print("[run_analysis]: Example usage: python run_analysis.py config-filename") 
  sys.exit(0) 

# config file name 
configFile = sys.argv[1] 
full_path  = "./input/json/" + configFile 
# read the JSON file to get number of iterations 
inData   = json.loads(open(full_path).read())

NUM_EPOCH = inData["num-iter"]
scriptName = "TestCalib_batch.C"
for i in xrange(1,NUM_EPOCH+1): 
   cmd = "root -q -b -l '{0}+({1},\"{2}\")'".format(scriptName,i,full_path)
   if(debug):  print(cmd)
   if(not debug): os.system(cmd)

cleanup() 
