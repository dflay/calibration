# !/usr/bin/python 

# A script to run the calibration analysis

import os 
import sys
import json
import csv
import datetime

from myfuncs import *
from parse import StringParser

from ana_funcs import run_analysis_swapDBShim
from ana_funcs import run_analysis_imposedGrad_xy
from ana_funcs import run_analysis_imposedGrad_z
from ana_funcs import run_analysis_results
from ana_funcs import read_runlist
from ana_funcs import readJSONKey
from ana_funcs import readJSONSubKey
from ana_funcs import cleanup

debug = False 

Nargs = len(sys.argv) 
if(Nargs<2): 
  print("[print_runlist]: Invalid input!")
  print("[print_runlist]: Example usage: python print_runlist.py run-period 1,4,7-10") 
  sys.exit(0) 

# run period (1 or 2, for example) 
runPeriod  = int(sys.argv[1])   
# read JSON config file 
json_prefix = os.getcwd() + "/input/json/run-{0}".format(runPeriod) 
filepath    = "{0}/config.json".format(json_prefix)
confData    = json.loads( open(filepath).read() )
filePrefix  = confData["file-prefix"] 

# string of probes that we want 
probeStr   = sys.argv[2] 
# turn that into a list 
myParser = StringParser(probeStr)
myParser.GenerateList() 
probe = myParser.fList 

aList   = [] 
runList = [] 

thePath = "{0}/{1}".format(json_prefix,filePrefix)

# swapping, Delta-B, shimmed gradients  
for entry in probe:
   probeNumber = int(entry)
   aList = read_runlist(thePath,probeNumber)
   for j in aList: 
      runList.append( int(j) )
   print runList
   print("========================")  
   # cleanup 
   # aList[:] = []  

print(runList) 

outpath = "runlist_run-{0}.txt".format(runPeriod) 

with open(outpath,'w') as outfile: 
   for entry in runList: 
      outfile.write("%d\n" % entry)

print("[print_runlist]: The data has been written to the file: {0}".format(outpath) )  

# # cleanup output data
# print("[print_runlist]: Cleaning up output files... ")  
# homeDir = os.getcwd()  
# os.chdir(outDir)
# # make csv dir
# if not os.path.exists("./csv"):  
#    os.mkdir("csv",0755) 
#    print("[print_runlist]: Directory csv created")
# else:  
#    print("[print_runlist]: Directory csv already exists")
# # move files to csv directory  
# os.system("mv *.csv csv")
# # move back to home directory 
# os.chdir(homeDir)
# print("[print_runlist]: Done.")   
# 
# cleanup() 
