# !/usr/bin/python 

# A script to run the calibration analysis

import os 
import sys
import json
import csv
import datetime

from myfuncs import *
from parse import StringParser

from ana_funcs import run_analysis_swapDB
from ana_funcs import run_analysis_shimGrad
from ana_funcs import run_analysis_imposedGrad_xy
from ana_funcs import run_analysis_imposedGrad_z
from ana_funcs import run_analysis_results
from ana_funcs import readJSONKey
from ana_funcs import readJSONSubKey
from ana_funcs import cleanup

debug       = False
cleanOutput = True  

Nargs = len(sys.argv) 
if(Nargs<2): 
  print("[run_analysis_batch]: Invalid input!")
  print("[run_analysis_batch]: Example usage: python run_analysis_batch.py run-period 1,4,7-10") 
  sys.exit(0) 

startTime = utcNow() 

# run period (1 or 2, for example) 
runPeriod  = int(sys.argv[1])   
# read JSON config file 
filepath = os.getcwd()+ "/input/json/run-{0}/config.json".format(runPeriod)
confData = json.loads( open(filepath).read() )

# string of probes that we want 
probeStr   = sys.argv[2] 
# turn that into a list 
myParser = StringParser(probeStr)
myParser.GenerateList() 
probe = myParser.fList 
 
deleteOutDir = False 

# gather analysis parameters
isBlind    = confData["blinding"]["enable"]
isSyst     = confData["syst"]["enable"]
blindLabel = confData["blinding"]["label"]
prodTag    = confData["prod-tag"]
nmrAnaTag  = confData["nmr-ana-tag"]

# get today's date 
today = datetime.datetime.today().strftime('%m-%d-%y')

outDir,plotDir,systDirNum = setupPaths(deleteOutDir,today,isBlind,blindLabel,isSyst)

# update the json object with the systematic directory number
# if running in systematic mode, we need to track this directory for storing output   
confData["syst"]["dir-num"] = int(systDirNum)
confData["ana-date"]        = today  

if(isBlind==1):
   print("[run_analysis_batch]: Using blind: {0}".format(blindLabel))

# swapping, Delta-B, shimmed gradients  
for entry in probe:
   probeNumber = int(entry) 
   run_analysis_swapDB(confData,probeNumber)
# need all TRLY dB processed before getting imposed gradients 
run_analysis_imposedGrad_xy(confData)
# dB/dz 
for entry in probe: 
   probeNumber = int(entry) 
   run_analysis_imposedGrad_z(confData,probeNumber) 
# shimmed gradients; need imposed gradients to be finished first in case we use alt method   
for entry in probe:  
   probeNumber = int(entry) 
   run_analysis_shimGrad(confData,probeNumber)
# now calculate all results 
for entry in probe:
   probeNumber = int(entry) 
   run_analysis_results(confData,probeNumber)

# now print some tables and collect everything into single output files
theDate    = datetime.datetime.today().strftime('%m-%d-%y')
scriptName = "MakeTables_prod.C"
cmd = "root -q -b -l '{0}+({1},\"{2}\",{3},{4})'".format(scriptName,runPeriod,theDate,int(isSyst),systDirNum)
if(debug):  print(cmd)
if(not debug): os.system(cmd)

# print config file to output to keep track of things 
outpath = outDir + "/config.json"
writeConfigFileProd_params(confData,outpath)

if(cleanOutput): 
   # cleanup output data
   print("[run_analysis_batch]: Cleaning up output files... ")  
   homeDir = os.getcwd()  
   os.chdir(outDir)
   # make csv dir
   if not os.path.exists("./csv"):  
      os.mkdir("csv",0755) 
      print("[run_analysis_batch]: Directory csv created")
   else:  
      print("[run_analysis_batch]: Directory csv already exists")
   # move files to csv directory  
   os.system("mv *.csv csv")
   os.system("mv csv/calibData* .")  # retrieve the csv output
   # move back to home directory 
   os.chdir(homeDir)
   print("[run_analysis_batch]: Done.")   

endTime = utcNow() 

elapsedTime     = float(endTime-startTime) 
elapsedTime_min = elapsedTime/60. 

print("[run_analysis_batch]: Elapsed time = {0:.3f} sec ({1:0.3f} min)".format(elapsedTime,elapsedTime_min))

# remove compilation remnants from working directory  
cleanup() 
