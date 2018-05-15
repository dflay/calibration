# !/usr/bin/python 

# A script to run the calibration analysis

import os 
import sys
import json
import csv

from myfuncs import * 

debug = False

json_prefix = "./input/json"

Nargs = len(sys.argv) 
if(Nargs!=2): 
  print("[run_analysis]: Invalid input!")
  print("[run_analysis]: Usage: python {0} calib-filename".format(sys.argv[0])) 
  sys.exit(0) 

calib_filename = sys.argv[1] 

filename = os.getcwd()+ "/input/json/{0}".format(calib_filename) 
# filename = os.getcwd()+ "/input/json/calib-alt_02-28-18.json"
# filename = os.getcwd()+"/input/json/calib_04-25-18.json"
inData   = json.loads(open(filename).read())

# gather analysis parameters 
theDate     = inData["date"]
probeNumber = inData["trly-probe"] # counting starts from zero   
isBlind     = inData["blinding"]
useP2PFit   = inData["p2p-fit"]

json_prefix = json_prefix + "/" + theDate
if not os.path.exists(json_prefix):
    os.makedirs(json_prefix) 

# a key list to use later 
keyList = []  
# useful constants 
gradLabel = ["rad"       ,"vert"     ,"azi"]
gradName  = ["norm-quad" ,"skew-quad","azi"]
axisName  = ["x"         ,"y"        ,"z"]

print("======================== ANALYSIS PARAMETERS ========================") 
print("date: {0}, blinding: {1}, trolley probe: {2}".format(theDate,isBlind,probeNumber) )
print("===============================================================") 

# PP scan info 
ppScan             = getPrimaryValue(inData,"pp-scan"     ,"enable") 
if(ppScan): 
   print("--> Will use PP scan data!")

# prepare TRLY mode input files 
trly_mode_prefix = "./input/runlists/{0}".format(theDate) 
deleteDir(trly_mode_prefix) # delete existing directory 
createDir(trly_mode_prefix) # recreate directory 
writeToFileTRLYMode(trly_mode_prefix,"interactive",inData)
writeToFileTRLYMode(trly_mode_prefix,"continuous" ,inData) 

isBlind = 0  
if( inData["blinding"] ):
   isBlind    = 1  
   print("*** Will blind the data ***") 

isFullAnalysis = 1 # true  

#------------------------- run PP analysis code -------------------------
# # Delta-B (PP)  
# finalLocation = 0         # shimmed location? no 
# fitData       = False     # fitting the PP data? 
# tag           = "none"    # clear the top-level key for the JSON object  
# scriptName    = "DeltaB_pp.C"
# for i in xrange(0,3):
#    tag     = "dB-pp_{0}-grad".format(gradLabel[i])
#    configPath = "{0}/delta-b_pp-{1}.json".format(json_prefix,axisName[i]) 
#    # make a run list  
#    keyList[:] = [] 
#    keyList.append("bare") 
#    keyList.append( "{0}".format(gradName[i]) )
#    keyList.append("bare-2") 
#    # write the input file 
#    writeConfigFile(inData,tag,keyList,isFullAnalysis,finalLocation,i,fitData,configPath)
#    cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,configPath)
#    if(debug):  print(cmd)
#    if(not debug): os.system(cmd)
# print("===============================================================")
# do final location  
finalLocation = 1         # shimmed location? yes 
fitData       = False     # fitting the PP data? 
tag           = "none"    # clear the top-level key for the JSON object  
scriptName    = "DeltaB_pp.C"
for i in xrange(0,3):
   tag     = "dB-pp_final-location_{0}-grad".format(gradLabel[i])
   configPath = "{0}/delta-b_pp-{1}.json".format(json_prefix,axisName[i]) 
   # make a run list  
   keyList[:] = [] 
   keyList.append("bare") 
   keyList.append( "{0}".format(gradName[i]) )
   keyList.append("bare-2") 
   # write the input file 
   writeConfigFile(inData,tag,keyList,isFullAnalysis,finalLocation,i,fitData,configPath)
   cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,configPath)
   if(debug):  print(cmd)
   if(not debug): os.system(cmd)
print("===============================================================") 
#------------------------- run TRLY analysis code -------------------------
# # Delta-B (TRLY)  
# finalLocation = 0         # shimmed location? no 
# fitData       = False     # fitting the PP data? 
# tag           = "none"    # clear the top-level key for the JSON object  
# scriptName    = "DeltaB_trly.C"
# for i in xrange(0,3):
#    tag     = "dB-trly_{0}-grad".format(gradLabel[i])
#    configPath = "{0}/delta-b_trly-{1}.json".format(json_prefix,axisName[i]) 
#    # make a run list  
#    keyList[:] = [] 
#    keyList.append("bare") 
#    keyList.append( "{0}".format(gradName[i]) )
#    keyList.append("bare-2") 
#    # write the input file 
#    writeConfigFile(inData,tag,keyList,isFullAnalysis,finalLocation,i,fitData,configPath)
#    cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,configPath)
#    if(debug):  print(cmd)
#    if(not debug): os.system(cmd) 
# print("===============================================================") 
# shimmed location   
finalLocation = 1         # shimmed location? no 
fitData       = False     # fitting the PP data? 
tag           = "none"    # clear the top-level key for the JSON object  
scriptName    = "DeltaB_trly.C"
for i in xrange(0,3):
   tag     = "dB-trly_final-location_{0}-grad".format(gradLabel[i])
   configPath = "{0}/delta-b_final-location_trly-{1}.json".format(json_prefix,axisName[i]) 
   # make a run list  
   keyList[:] = [] 
   keyList.append("bare") 
   keyList.append( "{0}".format(gradName[i]) )
   keyList.append("bare-2") 
   # write the input file 
   writeConfigFile(inData,tag,keyList,isFullAnalysis,finalLocation,i,fitData,configPath)
   cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,configPath)
   if(debug):  print(cmd)
   if(not debug): os.system(cmd) 
print("===============================================================") 
#------------------------- run gradient analysis code -------------------------
# Find imposed gradients: transverse 
finalLocation = 0              # shimmed location? no 
fitData       = True           # fitting TRLY data? yes  
tag           = "trans-grad"   # top-level JSON key 
# set up JSON path 
configPath = "{0}/get-trans-grad-trly.json".format(json_prefix) 
# make a run list 
keyList[:] = []                # clear entries  
keyList    = ["bare","norm-quad","skew-quad","bare-2"]  # run list 
writeConfigFile(inData,tag,keyList,isFullAnalysis,finalLocation,-1,fitData,configPath)   
scriptName = "GetTransverseGradients.C" 
cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,configPath)
if(debug):  print(cmd)
if(not debug): os.system(cmd)
print("===============================================================") 
# Find imposed gradients: azimuthal  
# finalLocation = 0              # shimmed location? no 
# fitData       = True           # fitting TRLY data? yes  
# tag           = "azi-grad"    # top-level JSON key 
# # set up JSON path 
# configPath = "{0}/get-azi-grad-trly.json".format(json_prefix) 
# # make a run list 
# keyList[:] = []                # clear entries  
# keyList    = ["bare","azi","bare-2"]  # run list 
# writeConfigFile(inData,tag,keyList,isFullAnalysis,finalLocation,-1,fitData,configPath)
# scriptName = "GetAziGradient.C" 
# cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,configPath)
# if(debug):  print(cmd)
# if(not debug): os.system(cmd)
# print("===============================================================") 
#------------------------- clean up files -------------------------
# clean up 
cmd = "python cleanup.py"
os.system(cmd) 

