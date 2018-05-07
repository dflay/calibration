# !/usr/bin/python 

# A script to run the calibration analysis

import os 
import sys
import json
import csv

from myfuncs import * 

debug = False

filename = os.getcwd()+"/input/json/delta-b_04-18-18.json"
inData   = json.loads(open(filename).read())

# gather analysis parameters 
theDate     = inData["date"]
probeNumber = inData["trly-probe"] - 1 # counting starts from zero   
isBlind     = inData["blinding"]

print("======================== ANALYSIS PARAMETERS ========================") 
print("date: {0}, blinding: {1}, trolley probe: {2}".format(theDate,isBlind,probeNumber) )
print("===============================================================") 

# get the fit functions we need 
fitFunc_trans = getFitFunction(inData,"trans-grad"  ,"fit"      )
fitFunc_azi   = getFitFunction(inData,"azi-grad"    ,"fit"      )

# prepare input files 
prepareRunListsDeltaBOnly(theDate,inData)  

isBlind = 0  
if( inData["blinding"] ):
   isBlind    = 1  
   print("*** Will blind the data ***") 

isFullAnalysis = 0  # false  
NAxes          = 2 

#------------------------- run PP analysis code -------------------------
# DeltaB (inital spot)  
# finalLocation = 0  
# scriptName    = "DeltaB_pp.C"
# for i in xrange(0,NAxes): 
#    cmd = "root -q -b -l '{0}+(\"{1}\",{2},{3},{4},{5})'".format(scriptName,theDate,i,finalLocation,isBlind,isFullAnalysis)
#    if(debug):  print(cmd)
#    if(not debug): os.system(cmd) 
# print("===============================================================") 
#------------------------- run TRLY analysis code -------------------------
# DeltaB (initial spot)  
scriptName = "DeltaB_trly.C"
finalLocation = 0  
for i in xrange(0,NAxes): 
   cmd = "root -q -b -l '{0}+(\"{1}\",{2},{3},{4},{5},{6})'".format(scriptName,theDate,probeNumber,i,finalLocation,isBlind,isFullAnalysis)
   if(debug):  print(cmd)
   if(not debug): os.system(cmd) 
print("===============================================================") 
# Find imposed gradients: transverse  
scriptName = "GetTransverseGradients.C" 
cmd = "root -q -b -l '{0}+(\"{1}\",\"{2}\",{3},{4},{5})'".format(scriptName,theDate,fitFunc_trans,probeNumber,isBlind,isFullAnalysis)
if(debug):  print(cmd)
if(not debug): os.system(cmd)
print("===============================================================") 
# Find imposed gradients: azimuthal  
# scriptName = "GetAziGradient.C" 
# cmd = "root -q -b -l '{0}+(\"{1}\",\"{2}\",{3},{4},{5})'".format(scriptName,theDate,fitFunc_azi,probeNumber,isBlind,isFullAnalysis)
# if(debug):  print(cmd)
# if(not debug): os.system(cmd)
# Now compute where to move the PP to align with the trolley  
# scriptName = "CalculatePPMovement.C"
# cmd = "root -q -b -l '{0}+(\"{1}\",{2})'".format(scriptName,theDate,probeNumber)
# if(debug):  print(cmd)
# if(not debug): os.system(cmd)
# print("===============================================================")
# clean up 
cmd = "python cleanup.py"
os.system(cmd) 

