# !/usr/bin/python 

# A script to run the calibration analysis

import os 
import sys
import json
import csv

from myfuncs import * 

debug = False

json_prefix = "./input/json"

# filename = os.getcwd()+"/input/json/calib_02-28-18.json"
filename = os.getcwd()+"/input/json/calib_04-25-18.json"
inData   = json.loads(open(filename).read())

# gather analysis parameters 
theDate     = inData["date"]
probeNumber = inData["trly-probe"] # counting starts from zero   
isBlind     = inData["blinding"]

print("======================== ANALYSIS PARAMETERS ========================") 
print("date: {0}, blinding: {1}, trolley probe: {2}".format(theDate,isBlind,probeNumber) )
print("===============================================================") 

# get the fit functions we need 
fitFunc_shim_trans = getPrimaryValue(inData,"shimmed-grad","trans-fit")
fitFunc_shim_azi   = getPrimaryValue(inData,"shimmed-grad","azi-fit"  )
fitFunc_trans      = getPrimaryValue(inData,"trans-grad"  ,"fit"      )
fitFunc_azi        = getPrimaryValue(inData,"azi-grad"    ,"fit"      )

# get P2P fit status 
useP2PFit          = inData["p2p-fit"]

# PP scan info 
ppScan             = getPrimaryValue(inData,"pp-scan"     ,"enable") 
fitFunc_ppScan_0   = getPrimaryValue(inData,"pp-scan-rad" ,"fit") 
fitFunc_ppScan_1   = getPrimaryValue(inData,"pp-scan-vert","fit") 
fitFunc_ppScan_2   = getPrimaryValue(inData,"pp-scan-azi" ,"fit") 

fitFunc_ppScan     = [fitFunc_ppScan_0,fitFunc_ppScan_1,fitFunc_ppScan_2] 

if(ppScan): 
   print("Will use PP scan data!")

p2pFitStatus = 0 # assume false 
if( useP2PFit ) :
   p2pFitStatus = 1  

# get the shimmed runs 
shimRun_pp   = getShimRun(inData,"pp-shimmed"  ) 
shimRun_trly = getShimRun(inData,"trly-shimmed") 

# prepare input files 
prepareRunLists(theDate,inData,ppScan)  

isBlind = 0  
if( inData["blinding"] ):
   isBlind    = 1  
   print("*** Will blind the data ***") 

isFullAnalysis = 1 # true  
# writeConfigFile(inData,"pp-shimmed-grad",keyList,isFullAnalysis,0,0,"test.json") 

#------------------------- run PP analysis code -------------------------
if(not ppScan):  
   # DeltaB (inital spot) 
   finalLocation = 0  
   scriptName    = "DeltaB_pp.C"
   for i in xrange(0,3): 
      cmd = "root -q -b -l '{0}+(\"{1}\",{2},{3},{4},{5},{6})'".format(scriptName,theDate,i,finalLocation,isBlind,isFullAnalysis,p2pFitStatus)
      if(debug):  print(cmd)
      if(not debug): os.system(cmd) 
   print("===============================================================") 
   # DeltaB (final spot)  
   finalLocation = 1  
   for i in xrange(0,3): 
      cmd = "root -q -b -l '{0}+(\"{1}\",{2},{3},{4},{5},{6})'".format(scriptName,theDate,i,finalLocation,isBlind,isFullAnalysis,p2pFitStatus)
      if(debug):  print(cmd)
      if(not debug): os.system(cmd) 
   print("===============================================================") 
else: 
   scriptName = "GetImposedGradient_pp.C"
   for i in xrange(0,3): 
      cmd = "root -q -b -l '{0}+(\"{1}\",\"{2}\",{3},{4},{5},{6},{7})'".format(scriptName,theDate,fitFunc_ppScan[i],probeNumber,i,isBlind,isFullAnalysis,p2pFitStatus)
      if(debug):  print(cmd)
      if(not debug): os.system(cmd)
   print("===============================================================") 
# Shimmed field 
scriptName = "GetShimmedField_pp.C"
cmd = "root -q -b -l '{0}+(\"{1}\",{2},{3})'".format(scriptName,theDate,isBlind,p2pFitStatus)
if(debug):  print(cmd)
if(not debug): os.system(cmd) 
print("===============================================================") 
#------------------------- run TRLY analysis code -------------------------
# DeltaB (initial spot)  
scriptName = "DeltaB_trly.C"
finalLocation = 0  
for i in xrange(0,3): 
   cmd = "root -q -b -l '{0}+(\"{1}\",{2},{3},{4},{5},{6},{7})'".format(scriptName,theDate,probeNumber,i,finalLocation,isBlind,isFullAnalysis,p2pFitStatus)
   if(debug):  print(cmd)
   if(not debug): os.system(cmd) 
print("===============================================================") 
# DeltaB (final spot)  
finalLocation = 1  
for i in xrange(0,3): 
   cmd = "root -q -b -l '{0}+(\"{1}\",{2},{3},{4},{5},{6},{7})'".format(scriptName,theDate,probeNumber,i,finalLocation,isBlind,isFullAnalysis,p2pFitStatus)
   if(debug):  print(cmd)
   if(not debug): os.system(cmd)
print("===============================================================") 
# Shimmed field (TRLY) 
scriptName = "GetShimmedField_trly.C"
cmd = "root -q -b -l '{0}+(\"{1}\",{2},{3},{4})'".format(scriptName,theDate,probeNumber,isBlind,p2pFitStatus)
if(debug):  print(cmd)
if(not debug): os.system(cmd)
print("===============================================================")
if(not ppScan): 
   #------------------------- run gradient analysis code -------------------------
   # Find imposed gradients: transverse  
   scriptName = "GetTransverseGradients.C" 
   cmd = "root -q -b -l '{0}+(\"{1}\",\"{2}\",{3},{4},{5})'".format(scriptName,theDate,fitFunc_trans,probeNumber,isBlind,isFullAnalysis)
   if(debug):  print(cmd)
   if(not debug): os.system(cmd)
   print("===============================================================") 
   # Find imposed gradients: azimuthal  
   scriptName = "GetAziGradient.C" 
   cmd = "root -q -b -l '{0}+(\"{1}\",\"{2}\",{3},{4},{5})'".format(scriptName,theDate,fitFunc_azi,probeNumber,isBlind,isFullAnalysis)
   if(debug):  print(cmd)
   if(not debug): os.system(cmd)
   print("===============================================================") 
   # Shimmed gradients: transverse 
   scriptName = "GetShimmedTransGrad.C" 
   cmd = "root -q -b -l '{0}+(\"{1}\",\"{2}\",{3},{4})'".format(scriptName,theDate,fitFunc_shim_trans,probeNumber,isBlind)
   if(debug):  print(cmd)
   if(not debug): os.system(cmd)
   print("===============================================================") 
   # Shimmed gradients: azimuthal  
   scriptName = "GetShimmedAziGrad.C" 
   cmd = "root -q -b -l '{0}+(\"{1}\",\"{2}\",{3},{4})'".format(scriptName,theDate,fitFunc_shim_azi,probeNumber,isBlind)
   if(debug):  print(cmd)
   if(not debug): os.system(cmd)
   print("===============================================================") 
else:  
   outpath = json_prefix + "/get-shim-grad-pp.json"
   finalLocation = 1
   fitData = True 
   keyList = []  
   for i in xrange(0,2):
      keyList[:] = [] # clear entries  
      if(i==0):
         keyList = ["pp","ppx","trly"]
      elif(i==1): 
         keyList = ["pp","ppy","trly"]
      elif(i==2): 
         keyList = ["pp","ppz","trly"]
      writeConfigFile(inData,"pp-shimmed-grad",keyList,isFullAnalysis,finalLocation,i,fitData,outpath)
      scriptName = "GetShimmedGrad_pp.C" 
      cmd = "root -q -b -l '{0}+()'".format(scriptName)
      if(debug):  print(cmd)
      if(not debug): os.system(cmd)
   # sys.exit(0) 
#------------------------- run calibration analysis code -------------------------
# Now do the calibration calculation 
keyList[:] = [] 
keyList = ["trly","pp"]
outpath = json_prefix + "/calibrate.json"
fitData = False
finalLocation = 1 
writeConfigFile(inData,"calib-runs",keyList,isFullAnalysis,finalLocation,-1,fitData,outpath)
sys.exit(0) 
scriptName = "Calibrate.C"
cmd = "root -q -b -l '{0}+()'".format(scriptName)
if(debug):  print(cmd)
if(not debug): os.system(cmd)
print("===============================================================")
#------------------------- clean up files -------------------------
# clean up 
cmd = "python cleanup.py"
os.system(cmd) 

