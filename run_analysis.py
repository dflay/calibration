# !/usr/bin/python 

# A script to run the calibration analysis

import os 
import sys
import json
import csv

from myfuncs import * 

debug = False

json_prefix = "./input/json"

filename = os.getcwd()+"/input/json/calib_02-28-18.json"
# filename = os.getcwd()+"/input/json/calib_04-25-18.json"
inData   = json.loads(open(filename).read())

# gather analysis parameters 
theDate     = inData["date"]
probeNumber = inData["trly-probe"] # counting starts from zero   
isBlind     = inData["blinding"]

json_prefix = json_prefix + "/" + theDate
if not os.path.exists(json_prefix):
    os.makedirs(json_prefix) 

# a key list to use later 
keyList = []  
# other useful constants 
gradLabel = ["rad"       ,"vert"     ,"azi"]
gradName  = ["norm-quad" ,"skew-quad","azi"]
axisName  = ["x"         ,"y"        ,"z"]

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
if(ppScan): 
   print("Will use PP scan data!")
   fitFunc_ppScan_0   = getPrimaryValue(inData,"pp-scan-rad" ,"fit") 
   fitFunc_ppScan_1   = getPrimaryValue(inData,"pp-scan-vert","fit") 
   fitFunc_ppScan_2   = getPrimaryValue(inData,"pp-scan-azi" ,"fit") 
   fitFunc_ppScan     = [fitFunc_ppScan_0,fitFunc_ppScan_1,fitFunc_ppScan_2]

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

#------------------------- run PP analysis code -------------------------
if(not ppScan): 
   # Delta-B (PP)  
   finalLocation = 0         # shimmed location? no 
   fitData       = False     # fitting the PP data? 
   tag           = "none"    # clear the top-level key for the JSON object  
   scriptName    = "DeltaB_pp.C"
   for i in xrange(0,3):
      tag     = "dB-pp_{0}-grad".format(gradLabel[i])
      outpath = "{0}/delta-b_pp-{1}.json".format(json_prefix,axisName[i]) 
      # make a run list  
      keyList[:] = [] 
      keyList.append("bare") 
      keyList.append( "{0}".format(gradName[i]) )
      keyList.append("bare-2") 
      # write the input file 
      writeConfigFile(inData,tag,keyList,isFullAnalysis,finalLocation,i,fitData,outpath)
      cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,outpath)
      if(debug):  print(cmd)
      if(not debug): os.system(cmd)
   print("===============================================================")
   # do final location  
   finalLocation = 1         # shimmed location? yes 
   fitData       = False     # fitting the PP data? 
   tag           = "none"    # clear the top-level key for the JSON object  
   scriptName    = "DeltaB_pp.C"
   for i in xrange(0,3):
      tag     = "dB-pp_final-location_{0}-grad".format(gradLabel[i])
      outpath = "{0}/delta-b_pp-{1}.json".format(json_prefix,axisName[i]) 
      # make a run list  
      keyList[:] = [] 
      keyList.append("bare") 
      keyList.append( "{0}".format(gradName[i]) )
      keyList.append("bare-2") 
      # write the input file 
      writeConfigFile(inData,tag,keyList,isFullAnalysis,finalLocation,i,fitData,outpath)
      cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,outpath)
      if(debug):  print(cmd)
      if(not debug): os.system(cmd)
   print("===============================================================") 
else: 
   finalLocation = 0         # shimmed location? no 
   fitData       = True      # fitting the PP data? 
   tag           = "none"    # clear the top-level key for the JSON object  
   scriptName    = "GetImposedGradient_pp.C"
   for i in xrange(0,3):
      tag     = "pp-scan-{0}".format(gradLabel[i])
      outpath = "{0}/get-imposed-grad-pp-{1}".format(json_prefix,axisName[i]) 
      # make a run list  
      keyList[:] = [] 
      keyList.append("bare") 
      keyList.append( "{0}".format(gradName[i]) )
      # write the input file 
      writeConfigFile(inData,tag,keyList,isFullAnalysis,finalLocation,i,fitData,outpath)
      cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,outpath)
      if(debug):  print(cmd)
      if(not debug): os.system(cmd)
   print("===============================================================")
# Shimmed field 
finalLocation = 1            # shimmed location? yes 
fitData       = False        # fitting PP data? no  
tag           = "pp-shimmed" # top-level JSON key 
# set up JSON path 
outpath = "{0}/pp-shimmed.json".format(json_prefix) 
# make a run list 
keyList[:] = []         # clear entries  
keyList.append("shim")  # run list 
writeConfigFile(inData,tag,keyList,isFullAnalysis,finalLocation,-1,fitData,outpath)
scriptName = "GetShimmedField_pp.C"
cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,outpath)
if(debug):  print(cmd)
if(not debug): os.system(cmd) 
print("===============================================================") 
#------------------------- run TRLY analysis code -------------------------
# Delta-B (TRLY)  
finalLocation = 0         # shimmed location? no 
fitData       = False     # fitting the PP data? 
tag           = "none"    # clear the top-level key for the JSON object  
scriptName    = "DeltaB_trly.C"
for i in xrange(0,3):
   tag     = "dB-trly_{0}-grad".format(gradLabel[i])
   outpath = "{0}/delta-b_trly-{1}.json".format(json_prefix,axisName[i]) 
   # make a run list  
   keyList[:] = [] 
   keyList.append("bare") 
   keyList.append( "{0}".format(gradName[i]) )
   keyList.append("bare-2") 
   # write the input file 
   writeConfigFile(inData,tag,keyList,isFullAnalysis,finalLocation,i,fitData,outpath)
   cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,outpath)
   if(debug):  print(cmd)
   if(not debug): os.system(cmd) 
print("===============================================================") 
# shimmed location   
finalLocation = 1         # shimmed location? no 
fitData       = False     # fitting the PP data? 
tag           = "none"    # clear the top-level key for the JSON object  
scriptName    = "DeltaB_trly.C"
for i in xrange(0,3):
   tag     = "dB-trly_final-location_{0}-grad".format(gradLabel[i])
   outpath = "{0}/delta-b_final-location_trly-{1}.json".format(json_prefix,axisName[i]) 
   # make a run list  
   keyList[:] = [] 
   keyList.append("bare") 
   keyList.append( "{0}".format(gradName[i]) )
   keyList.append("bare-2") 
   # write the input file 
   writeConfigFile(inData,tag,keyList,isFullAnalysis,finalLocation,i,fitData,outpath)
   cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,outpath)
   if(debug):  print(cmd)
   if(not debug): os.system(cmd) 
print("===============================================================") 
# Shimmed field (TRLY) 
finalLocation = 1              # shimmed location? yes 
fitData       = False          # fitting PP data? no  
tag           = "trly-shimmed" # top-level JSON key 
# set up JSON path 
outpath = "{0}/trly-shimmed.json".format(json_prefix) 
# make a run list 
keyList[:] = []         # clear entries  
keyList.append("shim")  # run list 
writeConfigFile(inData,tag,keyList,isFullAnalysis,finalLocation,-1,fitData,outpath)
scriptName = "GetShimmedField_trly.C"
cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,outpath)
if(debug):  print(cmd)
if(not debug): os.system(cmd)
print("===============================================================")
sys.exit(0)
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
   finalLocation = 1                 # shimmed location? yes 
   fitData       = True              # fitting PP data? yes  
   tag           = "trly-shimmed-grad_trans" # top-level JSON key 
   # set up JSON path 
   outpath = "{0}/get-shim-trans-grad-trly.json".format(json_prefix) 
   # make a run list 
   keyList[:] = []         # clear entries  
   keyList.append("trly")  # run list 
   writeConfigFile(inData,tag,keyList,isFullAnalysis,finalLocation,-1,fitData,outpath)
   scriptName = "GetShimmedTransGrad.C" 
   cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,outpath)
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
   # PP shimmed gradients  
   finalLocation = 1                 # shimmed location? yes 
   fitData       = True              # fitting PP data? yes  
   tag           = "pp-shimmed-grad" # top-level JSON key 
   for i in xrange(0,2):
      # set up JSON path 
      outpath = "{0}/get-shim-grad-pp-{1}.json".format(json_prefix,axisName[i]) 
      # make a run list 
      keyList[:] = [] # clear entries  
      keyList.append("pp")
      keyList.append( "pp{0}".format(axisName[i]) ) 
      keyList.append("trly") 
      writeConfigFile(inData,tag,keyList,isFullAnalysis,finalLocation,i,fitData,outpath)
      scriptName = "GetShimmedGrad_pp.C" 
      cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,outpath)
      if(debug):  print(cmd)
      if(not debug): os.system(cmd)
   print("===============================================================") 
   # TRLY shimmed gradients: transverse 
   finalLocation = 1                 # shimmed location? yes 
   fitData       = True              # fitting PP data? yes  
   tag           = "trly-shimmed-grad_trans" # top-level JSON key 
   # set up JSON path 
   outpath = "{0}/get-shim-trans-grad-trly.json".format(json_prefix) 
   # make a run list 
   keyList[:] = []         # clear entries  
   keyList.append("trly")  # run list 
   writeConfigFile(inData,tag,keyList,isFullAnalysis,finalLocation,-1,fitData,outpath)
   scriptName = "GetShimmedTransGrad.C" 
   cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,outpath)
   if(debug):  print(cmd)
   if(not debug): os.system(cmd)
   print("===============================================================") 
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
cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,outpath)
if(debug):  print(cmd)
if(not debug): os.system(cmd)
print("===============================================================")
#------------------------- clean up files -------------------------
# clean up 
cmd = "python cleanup.py"
os.system(cmd) 

