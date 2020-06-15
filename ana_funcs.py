# !/usr/bin/python 
# A script to run the calibration analysis

import os 
import sys
import json
import csv
import datetime

from myfuncs import * 

debug = False

#_______________________________________________________________________________
def cleanup(): 
   cmd = "rm *.d *.so *ACLiC*"
   os.system(cmd)
#_______________________________________________________________________________
def readJSONKey(filename,run_period,key): 
   # load data into JSON object  
   filepath = os.getcwd()+ "/input/json/run-{0}/{1}".format(run_period,filename)
   jData    = json.loads( open(filepath).read() )
   # return the value assocaited with the key 
   return jData[key] 
#_______________________________________________________________________________
def readJSONSubKey(filename,run_period,key,subkey): 
   # load data into JSON object  
   filepath = os.getcwd()+ "/input/json/run-{0}/{1}".format(run_period,filename)
   jData    = json.loads( open(filepath).read() )
   # return the value assocaited with the key and subkey  
   return jData[key][subkey] 
#_______________________________________________________________________________
def run_dB_analysis(deleteDir,run_period,calib_filename,trlyProbe): 
   filepath = os.getcwd()+ "/input/json/run-{0}/{1}".format(run_period,calib_filename) 
   inData   = json.loads( open(filepath).read() )
   # gather analysis parameters 
   runDate    = inData["date"]
   isBlind    = inData["blinding"]["enable"]
   blindLabel = inData["blinding"]["label"]
   # get today's date 
   today = datetime.datetime.today().strftime('%m-%d-%y')
   
   # delete the existing output data just in case it already exists 
   if(isBlind):
      outDir = "./output/blinded/" + blindLabel + "/" + today
   else: 
      outDir = "./output/unblinded/" + today
  
   if deleteDir:
      if os.path.exists(outDir):
         print( "Removing existing output from {0}".format(outDir) ) 
         deleteDir(outDir)  

   # set up the JSON directory to store config files 
   # that the ROOT macros use  
   json_prefix = "./input/json"
   json_prefix = json_prefix + "/ana" 
   if not os.path.exists(json_prefix):
       os.makedirs(json_prefix)

   # setup output directories 
   if(isBlind):
      outDir  = "./output/blinded/" + blindLabel + "/" + today
      plotDir = "./plots/blinded/" + blindLabel + "/" + today
   else: 
      outDir  = "./output/unblinded/" + today
      plotDir = "./plots/unblinded/"  + today

   if not os.path.exists(outDir):
       os.makedirs(outDir)
       print("[batchAnalysis]: Created directory: {0}".format(outDir) )

   if not os.path.exists(plotDir):
       os.makedirs(plotDir)
       print("[batchAnalysis]: Created directory: {0}".format(plotDir) )
   
   gradLabel = ["rad"       ,"vert"     ,"azi"]
   gradName  = ["norm-quad" ,"skew-quad","azi"]
   axisName  = ["x"         ,"y"        ,"z"]
   
   print("======================== DELTA-B ANALYSIS PARAMETERS ========================") 
   print("ana-date: {0}, blinding: {1}, trly-probe: {2}".format(runDate,isBlind,trlyProbe) )
   print("=====================================================================") 
      
   keyList = []

   #------------------------- run TRLY Delta-B calcs ----------------------------
   fitData       = False     # fitting the TRLY data? 
   tag           = "none"    # clear the top-level key for the JSON object  
   scriptName    = "DeltaB_trly_prod.C"
   for i in xrange(0,2):
      tag        = "dB-trly_{0}".format(axisName[i])
      configPath = "{0}/delta-b_trly-{1}.json".format(json_prefix,axisName[i])
      # make a run list  
      keyList[:] = []
      keyList.append("midas-runs")
      # write the input file 
      writeConfigFileProd_trlyDB(inData,tag,keyList,i,fitData,trlyProbe,configPath)
      cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,configPath)
      if(debug):  print(cmd)
      if(not debug): os.system(cmd)
   print("===============================================================")

   return 0
#_______________________________________________________________________________
def run_analysis_orig(deleteDir,configData,probeNumber):
   # original approach  
   json_prefix    = "./input/json"
   
   # calib_filename = "{0}_pr-{:02d}.json".format( fileNamePrefix,int(probeNumber) ) 
  
   # load data into JSON object  
   run_period     = configData["run-period"]  
   calib_filename = "{0}_pr-{1:02d}.json".format(configData["file-prefix"],probeNumber)  
   filepath       = os.getcwd() + "/input/json/run-{0}/{1}".format(run_period,calib_filename) 
   inData         = json.loads( open(filepath).read() )
   
   # gather analysis parameters
   isBlind    = configData["blinding"]["enable"]
   blindLabel = configData["blinding"]["label"]
   prodTag    = configData["prod-tag"]
   nmrAnaTag  = configData["nmr-ana-tag"] 
   runDate    = inData["date"]
   trlyProbe  = inData["trly-probe"] 
   
   # get today's date 
   today = datetime.datetime.today().strftime('%m-%d-%y')
   
   # delete the existing output data just in case it already exists 
   if(isBlind):
      outDir = "./output/blinded/" + blindLabel + "/" + today
   else: 
      outDir = "./output/unblinded/" + today
  
   if deleteDir:
      if os.path.exists(outDir):
         print( "Removing existing output from {0}".format(outDir) ) 
         deleteDir(outDir)  
   
   # set up the JSON directory to store config files 
   # that the ROOT macros use  
   json_prefix = json_prefix + "/ana"
   if not os.path.exists(json_prefix):
       os.makedirs(json_prefix)

   # setup output directories 
   if(isBlind):
      outDir  = "./output/blinded/" + blindLabel + "/" + today
      plotDir = "./plots/blinded/" + blindLabel + "/" + today
   else: 
      outDir  = "./output/unblinded/" + today
      plotDir = "./plots/unblinded/"  + today

   if not os.path.exists(outDir):
       os.makedirs(outDir)
       print("[batchAnalysis]: Created directory: {0}".format(outDir) )

   if not os.path.exists(plotDir):
       os.makedirs(plotDir)
       print("[batchAnalysis]: Created directory: {0}".format(plotDir) )
   
   gradLabel = ["rad"       ,"vert"     ,"azi"]
   gradName  = ["norm-quad" ,"skew-quad","azi"]
   axisName  = ["x"         ,"y"        ,"z"]
   
   print("======================== ANALYSIS PARAMETERS ========================") 
   print("run-date: {0}, blinding: {1}, trly-probe: {2}".format(runDate,isBlind,trlyProbe) )
   print("=====================================================================") 
      
   keyList = []
   
   #------------------------- run TRLY analysis code ---------------------------
   fitData       = False        # fitting the PP data?  
   tag           = "calib-run"  # top-level key for the JSON object 
   scriptName    = "Process_trly_prod.C"
   configPath    = "{0}/trly.json".format(json_prefix) 
   # keys needed for the analysis 
   keyList[:] = []
   keyList.append("midas-run")
   writeConfigFileProd(inData,configData,tag,keyList,-1,fitData,configPath) 
   cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,configPath)
   if(debug):  print(cmd)
   if(not debug): os.system(cmd)
   print("===============================================================")
   #------------------------- run PP analysis code ---------------------------
   fitData       = False        # fitting the PP data?  
   tag           = "calib-run"  # top-level key for the JSON object 
   scriptName    = "Process_pp_prod.C"
   configPath    = "{0}/pp.json".format(json_prefix) 
   # keys needed for the analysis 
   keyList[:] = []
   keyList.append("midas-run")
   writeConfigFileProd(inData,configData,tag,keyList,-1,fitData,configPath) 
   cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,configPath)
   if(debug):  print(cmd)
   if(not debug): os.system(cmd)
   print("===============================================================")
   #------------------------- run PP Delta-B calcs ------------------------------
   fitData       = False     # fitting the PP data? 
   tag           = "none"    # clear the top-level key for the JSON object  
   scriptName    = "DeltaB_pp_prod.C"
   for i in xrange(0,3):
      tag        = "dB-pp_{0}".format(axisName[i])
      configPath = "{0}/delta-b_pp-{1}.json".format(json_prefix,axisName[i])
      # make a run list  
      keyList[:] = []
      keyList.append("midas-runs")
      # write the input file 
      writeConfigFileProd(inData,configData,tag,keyList,i,fitData,configPath)
      cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,configPath)
      if(debug):  print(cmd)
      if(not debug): os.system(cmd)
   print("===============================================================")
   #------------------------- run TRLY Delta-B calcs ----------------------------
   fitData       = False     # fitting the TRLY data? 
   tag           = "none"    # clear the top-level key for the JSON object  
   scriptName    = "DeltaB_trly_prod.C"
   for i in xrange(0,3):
      tag        = "dB-trly_{0}".format(axisName[i])
      configPath = "{0}/delta-b_trly-{1}.json".format(json_prefix,axisName[i])
      # make a run list  
      keyList[:] = []
      keyList.append("midas-runs")
      # write the input file 
      writeConfigFileProd(inData,configData,tag,keyList,i,fitData,configPath)
      cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,configPath)
      if(debug):  print(cmd)
      if(not debug): os.system(cmd)
   print("===============================================================")
   #------------------------- run local shim scan calcs ---------------------
   # Shimmed gradient: PP scan analysis  
   fitData = True                 # fitting PP data? yes  
   tag     = "local-shim-scan"    # top-level JSON key 
   for i in xrange(0,3): 
      # set up JSON path 
      configPath = "{0}/local-shim-grad_pp{1}.json".format(json_prefix,axisName[i])
      # make a run list  
      keyList[:] = []
      keyList.append( "midas-runs"                )
      keyList.append( "pp{0}".format(axisName[i]) ) 
      writeConfigFileProd_ShimScan(inData,configData,tag,keyList,i,fitData,configPath)
      scriptName = "LocalScanGrad_pp_prod.C"
      cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,configPath)
      if(debug):  print(cmd)
      if(not debug): os.system(cmd)
   print("===============================================================")
   #------------------------- run imposed grad calcs ---------------------
   # Imposed gradients (xy)   
   fitData = False     # we are in fact fitting data, but will determine fit function on the fly   
   tag     = "NONE"    # top-level JSON key 
   # set up JSON path 
   configPath = "{0}/imposed-grad-xy.json".format(json_prefix,axisName[i])
   # make a run list  
   keyList[:] = []
   keyList.append("NONE")
   writeConfigFileProd(inData,configData,tag,keyList,-1,fitData,configPath)
   scriptName = "ImposedGrad_xy_prod.C"
   cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,configPath)
   if(debug):  print(cmd)
   if(not debug): os.system(cmd)
   print("===============================================================")
   #------------------------- run misalignment calcs ----------------------------
   fitData       = False
   tag           = "NONE"
   scriptName    = "Misalignment_prod.C"
   configPath    = "{0}/misalign.json".format(json_prefix) 
   keyList[:]    = []
   keyList.append("NONE") 
   writeConfigFileProd(inData,configData,tag,keyList,-1,fitData,configPath)
   cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,configPath)
   if(debug):  print(cmd)
   if(not debug): os.system(cmd)
   #------------------------- run calib analysis code ---------------------------
   print("===============================================================")
   fitData       = False
   tag           = "calib-run"
   scriptName    = "CalibrateMultiSwap_prod.C"
   configPath    = "{0}/calib.json".format(json_prefix) 
   # make a run list  
   keyList[:] = []
   keyList.append("midas-run")
   writeConfigFileProd(inData,configData,tag,keyList,-1,fitData,configPath) 
   cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,configPath)
   if(debug):  print(cmd)
   if(not debug): os.system(cmd)
   #------------------------- process all results ---------------------------
   print("===============================================================")
   fitData       = False
   tag           = "calib-run"
   scriptName    = "ProcessResults_prod.C"
   configPath    = "{0}/process-results.json".format(json_prefix) 
   # make a run list  
   keyList[:] = []
   keyList.append("midas-run")
   writeConfigFileProd(inData,configData,tag,keyList,-1,fitData,configPath) 
   cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,configPath)
   if(debug):  print(cmd)
   if(not debug): os.system(cmd)
   
   #------------------------- clean up files -------------------------
   # clean up
   cleanup()  

   return 0 
#_______________________________________________________________________________
def run_analysis_swapDB(configData,probeNumber): 
  
   # load data into JSON object  
   run_period     = configData["run-period"]  
   calib_filename = "{0}_pr-{1:02d}.json".format(configData["file-prefix"],probeNumber)  
   filepath       = os.getcwd() + "/input/json/run-{0}/{1}".format(run_period,calib_filename) 
   inData         = json.loads( open(filepath).read() )
   
   # gather analysis parameters
   isBlind      = configData["blinding"]["enable"]
   blindLabel   = configData["blinding"]["label"]
   prodTag      = configData["prod-tag"]
   nmrAnaTag    = configData["nmr-ana-tag"] 
   oscStatus    = configData["osc-cor"]["enable"] 
   oscType      = configData["osc-cor"]["type"] 
   runDate      = inData["date"]
   trlyProbe    = inData["trly-probe"] 
   zScanRunList = [] 
   zScanRunList = inData["trly_z-scan"]["midas-runs"] 

   useTRLYScan = False
   # FIXME: Need to iron out this code first...  
   # determine if we have a trly local z scan run to use
   # zRun = int(zScanRunList[0])  
   # if zRun>0:
   #   useTRLYScan = True 

   # set up the JSON directory to store config files 
   # that the ROOT macros use  
   json_prefix = "./input/json/ana"
   if not os.path.exists(json_prefix):
       os.makedirs(json_prefix)
      
   gradLabel = ["rad"       ,"vert"     ,"azi"]
   gradName  = ["norm-quad" ,"skew-quad","azi"]
   axisName  = ["x"         ,"y"        ,"z"]
   
   print("======================== ANALYSIS PARAMETERS ========================") 
   print("trolley probe:  {0}".format(trlyProbe) )
   print("run date:       {0}".format(runDate) )
   print("blinding:       {0}, enable = {1}".format(blindLabel,isBlind) )
   print("production tag: {0}".format(prodTag) )
   print("NMR-ANA tag:    {0}".format(nmrAnaTag) )
   print("osc cor:        {0}, enable = {1}".format(oscType,oscStatus) )
   if useTRLYScan: 
      print("TRLY local z scan: Using run {0}".format(zRun))
   print("=====================================================================") 
      
   keyList = []
 
   if(oscType=="standard"): 
      #------------------------- run TRLY analysis code ---------------------------
      fitData       = False        # fitting the PP data?  
      tag           = "calib-run"  # top-level key for the JSON object 
      scriptName    = "Process_trly_prod.C"
      configPath    = "{0}/trly.json".format(json_prefix) 
      # keys needed for the analysis 
      keyList[:] = []
      keyList.append("midas-run")
      writeConfigFileProd(inData,configData,tag,keyList,-1,fitData,configPath) 
      cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,configPath)
      if(debug):  print(cmd)
      if(not debug): os.system(cmd)
      print("===============================================================")
      #------------------------- run PP analysis code ---------------------------
      fitData       = False        # fitting the PP data?  
      tag           = "calib-run"  # top-level key for the JSON object 
      scriptName    = "Process_pp_prod.C"
      configPath    = "{0}/pp.json".format(json_prefix) 
      # keys needed for the analysis 
      keyList[:] = []
      keyList.append("midas-run")
      writeConfigFileProd(inData,configData,tag,keyList,-1,fitData,configPath) 
      cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,configPath)
      if(debug):  print(cmd)
      if(not debug): os.system(cmd)
      print("===============================================================")
   elif(oscType=="hybrid"): 
      #------------------------- run TRLY/PP hybrid analysis code ---------------------------
      # use the time of the first event of a given trial as the reference of the oscillation correction 
      fitData       = False        # fitting the PP data?  
      tag           = "calib-run"  # top-level key for the JSON object 
      scriptName    = "Process_pptr_prod.C"
      configPath    = "{0}/pp-trly-hybrid.json".format(json_prefix) 
      # keys needed for the analysis 
      keyList[:] = []
      keyList.append("midas-run")
      writeConfigFileProd(inData,configData,tag,keyList,-1,fitData,configPath) 
      cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,configPath)
      if(debug):  print(cmd)
      if(not debug): os.system(cmd)
      print("===============================================================")
   else: 
      print("[run_analysis_batch]: Invalid oscillation type {0}".format(oscType))

   #------------------------- run PP Delta-B calcs ------------------------------
   fitData       = False     # fitting the PP data? 
   tag           = "none"    # clear the top-level key for the JSON object  
   scriptName    = "DeltaB_pp_prod.C"
   for i in xrange(0,3):
      tag        = "dB-pp_{0}".format(axisName[i])
      configPath = "{0}/delta-b_pp-{1}.json".format(json_prefix,axisName[i])
      # make a run list  
      keyList[:] = []
      keyList.append("midas-runs")
      # write the input file 
      writeConfigFileProd(inData,configData,tag,keyList,i,fitData,configPath)
      cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,configPath)
      if(debug):  print(cmd)
      if(not debug): os.system(cmd)
   print("===============================================================")
   #------------------------- run TRLY Delta-B calcs ----------------------------
   fitData       = False     # fitting the TRLY data? 
   tag           = "none"    # clear the top-level key for the JSON object  
   scriptName    = "DeltaB_trly_prod.C"
   for i in xrange(0,3):
      tag        = "dB-trly_{0}".format(axisName[i])
      configPath = "{0}/delta-b_trly-{1}.json".format(json_prefix,axisName[i])
      # make a run list  
      keyList[:] = []
      keyList.append("midas-runs")
      # write the input file 
      writeConfigFileProd(inData,configData,tag,keyList,i,fitData,configPath)
      cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,configPath)
      if(debug):  print(cmd)
      if(not debug): os.system(cmd)
   print("===============================================================")
   #------------------------- clean up files -------------------------
   cleanup()  

   return 0
#_______________________________________________________________________________
def run_analysis_shimGrad(configData,probeNumber): 
  
   # load data into JSON object  
   run_period     = configData["run-period"]  
   calib_filename = "{0}_pr-{1:02d}.json".format(configData["file-prefix"],probeNumber)  
   filepath       = os.getcwd() + "/input/json/run-{0}/{1}".format(run_period,calib_filename) 
   inData         = json.loads( open(filepath).read() )
   
   # gather analysis parameters
   isBlind      = configData["blinding"]["enable"]
   blindLabel   = configData["blinding"]["label"]
   prodTag      = configData["prod-tag"]
   nmrAnaTag    = configData["nmr-ana-tag"] 
   oscStatus    = configData["osc-cor"]["enable"] 
   oscType      = configData["osc-cor"]["type"] 
   runDate      = inData["date"]
   trlyProbe    = inData["trly-probe"] 
   zScanRunList = [] 
   zScanRunList = inData["trly_z-scan"]["midas-runs"] 

   useTRLYScan = False
   # FIXME: Need to iron out this code first...  
   # determine if we have a trly local z scan run to use
   # zRun = int(zScanRunList[0])  
   # if zRun>0:
   #   useTRLYScan = True 

   # set up the JSON directory to store config files 
   # that the ROOT macros use  
   json_prefix = "./input/json/ana"
   if not os.path.exists(json_prefix):
       os.makedirs(json_prefix)
      
   gradLabel = ["rad"       ,"vert"     ,"azi"]
   gradName  = ["norm-quad" ,"skew-quad","azi"]
   axisName  = ["x"         ,"y"        ,"z"]
   
   # print("======================== ANALYSIS PARAMETERS ========================") 
   # print("trolley probe:  {0}".format(trlyProbe) )
   # print("run date:       {0}".format(runDate) )
   # print("blinding:       {0}, enable = {1}".format(blindLabel,isBlind) )
   # print("production tag: {0}".format(prodTag) )
   # print("NMR-ANA tag:    {0}".format(nmrAnaTag) )
   # print("osc cor:        {0}, enable = {1}".format(oscType,oscStatus) )
   # if useTRLYScan: 
   #    print("TRLY local z scan: Using run {0}".format(zRun))
   # print("=====================================================================") 
      
   keyList = []

   #------------------------- run local shim scan calcs ---------------------
   # Shimmed gradient: PP scan analysis 
   fitData = True                 # fitting PP data? yes  
   tag     = ""    # top-level JSON key 
   for i in xrange(0,3): 
      # set up JSON path 
      configPath = "{0}/local-shim-grad-{1}.json".format(json_prefix,axisName[i])
      # make a run list  
      keyList[:] = []
      keyList.append( "midas-runs" )
      if(useTRLYScan and i==2): 
         # need the trly scan  
         tag        = "trly_z-scan"     # top-level JSON key 
         scriptName = "LocalScanGrad_tr_prod.C"
      else:
         # using PP scan data  
         tag        = "local-shim-scan" # top-level JSON key 
         scriptName = "LocalScanGrad_pp_prod.C"
         keyList.append( "pp{0}".format(axisName[i]) ) 
      writeConfigFileProd_ShimScan(inData,configData,tag,keyList,i,fitData,configPath)
      cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,configPath)
      if(debug):  print(cmd)
      if(not debug): os.system(cmd)
   # fitData = True                 # fitting PP data? yes  
   # tag     = "local-shim-scan"    # top-level JSON key 
   # for i in xrange(0,3): 
   #    # set up JSON path 
   #    configPath = "{0}/local-shim-grad_pp{1}.json".format(json_prefix,axisName[i])
   #    # make a run list  
   #    keyList[:] = []
   #    keyList.append( "midas-runs"                )
   #    keyList.append( "pp{0}".format(axisName[i]) ) 
   #    writeConfigFileProd_ShimScan(inData,configData,tag,keyList,i,fitData,configPath)
   #    scriptName = "LocalScanGrad_pp_prod.C"
   #    cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,configPath)
   #    if(debug):  print(cmd)
   #    if(not debug): os.system(cmd)
   print("===============================================================")
   #------------------------- clean up files -------------------------
   cleanup()  

   return 0 

#_______________________________________________________________________________
def run_analysis_results(configData,probeNumber): 
   json_prefix    = "./input/json"
  
   # load data into JSON object  
   run_period     = configData["run-period"]  
   calib_filename = "{0}_pr-{1:02d}.json".format(configData["file-prefix"],probeNumber)  
   filepath       = os.getcwd() + "/input/json/run-{0}/{1}".format(run_period,calib_filename) 
   inData         = json.loads( open(filepath).read() )
   
   # set up the JSON directory to store config files 
   # that the ROOT macros use  
   runDate     = inData["date"]
   json_prefix = json_prefix + "/ana" 
   if not os.path.exists(json_prefix):
       os.makedirs(json_prefix)
   
   keyList = []

   print("======================== CALCULATING ANALYSIS RESULTS ========================") 
      
   #------------------------- run misalignment calcs ----------------------------
   fitData       = False
   tag           = "NONE"
   scriptName    = "Misalignment_prod.C"
   configPath    = "{0}/misalign.json".format(json_prefix) 
   keyList[:]    = []
   keyList.append("NONE") 
   writeConfigFileProd(inData,configData,tag,keyList,-1,fitData,configPath)
   cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,configPath)
   if(debug):  print(cmd)
   if(not debug): os.system(cmd)
   #------------------------- run calib analysis code ---------------------------
   print("===============================================================")
   fitData       = False
   tag           = "calib-run"
   scriptName    = "CalibrateMultiSwap_prod.C"
   configPath    = "{0}/calib.json".format(json_prefix) 
   # make a run list  
   keyList[:] = []
   keyList.append("midas-run")
   writeConfigFileProd(inData,configData,tag,keyList,-1,fitData,configPath) 
   cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,configPath)
   if(debug):  print(cmd)
   if(not debug): os.system(cmd)
   #------------------------- process all results ---------------------------
   print("===============================================================")
   fitData       = False
   tag           = "calib-run"
   scriptName    = "ProcessResults_prod.C"
   configPath    = "{0}/process-results.json".format(json_prefix) 
   # make a run list  
   keyList[:] = []
   keyList.append("midas-run")
   writeConfigFileProd(inData,configData,tag,keyList,-1,fitData,configPath) 
   cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,configPath)
   if(debug):  print(cmd)
   if(not debug): os.system(cmd)
   
   #------------------------- clean up files -------------------------
   # clean up
   cleanup()  

   return 0

#_______________________________________________________________________________
def run_analysis_imposedGrad_xy(configData):
   # calculate imposed gradients; we rely on previously processed dBq values from trolley  
   json_prefix    = "./input/json/imp-grad"
  
   # load data into JSON object  
   # run_period     = configData["run-period"]  
   # calib_filename = "{0}_pr-{1:02d}.json".format(configData["file-prefix"],probeNumber)  
   # filepath       = os.getcwd() + "/input/json/run-{0}/{1}".format(run_period,calib_filename) 
   # inData         = json.loads( open(filepath).read() )

   fitDim = configData["imp-grad"]["fit-dimension"] 
   
   # set up the JSON directory to store config files 
   # that the ROOT macros use  
   json_prefix = json_prefix
   if not os.path.exists(json_prefix):
       os.makedirs(json_prefix)
      
   keyList = []

   #------------------------- run imposed grad calcs ---------------------
   # Imposed gradients (xy)   
   fitData = False     # we are in fact fitting data, but will determine fit function on the fly   
   tag     = "NONE"    # top-level JSON key 
   # set up JSON path 
   configPath = "{0}/imposed-grad-xy.json".format(json_prefix)
   # make a run list  
   keyList[:] = []
   keyList.append("NONE")
   writeConfigFileProd_imposedGrad(configData,tag,keyList,-1,fitData,configPath)
   if(fitDim==1): 
      scriptName = "ImposedGrad_xy_prod.C"
   elif(fitDim==2): 
      scriptName = "ImposedGrad_xy_2D_prod.C"
   cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,configPath)
   if(debug):  print(cmd)
   if(not debug): os.system(cmd)
   print("===============================================================")
   # clean up
   cleanup()  

   return 0
#_______________________________________________________________________________
def run_analysis_imposedGrad_z(configData,probeNumber):
   # calculate imposed gradients; we rely on previously processed dBq values from trolley  
   json_prefix    = "./input/json/imp-grad"
  
   # load data into JSON object  
   run_period     = configData["run-period"]  
   calib_filename = "{0}_pr-{1:02d}.json".format(configData["file-prefix"],probeNumber)  
   filepath       = os.getcwd() + "/input/json/run-{0}/{1}".format(run_period,calib_filename) 
   inData         = json.loads( open(filepath).read() )
   
   # set up the JSON directory to store config files 
   # that the ROOT macros use  
   json_prefix = json_prefix
   if not os.path.exists(json_prefix):
       os.makedirs(json_prefix)
      
   keyList = []

   #------------------------- run imposed grad calcs ---------------------
   # Imposed gradients ()   
   fitData = False          # we are in fact fitting data, but will determine fit function on the fly   
   tag     = "imp-grad-z"   # top-level JSON key 
   # set up JSON path 
   configPath = "{0}/imposed-grad-z.json".format(json_prefix)
   # make a run list  
   keyList[:] = []
   keyList.append("midas-runs")
   writeConfigFileProd(inData,configData,tag,keyList,-1,fitData,configPath)
   scriptName = "ImposedGrad_z_prod.C"
   cmd = "root -q -b -l '{0}+(\"{1}\")'".format(scriptName,configPath)
   if(debug):  print(cmd)
   if(not debug): os.system(cmd)
   print("===============================================================")
   # clean up
   cleanup()  

   return 0
#_______________________________________________________________________________
def read_runlist(filePrefix,probeNumber): 
   # read all runs for a given probe
   # return a list of runs 
   inpath = "{0}-{1:02d}.json".format(filePrefix,probeNumber) 
   print("Reading data from: {0}".format(inpath)) 
   inData = json.loads( open(inpath).read() ) 

   keyList = [] 
   keyList.append("dB-pp_x") 
   keyList.append("dB-pp_y") 
   keyList.append("dB-pp_z") 
   if(probeNumber==1):
      # dB(x,y) for trolley is the same for all probes; do this once 
      keyList.append("dB-trly_x") 
      keyList.append("dB-trly_y") 
   keyList.append("dB-trly_z") 
   keyList.append("calib-run") 
   keyList.append("local-shim-scan")
   if(probeNumber==1): 
      keyList.append("trly_z-scan") 
      keyList.append("imp-grad-z") 

   subList = [] 
   myList = [] 
   theRun=0
   for entry in keyList: 
      if(entry=="calib-run"): 
         theRun = int(inData[entry]["midas-run"])   
         # make sure the run is valid  
         if(theRun>0): 
            myList.append(theRun)
      else: 
         subList = inData[entry]["midas-runs"]  
         for j in subList: 
            theRun = int(j) 
            if(theRun>0): 
               myList.append(j) 
      # cleanup 
      subList[:] = [] 

   return myList  

