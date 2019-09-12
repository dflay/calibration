# !/usr/bin/python 
# A script that tests out JSON file reading  

import os 
import sys
import json
import csv
import shutil

#_______________________________________________________________________________
def writeConfigFile(data,tag,keyList,isFullAnalysis,isFinalLocation,axis,fitData,outpath): 
   # write a JSON file for the ROOT script defined by tag 
   # get the run list 
   runList   = []
   labelList = []  
   subList   = []
   i=0
   for key in keyList: 
      if(key=="ppx" or key=="ppy" or key=="ppz"):
         # for PP NMR-DAQ runs  
         subList = data[tag][key] # this is a list!  
         for subRun in subList: 
            runList.append(subRun) 
            labelList.append("sr{0}".format(i+1))
            i = i + 1 
      else: 
         # regular MIDAS runs 
         runList.append(data[tag][key])   
         labelList.append(key)  
   outData = {} 
   outData['date']       = data['date']
   outData['type']       = data['type'] 
   outData['blinding']   = data['blinding'] 
   outData['trly-probe'] = data['trly-probe'] 
   outData['p2p-fit']    = data['p2p-fit'] 
   if(fitData): 
      outData['fit'] = data[tag]['fit']  
   else: 
      outData['fit'] = "NONE"
   outData['full-ana']   = isFullAnalysis 
   outData['final-loc']  = isFinalLocation 
   outData['nruns']      = len(runList) 
   outData['run-list']   = runList   
   outData['run-label']  = labelList  
   outData['fxpr-set']   = data['fxpr-set'] 
   outData['axis']       = axis 
   # print("{0}".format(outData) )
   outfile = open(outpath,'w')
   json.dump(outData,outfile)
   outfile.close()
#_______________________________________________________________________________
def writeConfigFileProd_trlyDB(data,tag,keyList,axis,fitData,trlyProbe,outpath): 
   # write a JSON file for the ROOT script defined by tag 
   # get the run list 
   runList   = []
   labelList = []  
   subList   = []
   i=0
   for key in keyList: 
      if(key=="ppx" or key=="ppy" or key=="ppz" or key=="midas-runs"):
         # for PP NMR-DAQ runs  
         subList = data[tag][key] # this is a list!  
         for subRun in subList: 
            runList.append(subRun) 
            labelList.append("sr{0}".format(i+1))
            i = i + 1 
      elif(key=="NONE"): 
         print("No runs to use!") 
      else: 
         # single MIDAS runs 
         runList.append(data[tag][key])   
         labelList.append(key) 
   outData = {} 
   # config file data 
   outData['type']                   = confData['type'] 
   outData['blinding']               = confData['blinding'] 
   outData['run-period']             = confData['run-period'] 
   outData['prod-tag']               = confData['prod-tag'] 
   outData['nmr-ana-tag']            = confData['nmr-ana-tag']
   outData['fxpr-set']               = confData['fxpr-set']
   outData['free-proton-cor']        = confData['free-proton-cor'] 
   outData['load-trly-swap-times']   = confData['load-trly-swap-times'] 
   outData['use-aba-time-weight']    = confData['use-aba-time-weight']
   outData['use-trly-temp-cor']      = confData['use-trly-temp-cor'] 
   outData['use-osc-cor']            = confData['use-osc-cor'] 
   outData['cut-file']               = confData['cut-file'] 
   outData['num-events-to-avg']      = confData['num-events-to-avg']  
   outData['num-events-time-window'] = confData['num-events-time-window'] 
   outData['fxpr-remove-drift']      = confData['fxpr-remove-drift']  
   outData['use-misalign-cor']       = confData['use-misalign-cor']  
    
   # probe-specific data  
   outData['date']                   = data['date']
   outData['trly-probe']             = int(trlyProbe) 

   if(fitData): 
      outData['fit'] = data[tag]['fit']  
   else: 
      outData['fit'] = "NONE"

   outData['nruns']      = len(runList) 
   outData['run-list']   = runList   
   outData['run-label']  = labelList  
   outData['axis']       = axis 
   # print("{0}".format(outData) )
   outfile = open(outpath,'w')
   json.dump(outData,outfile)
   outfile.close()
#_______________________________________________________________________________
def writeConfigFileProd_ShimScan(data,confData,tag,keyList,axis,fitData,outpath):
   # special setup because sometimes we have multiple MIDAS runs for a scan...  
   # write a JSON file for the ROOT script defined by tag 
   # get the run list 
   runList   = []
   labelList = []  
   subList   = []
   i=0
   for key in keyList: 
      if(key=="ppx" or key=="ppy" or key=="ppz" or key=="midas-runs"):
         # for a list of PP NMR-DAQ runs or MIDAS runs  
         subList = data[tag][key] # this is a list!  
         for subRun in subList: 
            runList.append(subRun) 
            if(key=="midas-runs"):
                labelList.append("srm{0}".format(i+1))
            else: 
                labelList.append("sr{0}".format(i+1))
            i = i + 1 
      elif(key=="NONE"): 
         print("No runs to use!")

   # extract the MIDAS run we want
   # first collect the midas runs into a separate list  
   midasRunList  = []
   nmrAnaRunList = []
   N = len(runList)
   for i in xrange(0,N): 
      theLabel = "srm{0}".format(i+1) 
      if(labelList[i]==theLabel): 
         midasRunList.append(runList[i])
      else: 
         nmrAnaRunList.append(runList[i]) 
   # build final runlist 
   fRunList   = []
   fLabelList = []
   # choose midas run that corresponds to the axis we're using
   NM = len(midasRunList) 
   if(NM==1): 
      # if we have a single MIDAS run, use that 
      fRunList.append(midasRunList[0])
      fLabelList.append("midas-run")  
   else: 
      # otherwise, we choose the index that corresponds to the axis we want  
      fRunList.append(midasRunList[axis])
      fLabelList.append("midas-run")  
   # add all NMR-ANA subruns -- these are determined from the key ppx, ppy, ppz
   # and hence should match the MIDAS run chosen  
   for sr in nmrAnaRunList: 
      fRunList.append(sr) 
   # create labels for the final run list   
   NRL = len(fRunList)
   for i in xrange(1,NRL): 
      fLabelList.append( "sr{0}".format(i+1) ) 
 
   # save to output json object  
   outData = {} 
   # config file data 
   outData['type']                   = confData['type'] 
   outData['blinding']               = confData['blinding'] 
   outData['run-period']             = confData['run-period'] 
   outData['prod-tag']               = confData['prod-tag'] 
   outData['nmr-ana-tag']            = confData['nmr-ana-tag']
   outData['fxpr-set']               = confData['fxpr-set']
   outData['free-proton-cor']        = confData['free-proton-cor'] 
   outData['load-trly-swap-times']   = confData['load-trly-swap-times'] 
   outData['use-aba-time-weight']    = confData['use-aba-time-weight']
   outData['use-trly-temp-cor']      = confData['use-trly-temp-cor'] 
   outData['use-osc-cor']            = confData['use-osc-cor'] 
   outData['cut-file']               = confData['cut-file'] 
   outData['num-events-to-avg']      = confData['num-events-to-avg']  
   outData['num-events-time-window'] = confData['num-events-time-window']  
   outData['fxpr-remove-drift']      = confData['fxpr-remove-drift']  
   outData['use-misalign-cor']       = confData['use-misalign-cor']  
    
   # probe-specific data  
   outData['date']                 = data['date']
   outData['trly-probe']           = data['trly-probe'] 
   outData['trly_azi-angle']       = data['trly_azi-angle'] 
   outData['dBz-current']          = data['dBz-current'] 

   if(axis==0): 
      outData['load-pp-scc-times']   = data['dB-pp_x']['load-times']
      outData['load-trly-scc-times'] = data['dB-trly_x']['load-times']
   elif(axis==1): 
      outData['load-pp-scc-times']   = data['dB-pp_y']['load-times']
      outData['load-trly-scc-times'] = data['dB-trly_y']['load-times']
   elif(axis==2):  
      outData['load-pp-scc-times']   = data['dB-pp_z']['load-times']
      outData['load-trly-scc-times'] = data['dB-trly_z']['load-times']

   if(fitData): 
      outData['fit'] = data[tag]['fit']  
   else: 
      outData['fit'] = "NONE"

   outData['nruns']      = len(fRunList) 
   outData['run-list']   = fRunList   
   outData['run-label']  = fLabelList  
   outData['axis']       = axis 
   # print("{0}".format(outData) )
   outfile = open(outpath,'w')
   json.dump(outData,outfile)
   outfile.close()
#_______________________________________________________________________________
def writeConfigFileProd_params(confData,outpath): 
   # write a JSON file that has all the analysis parameters relevant for the analysis 
   # fill output json object/dictionary   
   outData = {} 
   outData["run-period"]             = confData["run-period"] 
   outData["blinding"]               = confData["blinding"] 
   outData["prod-tag"]               = confData["prod-tag"] 
   outData["nmr-ana-tag"]            = confData["nmr-ana-tag"] 
   outData["cut-file"]               = confData["cut-file"] 
   outData["free-proton-cor"]        = confData["free-proton-cor"] 
   outData["load-trly-swap-times"]   = confData["load-trly-swap-times"] 
   outData["use-aba-time-weight"]    = confData["use-aba-time-weight"] 
   outData["use-trly-temp-cor"]      = confData["use-trly-temp-cor"] 
   outData["use-osc-cor"]            = confData["use-osc-cor"] 
   outData["fxpr-remove-drift"]      = confData["fxpr-remove-drift"] 
   outData["use-misalign-cor"]       = confData["use-misalign-cor"] 
   outData["num-events-to-avg"]      = confData["num-events-to-avg"] 
   outData["num-events-time-window"] = confData["num-events-time-window"]   
   outData["fxpr-set"]               = confData["fxpr-set"]   
   # write to file  
   if os.path.isfile(outpath): 
      print("[writeConfigFileProd_params]: File {0} exists, deleting first".format(outpath) )
      os.remove(outpath) 
   with open(outpath,'a') as outfile: 
       json.dump(outData,outfile,indent=2)
   print("[writeConfigFileProd_params]: File {0} written".format(outpath)) 
#_______________________________________________________________________________
def writeConfigFileProd_imposedGrad(confData,tag,keyList,axis,fitData,outpath): 
   # write a JSON file for the ROOT script defined by tag 
   # build the output json object 
   outData = {}
   # config file data 
   outData['type']                   = confData['type'] 
   outData['blinding']               = confData['blinding'] 
   outData['run-period']             = confData['run-period']
   outData['use-osc-cor']            = confData['use-osc-cor'] 
   # outData['prod-tag']               = confData['prod-tag'] 
   # outData['nmr-ana-tag']            = confData['nmr-ana-tag']
   # outData['fxpr-set']               = confData['fxpr-set']
   # outData['free-proton-cor']        = confData['free-proton-cor'] 
   # outData['load-trly-swap-times']   = confData['load-trly-swap-times'] 
   # outData['use-aba-time-weight']    = confData['use-aba-time-weight']
   # outData['use-trly-temp-cor']      = confData['use-trly-temp-cor'] 
   # outData['cut-file']               = confData['cut-file'] 
   # outData['num-events-to-avg']      = confData['num-events-to-avg']  
   # outData['num-events-time-window'] = confData['num-events-time-window']  
   # outData['fxpr-remove-drift']      = confData['fxpr-remove-drift']  
   # outData['use-misalign-cor']       = confData['use-misalign-cor']  

   # print("{0}".format(outData) )
   outfile = open(outpath,'w')
   json.dump(outData,outfile)
   outfile.close()
#_______________________________________________________________________________
def writeConfigFileProd(data,confData,tag,keyList,axis,fitData,outpath): 
   # write a JSON file for the ROOT script defined by tag 
   # get the run list 
   runList   = []
   labelList = []  
   subList   = []
   i=0
   for key in keyList: 
      if(key=="ppx" or key=="ppy" or key=="ppz" or key=="midas-runs"):
         # for a list of PP NMR-DAQ runs or MIDAS runs  
         subList = data[tag][key] # this is a list!  
         for subRun in subList: 
            runList.append(subRun) 
            labelList.append("sr{0}".format(i+1))
            i = i + 1 
      elif(key=="NONE"): 
         print("No runs to use!") 
      else: 
         # single MIDAS runs 
         runList.append(data[tag][key])   
         labelList.append(key) 
   # now to build the output json object 
   outData = {}
   # config file data 
   outData['type']                   = confData['type'] 
   outData['blinding']               = confData['blinding'] 
   outData['run-period']             = confData['run-period'] 
   outData['prod-tag']               = confData['prod-tag'] 
   outData['nmr-ana-tag']            = confData['nmr-ana-tag']
   outData['fxpr-set']               = confData['fxpr-set']
   outData['free-proton-cor']        = confData['free-proton-cor'] 
   outData['load-trly-swap-times']   = confData['load-trly-swap-times'] 
   outData['use-aba-time-weight']    = confData['use-aba-time-weight']
   outData['use-trly-temp-cor']      = confData['use-trly-temp-cor'] 
   outData['use-osc-cor']            = confData['use-osc-cor'] 
   outData['cut-file']               = confData['cut-file'] 
   outData['num-events-to-avg']      = confData['num-events-to-avg']  
   outData['num-events-time-window'] = confData['num-events-time-window']  
   outData['fxpr-remove-drift']      = confData['fxpr-remove-drift']  
   outData['use-misalign-cor']       = confData['use-misalign-cor']  
    
   # probe-specific data  
   outData['date']                 = data['date']
   outData['trly-probe']           = data['trly-probe']
   outData['trly_azi-angle']       = data['trly_azi-angle'] 
   outData['dBz-current']          = data['dBz-current'] 

   if(axis==0): 
      outData['load-pp-scc-times']   = data['dB-pp_x']['load-times']
      outData['load-trly-scc-times'] = data['dB-trly_x']['load-times']
   elif(axis==1): 
      outData['load-pp-scc-times']   = data['dB-pp_y']['load-times']
      outData['load-trly-scc-times'] = data['dB-trly_y']['load-times']
   elif(axis==2):  
      outData['load-pp-scc-times']   = data['dB-pp_z']['load-times']
      outData['load-trly-scc-times'] = data['dB-trly_z']['load-times']

   if(fitData): 
      outData['fit'] = data[tag]['fit']  
   else: 
      outData['fit'] = "NONE"

   outData['nruns']      = len(runList) 
   outData['run-list']   = runList   
   outData['run-label']  = labelList  
   outData['axis']       = axis 
   # print("{0}".format(outData) )
   outfile = open(outpath,'w')
   json.dump(outData,outfile)
   outfile.close()
#_______________________________________________________________________________
def writeConfigFileSingleRun(data,tag,label,runNumber,isFullAnalysis,isFinalLocation,axis,fitData,outpath): 
   # write a JSON file for the ROOT script defined by tag 
   # get the run list 
   runList   = []
   labelList = []  
   runList.append(runNumber)   
   labelList.append(label)  
   outData = {} 
   outData['date']       = data['date']
   outData['type']       = data['type'] 
   outData['blinding']   = data['blinding'] 
   outData['trly-probe'] = data['trly-probe'] 
   outData['p2p-fit']    = data['p2p-fit'] 
   if(fitData): 
      outData['fit'] = data[tag]['fit']  
   else: 
      outData['fit'] = "NONE"
   outData['full-ana']   = isFullAnalysis 
   outData['final-loc']  = isFinalLocation 
   outData['nruns']      = len(runList) 
   outData['run-list']   = runList   
   outData['run-label']  = labelList  
   outData['fxpr-set']   = data['fxpr-set']  
   outData['axis']       = axis 
   # print("{0}".format(outData) )
   outfile = open(outpath,'w')
   json.dump(outData,outfile)
   outfile.close()
#_______________________________________________________________________________
def writeConfigFileSimple(inData,tag,outpath): 
   # write a JSON file for a ROOT script 
   # obtain a run list defined by tag, write to outpath  
   # get the run list 
   runList   = inData[tag]
   NRUN = len(runList)
   labelList = [] 
   for i in xrange(0,NRUN): 
      labelList.append("{0}".format(i+1)) 
   # build output object  
   outData = {} 
   outData['date']       = inData['date']
   outData['blinding']   = inData['blinding']
   outData['type']       = inData['type'] 
   outData['trly-probe'] = inData['trly-probe'] 
   outData['nruns']      = len(runList) 
   outData['run-list']   = runList   
   outData['run-label']  = labelList 
   outData['fxpr-set']   = inData['fxpr-set']  
   # print("{0}".format(outData) )
   outfile = open(outpath,'w')
   json.dump(outData,outfile)
   outfile.close() 
#_______________________________________________________________________________
def getPrimaryValue(data,tag,key1): 
   return data[tag][key1] 
#_______________________________________________________________________________
def getSecondaryValue(data,tag,key1,key2): 
   return data[tag][key1][key2]
#_______________________________________________________________________________
def getShimRun(data,tag): 
   return data[tag]["shim"]
#_______________________________________________________________________________
def writeToFileShimGrad_pp(prefix,date,key,data): 
   fileName = "{0}/pp-shimmed-scan_{1}_{2}.csv".format(prefix,key,date)
   ppRun    = data["pp-shimmed-grad"]["pp"]
   trlyRun  = data["pp-shimmed-grad"]["trly"]
   runList  = data["pp-shimmed-grad"][key] 
   linePP   = "{0},pp,0".format(ppRun) 
   lineTRLY = "{0},trly,0".format(trlyRun) 
   f        = open(fileName,"w+")
   f.write(linePP+"\n")
   i = 0  
   for run in runList:
      i    = i + 1 
      line = "{0},sr{1},0".format(run,i) 
      f.write(line+"\n")
   f.write(lineTRLY+"\n") 
   f.close()
#_______________________________________________________________________________
def writeToFileShimGrad_pps(prefix,date,data): 
   # write out all the files we need for shimmed field scans with the PP 
   # first figure out which axes we have data for
   scanAxes = getPrimaryValue(data,"pp-shimmed-grad","axes")
   if( scanAxes == "x" ):
      writeToFileShimGrad_pp(prefix,date,"ppx",data)
   elif( scanAxes == "y" ):
      writeToFileShimGrad_pp(prefix,date,"ppy",data)
   elif( scanAxes == "z" ):
      writeToFileShimGrad_pp(prefix,date,"ppz",data)
   elif( scanAxes == "xy" ):
      writeToFileShimGrad_pp(prefix,date,"ppx",data)
      writeToFileShimGrad_pp(prefix,date,"ppy",data)
   elif( scanAxes == "xz" ):
      writeToFileShimGrad_pp(prefix,date,"ppx",data)
      writeToFileShimGrad_pp(prefix,date,"ppz",data)
   elif( scanAxes == "yz" ):
      writeToFileShimGrad_pp(prefix,date,"ppy",data)
      writeToFileShimGrad_pp(prefix,date,"ppz",data)
   elif( scanAxes == "xyz" ):
      writeToFileShimGrad_pp(prefix,date,"ppy",data)
      writeToFileShimGrad_pp(prefix,date,"ppy",data)
      writeToFileShimGrad_pp(prefix,date,"ppz",data)
#_______________________________________________________________________________
def writeToFileGrad_pps(prefix,date,tag,gradType,data): 
   fileName = "{0}/pp-scan_{1}-grad_{2}.csv".format(prefix,gradType,date)
   # print("filename = {0}, tag = {1}, date = {2}, gradType = {3}".format(fileName,tag,date,gradType) )
   runList  = data[tag]
   f        = open(fileName,"w+")
   line0    = "{0},{1},0".format(runList["bare"],"bare") 
   line1    = "{0},{1},0.1".format(runList[gradType],gradType) 
   f.write(line0+"\n")
   f.write(line1+"\n")
   f.close()
#_______________________________________________________________________________
def writeToFileGrad(prefix,date,tag,gradType,data,writeLastLine): 
   fileName = "{0}/{1}_{2}.csv".format(prefix,tag,date)
   # print("filename = {0}, tag = {1}, date = {2}, gradType = {3}".format(fileName,tag,date,gradType) )
   runList  = data[tag]
   f        = open(fileName,"w+")
   line0    = "{0},{1},0".format(runList["bare"],"bare") 
   line1    = "{0},{1},0.1".format(runList[gradType],gradType) 
   f.write(line0+"\n")
   f.write(line1+"\n")
   if(writeLastLine): 
      line2 = "{0},{1},0".format(runList["bare-2"],"bare-2") 
      f.write(line2+"\n")
   f.close() 
#_______________________________________________________________________________
def writeToFileTransGrad(prefix,date,data): 
   fileName = "{0}/trans-grad_{1}.csv".format(prefix,date)
   runList  = data["trans-grad"]
   f        = open(fileName,"w+")
   line0    = "{0},{1},0".format(runList["bare"]       ,"bare"     ) 
   line1    = "{0},{1},0.1".format(runList["norm-quad"],"norm-quad") 
   line2    = "{0},{1},0.1".format(runList["skew-quad"],"skew-quad") 
   line3    = "{0},{1},0".format(runList["bare-2"]     ,"bare-2"   ) 
   f.write(line0+"\n")
   f.write(line1+"\n")
   f.write(line2+"\n")
   f.write(line3+"\n")
   f.close() 
#_______________________________________________________________________________
def writeToFileShim(prefix,date,tag,data): 
   fileName = "{0}/{1}_{2}.csv".format(prefix,tag,date)
   runList  = data[tag]
   f        = open(fileName,"w+")
   line     = "{0},{1},0".format(runList["shim"],"shim") 
   f.write(line+"\n")
   f.close()

#_______________________________________________________________________________
def writeToFileShimGrad(prefix,date,tag,data): 
   fileName = "{0}/trly-shimmed_{1}_{2}.csv".format(prefix,tag,date)
   runList  = data["shimmed-grad"]
   f        = open(fileName,"w+")
   line     = "{0},{1},0".format(runList[tag],"shim") 
   f.write(line+"\n")
   f.close()
#_______________________________________________________________________________
def writeToFileTRLYMode(prefix,tag,data): 
   fileName = "{0}/trly-mode.csv".format(prefix) 
   runList  = data[tag] 
   f        = open(fileName,"a")
   for entry in runList: 
      line = "{0},{1}".format(entry,tag) 
      f.write(line+"\n") 
   f.close()  
#_______________________________________________________________________________
def prepareRunLists(theDate,data,ppScan): 
   prefix          = "./input/runlists/{0}".format(theDate) 
   gradType        = ["norm-quad"  ,"skew-quad"   ,"azi"]
   fileListPPs     = ["pp-scan-rad","pp-scan-vert","pp-scan-azi"]
   fileListTRLY    = ["dB-trly_rad-grad","dB-trly_vert-grad","dB-trly_azi-grad"]
   fileListTRLY_fl = ["dB-trly_final-location_rad-grad","dB-trly_final-location_vert-grad","dB-trly_final-location_azi-grad"]
   fileListPP      = ["dB-pp_rad-grad"  ,"dB-pp_vert-grad"  ,"dB-pp_azi-grad"  ]
   fileListPP_fl   = ["dB-pp_final-location_rad-grad"  ,"dB-pp_final-location_vert-grad","dB-pp_final-location_azi-grad"]
   # make a directory to store everything
   if not os.path.exists(prefix):
      # doesn't already exist, let's create it  
      os.makedirs(prefix) 
   else:
      # directory exists, remove it then create 
      shutil.rmtree(prefix) 
      os.makedirs(prefix)
   if(ppScan): 
      # dB and imposed gradients for the PP
      i=0
      for fn in fileListPPs: 
         writeToFileGrad_pps(prefix,theDate,fn,gradType[i],data) 
         i = i + 1
      # shimmed field scans
      writeToFileShimGrad_pps(prefix,theDate,data) 
   else: 
      # dB for the PP
      i=0
      for fn in fileListPP: 
         writeToFileGrad(prefix,theDate,fn,gradType[i],data,True) 
         i = i + 1
      # transverse gradients (imposed) 
      writeToFileTransGrad(prefix,theDate,data)  
      # azimuthal gradients (imposed)
      writeToFileGrad(prefix,theDate,"azi-grad","azi",data,True) 
      # gradients in shimmed field 
      writeToFileShimGrad(prefix,theDate,"trans",data)  
      writeToFileShimGrad(prefix,theDate,"azi"  ,data)  
   # dB for the TRLY 
   i=0
   for fn in fileListTRLY: 
      writeToFileGrad(prefix,theDate,fn,gradType[i],data,True) 
      i = i + 1
   # dB at final location for PP 
   i=0
   for fn in fileListPP_fl: 
      writeToFileGrad(prefix,theDate,fn,gradType[i],data,True) 
      i = i + 1
   # dB at final location for TRLY
   i=0
   for fn in fileListTRLY_fl: 
      writeToFileGrad(prefix,theDate,fn,gradType[i],data,True) 
      i = i + 1
   # shimmed data 
   writeToFileShim(prefix,theDate,"pp-shimmed"  ,data)  
   writeToFileShim(prefix,theDate,"trly-shimmed",data) 
   # trolley mode data 
   writeToFileTRLYMode(prefix,"interactive",data)  
   writeToFileTRLYMode(prefix,"continuous" ,data)  

def prepareRunListsDeltaBOnly(theDate,data,ppScan): 
   prefix          = "./input/runlists/{0}".format(theDate) 
   gradType        = ["norm-quad"  ,"skew-quad"   ,"azi"]
   fileListTRLY    = ["dB-trly_rad-grad","dB-trly_vert-grad","dB-trly_azi-grad"]
   fileListPP      = ["dB-pp_rad-grad"  ,"dB-pp_vert-grad"  ,"dB-pp_azi-grad"  ]
   # make a directory to store everything
   if not os.path.exists(prefix):
      # doesn't already exist, let's create it  
      os.makedirs(prefix) 
   else:
      # directory exists, remove it then create 
      shutil.rmtree(prefix) 
      os.makedirs(prefix) 
   i=0
   if(ppScan): 
      # dB for the PP
      for fn in gradType: 
         writeToFileGrad_pps(prefix,theDate,"pp-scan",gradType[i],data) 
         i = i + 1
   else: 
      # dB for the PP
      for fn in fileListPP: 
         writeToFileGrad(prefix,theDate,fn,gradType[i],data,True) 
         i = i + 1
      # transverse gradients (imposed) 
      writeToFileTransGrad(prefix,theDate,data)  
      # azimuthal gradients (imposed)
      writeToFileGrad(prefix,theDate,"azi-grad","azi",data)
   # dB for the TRLY 
   i=0
   for fn in fileListTRLY: 
      writeToFileGrad(prefix,theDate,fn,gradType[i],data,True) 
      i = i + 1
   # trolley mode data 
   writeToFileTRLYMode(prefix,"interactive",data)  
   writeToFileTRLYMode(prefix,"continuous" ,data)  
#_______________________________________________________________________________
def deleteDir(path): 
   shutil.rmtree(path) 
#_______________________________________________________________________________
def createDir(path): 
   os.makedirs(path)
#_______________________________________________________________________________
def setupPaths(deleteDir,today,isBlind,blindLabel): 
   # set up directories and paths to inputs and outputs 
   if(isBlind):
      outDir = "./output/blinded/" + blindLabel + "/" + today
   else:
      outDir = "./output/unblinded/" + today
   
   if deleteDir:
      if os.path.exists(outDir):
         print( "Removing existing output from {0}".format(outDir) )
         deleteDir(outDir)
  
   # setup output directories 
   if(isBlind):
      outDir  = "./output/blinded/" + blindLabel + "/" + today
      plotDir = "./plots/blinded/" + blindLabel + "/" + today
   else:
      outDir  = "./output/unblinded/" + today
      plotDir = "./plots/unblinded/"  + today
   
   if not os.path.exists(outDir):
       os.makedirs(outDir)
       print("[run_analysis]: Created directory: {0}".format(outDir) )
   
   if not os.path.exists(plotDir):
       os.makedirs(plotDir)
       print("[run_analysis]: Created directory: {0}".format(plotDir) )

   return outDir,plotDir

