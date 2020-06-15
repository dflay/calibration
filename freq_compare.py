#! /usr/bin/python 

# Compare NMR-ANA output against simulation data 

import csv
import json
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as md
import numpy as np

from parse import StringParser

#_______________________________________________________________________________
def getDiff(axis,axisErr,newAxis,newAxisErr,df1,df2,df3):
   # take differences of specific columns from different datafames 
   # get lists of data
   x  = df1[axis].tolist()
   ex = df1[axisErr].tolist()
   y  = df2[axis].tolist()
   ey = df2[axisErr].tolist()

   z  = []
   ez = []

   # compute differences 
   getDiff_lists(x,ex,y,ey,z,ez)

   # add to data frame 
   colNum = len(df3.columns)
   df3.insert(colNum,newAxis   ,z)
   df3.insert(colNum,newAxisErr,ez)

   return
#_______________________________________________________________________________
def getList(probeKey,colKey,df):
   # grab rows that match a set of keys
   p = df['probe'].tolist() 
   x = df[colKey].tolist()
   y = [] # output list 
   N = len(p)
   for i in xrange(0,N):
      arg = int(p[i]) 
      if arg==probeKey:
        # found the matching probe number, fill output list
        y.append(x[i]) 
   return y 
#_______________________________________________________________________________
def getDiff_lists(x,y,z):
   # take differences of two lists
   N = len(x)
   for i in xrange(0,N):
      # protect against weird strings 
      xf  = float(x[i])
      yf  = float(y[i])
      # compute differences  
      arg = yf-xf
      # fill list 
      z.append(arg)
   return
#_______________________________________________________________________________
def getDiff_lists_werr(x,xe,y,ye,z,ze):
   # take differences of two lists
   N = len(x)
   for i in xrange(0,N):
      # protect against weird strings 
      xf  = float(x[i])
      xfe = float(xe[i])
      yf  = float(y[i])
      yfe = float(ye[i])
      # compute differences  
      arg = yf-xf
      arg_err = math.sqrt( xfe*xfe + yfe*yfe )
      # fill lists 
      z.append(arg)
      ze.append(arg_err)
   return
#_______________________________________________________________________________
def getStats(colName,df):
   # get mean and standard deviation of a column from dataframe df
   # first convert to float just to make sure
   x     = df[colName].tolist()
   xf    = []
   N = len(x)
   for entry in x:
      xf.append( float(entry) )
   mean  = np.mean(xf)
   stdev = np.std(xf)
   return mean,stdev
#_______________________________________________________________________________

# gather command-line input 
NA = len(sys.argv)
if(NA!=3):
   print("Usage:   python {0} run-str sub-dir".format(sys.argv[0]) )
   print("Example: python {0} 2000-2030 my-dir".format(sys.argv[0]) )
   sys.exit(0)

# get list of runs 
runStr   = sys.argv[1]
myParser = StringParser(runStr)
myParser.GenerateList()
runList = myParser.fList

# directory where runs are stored
dirName = sys.argv[2]

headerNames = ["trace","channel","time_ns","zc","nc","ampl","rms_noise","t2_time","temp_ohms","freq_LO","freq_pi2",
               "freq_mid","freq_lin","freq_lsq","freq_mid_ph","freq_lin_ph","freq_lsq_ph"]

prefix = "./input/NMR-ANA/{0}".format(dirName)

csv_path = "{0}/run-{1:05d}/results.csv".format(prefix,int(runList[0]) )
print("plot_nmrAna: Reading in {0}".format(csv_path) )
# df = pd.read_csv(csv_path,header="None") 
df = pd.read_csv(csv_path,names=headerNames) 
# run and probe columns
theProbe = 1
rList = [] # run 
pList = [] # probe list 
NS = len(df.index)
for j in xrange(0,NS): 
   rList.append(runList[0])
   pList.append(theProbe)  

NR = len(runList)
for i in xrange(1,NR):
   csv_path = "{0}/run-{1:05d}/results.csv".format(prefix,int(runList[i]) )
   print("plot_nmrAna: Reading in {0}".format(csv_path) )
   df_new = pd.read_csv(csv_path,names=headerNames)
   df     = df.append(df_new)
   NS = len(df_new.index)
   print("N shots = {0}".format(NS) )
   theProbe = theProbe + 1 
   for j in xrange(0,NS): 
      rList.append(runList[i])
      pList.append(theProbe)

# print rList
NRR = len(rList) 
ND  = len(df.index)
if(ND!=NRR):
   print("[plot_nmrAna]: ERROR! runList size: {0}, rList size: {1}, dataFrame size: {2}".format(NR,NRR,ND) )
   sys.exit(0)

# add run list as a column
df.insert(0,"nmrAna_run",rList) 
df.insert(0,"probe"     ,pList) 

# convert time to sec 
df['time_ns'] = df['time_ns']/1E+9 
df.rename(columns = {'time_ns':'time'}, inplace=True)

# convert times to time stamps
df['time'] = df['time'].astype("datetime64[s]") 

# print data
print df

# read in simulation truth data
csv_path = "./input/simulation/run-1_sim_info.csv" 
df_tru = pd.read_csv(csv_path)

print df_tru

# create a new data frame with differences 
data_diff = pd.DataFrame() # new dataframe of differences 
probeList = df['probe'].tolist()
data_diff.insert(0,'probe',pList) 

# difference of DF frequencies and simulated avg freq 
df_freq  = df['freq_lsq_ph'].tolist()
avg_freq = df_tru['avg_freq'].tolist()
rh_freq  = df_tru['rh_freq'].tolist()

# add DF results to truth table 
nc = len(df_tru.columns) 
df_tru.insert(nc,'df_freq',df_freq)

outpath = "sim_run-1_comp.csv"
df_tru.to_csv(outpath,sep=',',index=False,float_format='%.3f')

df_avg_diff = []  
getDiff_lists(avg_freq,df_freq,df_avg_diff) 
df_rh_diff = []  
getDiff_lists(rh_freq,df_freq,df_rh_diff) 

data_diff.insert(1,'df_avg_diff',df_avg_diff)
data_diff.insert(2,'df_rh_diff' ,df_rh_diff)

print data_diff

NROW = 9 
NCOL = 2 
fig, ax = plt.subplots(nrows=NROW, ncols=NCOL)
plt.subplots_adjust(left=None, bottom=None, right=None, top=0.9, wspace=None, hspace=None)

NTR = 7
trace = []
for i in xrange(0,NTR):
   trace.append(i+1) 

ri=0
ci=0

for i in xrange(0,17):
   x   = getList(i+1,'df_avg_diff',data_diff)
   print("Probe {0}".format(i+1))
   M = len(x)
   outStr = "{0:.3f}".format(x[0]) 
   for j in xrange(1,M):
     outStr = "{0}, {1:.3f}".format(outStr,x[j]) 
   print outStr
   # add to plot
   mu  = np.mean(x)  
   std = np.std(x)  
   ylo = mu - 3.*std  
   yhi = mu + 3.*std 
   if(i<9):
      ri = i
      ci = 0
   else: 
      ri = i-9
      ci = 1 
   ax[ri][ci].plot(trace,x)
   ax[ri][ci].set_title("Probe {0:02d}".format(i+1) )
   ax[ri][ci].set_ylim(ylo,yhi)
   # clean up  
   x[:] = [] 

# for ax in fig.get_axes():
#     ax.label_outer()

fig2, ax2 = plt.subplots(nrows=NROW, ncols=NCOL)
plt.subplots_adjust(left=None, bottom=None, right=None, top=0.9, wspace=None, hspace=None)

for i in xrange(0,17):
   x   = getList(i+1,'df_rh_diff',data_diff)
   print("Probe {0}".format(i+1))
   M = len(x)
   outStr = "{0:.3f}".format(x[0]) 
   for j in xrange(1,M):
     outStr = "{0}, {1:.3f}".format(outStr,x[j]) 
   print outStr
   # add to plot
   mu  = np.mean(x)  
   std = np.std(x)  
   ylo = mu - 3.*std  
   yhi = mu + 3.*std 
   if(i<9):
      ri = i
      ci = 0
   else: 
      ri = i-9
      ci = 1 
   ax2[ri][ci].plot(trace,x)
   ax2[ri][ci].set_title("Probe {0:02d}".format(i+1) )
   ax2[ri][ci].set_ylim(ylo,yhi)
   # clean up  
   x[:] = [] 

# for ax in fig.get_axes():
#     ax.label_outer()

fig2.tight_layout()

plt.show()
  
