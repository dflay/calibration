#! /usr/bin/python 

# plot PP-DAQ data sets 

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
def convTemp(tl,rl,res):
   # first args are from the tempearture table, last is the resistance we have 
   index = np.searchsorted(rl,res,side='right')
   i = index[0]
   # linear interpolation 
   b         = (res-rl[i-1])/(rl[i]-rl[i-1])
   temp_degC = tl[i-1] + b*(tl[i]-tl[i-1])
   return temp_degC
#_______________________________________________________________________________

# gather command-line input 
NA = len(sys.argv)
if(NA!=3):
   print("Usage:   python {0} run-str sub-dir".format(sys.argv[0]) )
   print("Example: python {0} 2000-2030 my-dir".format(sys.argv[0]) )
   sys.exit(0)

# get list of runs 
runStr = sys.argv[1]
myParser = StringParser(runStr)
myParser.GenerateList()
runList = myParser.fList

# directory where runs are stored
dirName = sys.argv[2]

headerNames = ["trace","channel","time_ns","zc","nc","ampl","rms_noise","t2_time","temp","freq_LO","freq_pi2",
               "freq_mid","freq_lin","freq_lsq","freq_mid_ph","freq_lin_ph","freq_lsq_ph"]

prefix = "./input/NMR-ANA/{0}".format(dirName)

csv_path = "{0}/run-{1:05d}/results.csv".format(prefix,int(runList[0]) )
# df = pd.read_csv(csv_path,header="None") 
df = pd.read_csv(csv_path,names=headerNames) 
# run column
rList = []
for j in xrange(0,10): 
   rList.append(runList[0]) 

NR = len(runList)
for i in xrange(1,NR):
   csv_path = "{0}/run-{1:05d}/results.csv".format(prefix,int(runList[i]) )
   # print("plot_nmrAna: Reading in {0}".format(csv_path) )
   for j in xrange(0,10): 
      rList.append(runList[i]) 
   df_new = pd.read_csv(csv_path,names=headerNames)
   df = df.append(df_new)

# print rList
NRR = len(rList) 
ND  = len(df.index)
if(ND!=NRR):
   print("[plot_nmrAna]: ERROR! runList size: {0}, rList size: {1}, dataFrame size: {2}".format(NR,NRR,ND) )
   sys.exit(0)

# get PT1000 conversion table 
df_temp_table = pd.read_csv("/home/dflay/Workplace/dflay/temp-tables/PT1000.csv")  

# add run list as a column
df.insert(0,"nmrAna_run",rList) 

# convert time to sec 
df['time_ns'] = df['time_ns']/1E+9 
df.rename(columns = {'time_ns':'time'}, inplace=True)

# convert temperature from Ohms to deg C
# get list of temperatures 
theTMP_ohms = df['temp'].tolist() # in Ohms
TMP = []
NT = len(theTMP_ohms)
# calculate value in deg C
for i in xrange(0,NT):  
   tmp = convTemp(df_temp_table['Temperature'],df_temp_table['Resistance (Ohms)'],theTMP_ohms[i])
   # print("{0:.3f}, {1:.3f}".format(theTMP_ohms[i],tmp) )
   TMP.append(tmp)
# add a new column with temp in deg C 
df.insert(0,"temp_degC",TMP) 

# convert times to time stamps
df['time'] = df['time'].astype("datetime64[s]") 

# print data
print df

# create the figure
NROW = 3
NCOL = 1 
fig, ax = plt.subplots(nrows=NROW, ncols=NCOL)
ax[0].set_title("PP Data")

yAxisName  = ["freq_lsq_ph","temp_degC","ampl"]
yAxisLabel = ["Frequency (Hz)","Temperature (deg C)","Amplitude (V)"]

# create the plot
# arguments of subplot are: row-num col-num sub-canvas-num

xAxis = 'time'
NP = len(yAxisName)
LST = []
for i in xrange(0,NP):
   canvasIndex = 311 + i
   plt.subplot(canvasIndex)
   # determine y axis 
   yAxis = yAxisName[i]
   # add data to plot 
   plt.plot(df[xAxis],df[yAxis])
   # get limits 
   LST = df[yAxis].tolist()
   mean  = np.mean(LST)  
   stdev = np.std(LST)
   ylo = mean - 3.*stdev  
   yhi = mean + 3.*stdev 
   print("axis {0}, mean = {1:.3f}, stdev = {2:.3f}, ylo = {3:.3f}, yhi = {4:.3f}".format(i,mean,stdev,ylo,yhi) )
   if(yAxisName[i]=='temp_degC'):
      ylo = mean - 10.*stdev 
      yhi = mean + 10.*stdev 
   elif(yAxisName[i]=='ampl'):
      yli = mean - 0.2 
      yhi = mean + 0.2 
   # set graph limits
   plt.ylim(bottom=ylo,top=yhi)
   # set up labels -- DOES NOT WORK
   plt.xlabel('')
   plt.ylabel(yAxisLabel[i])
   # ax[i].label_outer() # for stacked plots, only label the outer ones
   # format the time axis -- DOES NOT WORK
   xfmt = md.DateFormatter('%Y-%m-%d %H:%M:%S')
   ax[i].xaxis.set_major_formatter(xfmt)
   # puts the date on an angle -- THIS WORKS 
   plt.gcf().autofmt_xdate()

plt.show()
 
