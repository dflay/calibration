# read in calibration coeffs from other analyzers and create a json object 

import csv
import json
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as md
import numpy as np

#_______________________________________________________________________________
def getDiff(axis,axisErr,newAxis,newAxisErr,df1,df2,df3):

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
   colNum = len(df3.columns()) 
   df3.append(colNum,newAxis   ,z) 
   df3.append(colNum,newAxisErr,ez) 

   return
#_______________________________________________________________________________
def getDiff_lists(x,xe,y,ye,z,ze):

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


# create file paths
csv_path_rh = "./input/ran-hong/run-1_04-03-20.csv" 
csv_path_bl = "./input/bingzhi-li/run-1_04-05-20.csv" 
csv_path_df = "./output/blinded/flay/04-04-20/run-1/calibData_04-04-20.csv"
​
# create a pandas dataframe, reading in the csv file  
print("Reading data from: {0}".format(csv_path_rh)) 
data_rh = pd.read_csv(csv_path_rh,index_col=False) # index_col = False when you don't have an index column
print("Reading data from: {0}".format(csv_path_bl)) 
data_bl = pd.read_csv(csv_path_bl,index_col=False) # index_col = False when you don't have an index column
print("Reading data from: {0}".format(csv_path_df)) 
data_df = pd.read_csv(csv_path_df,index_col=False) # index_col = False when you don't have an index column
​
# colors for everyone 
color = ["blue","red","#20B010"]
​
# calculate differences 
probeList = data_df['Probe'].tolist() # a list of the probe numbers to add to the new data frame 
data_diff = pd.DataFrame() # new dataframe of differences 
data_diff.insert(0,"Probe",probeList)
​
# take differences relative to DF
getDiff("calibCoeff","calibCoeffErr","cc_rhdf","cc_rhdf_err",data_df,data_rh,data_diff)
getDiff("calibCoeff","calibCoeffErr","cc_bldf","cc_bldf_err",data_df,data_bl,data_diff)
​
print data_diff
​
# compare calibration coefficients (no corrections)
# probably need a subplot here -- bottom panel for differences 
currentAxis = plt.gca() # grab current axis 
data_df.plot(kind="scatter", x="Probe", y="calibCoeff", yerr="calibCoeffErr", marker="s", s=30, color=color[0], ax=currentAxis)
data_rh.plot(kind="scatter", x="Probe", y="calibCoeff", yerr="calibCoeffErr", marker="o", s=30, color=color[1], ax=currentAxis)
data_bl.plot(kind="scatter", x="Probe", y="calibCoeff", yerr="calibCoeffErr", marker="v", s=30, color=color[2], ax=currentAxis)
currentAxis.legend(["DF","RH","BL"])
plt.show()


