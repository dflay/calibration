# Less extensive comparisons of data 
# focus on calibration coeffs  

import csv
import json
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as md
import numpy as np
import math

#_______________________________________________________________________________
def getDiff_2df(axis,axisErr,newAxis,newAxisErr,df1,df2,df3):
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
def getDiff(a1,a1Err,a2,a2Err,newAxis,newAxisErr,df,df_out):

   # get lists of data
   x  = df[a1].tolist()
   ex = df[a1Err].tolist()
   y  = df[a2].tolist()
   ey = df[a2Err].tolist()

   z  = [] 
   ez = []

   # compute differences 
   getDiff_lists(x,ex,y,ey,z,ez)

   # add to data frame 
   colNum = len(df_out.columns) 
   df_out.insert(colNum,newAxis   ,z) 
   df_out.insert(colNum,newAxisErr,ez) 

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

plt.tight_layout()

# create file paths

# create a pandas dataframe, reading in the csv file  
csv_path = "./output/blinded/flay/04-07-20/run-1/calibData_04-07-20.csv"
print("Reading data from: {0}".format(csv_path)) 
data_df1 = pd.read_csv(csv_path,index_col=False) # index_col = False when you don't have an index column

csv_path = "./input/ran-hong/run-1_04-11-20.csv"
print("Reading data from: {0}".format(csv_path)) 
data_df2 = pd.read_csv(csv_path,index_col=False) # index_col = False when you don't have an index column

csv_path = "./input/bingzhi-li/run-1_04-06-20.csv"
print("Reading data from: {0}".format(csv_path)) 
data_df3 = pd.read_csv(csv_path,index_col=False) # index_col = False when you don't have an index column

# marker parameters  
color  = ["blue","red","#20B010"]
mStyle = ["s","o","v"]
mSize  = 80 

# calculate differences 
probeList = data_df1['Probe'].tolist() # a list of the probe numbers to add to the new data frame 
data_diff = pd.DataFrame() # new dataframe of differences 
data_diff.insert(0,"Probe",probeList)

# a comparison dataframe 
data_comp = pd.DataFrame() # new dataframe of differences 
data_comp.insert(0,"Probe",probeList)

cc1  = data_df1['calibCoeff'].tolist()
cc1e = data_df1['calibCoeffErr'].tolist()
data_comp.insert(1,'calibCoeff_1',cc1)

ccc1  = data_df1['calibCoeff_cor'].tolist()
ccc1e = data_df1['calibCoeffErr_cor'].tolist()
data_comp.insert(1,'calibCoeff_cor_1',ccc1)
 
cc2  = data_df2['calibCoeff'].tolist()
cc2e = data_df2['calibCoeffErr'].tolist()
data_comp.insert(1,'calibCoeff_2',cc1)

ccc2  = data_df2['calibCoeff_cor'].tolist()
ccc2e = data_df2['calibCoeffErr_cor'].tolist()
data_comp.insert(1,'calibCoeff_cor_2',ccc1) 

cc3  = data_df3['calibCoeff'].tolist()
cc3e = data_df3['calibCoeffErr'].tolist()
data_comp.insert(1,'calibCoeff_3',cc1)

ccc3  = data_df3['calibCoeff_cor'].tolist()
ccc3e = data_df3['calibCoeffErr_cor'].tolist()
data_comp.insert(1,'calibCoeff_cor_3',ccc1)

# estimate of uncertainty on DF vs RH 
NP  = len(cc1) 
ee2  = [] 
cee2 = [] 
ee3  = [] 
cee3 = [] 
for i in xrange(0,NP):
   arg = math.sqrt( cc1e[i]*cc1e[i] + cc2e[i]*cc2e[i] ) 
   ee2.append(arg) 
   arg = math.sqrt( ccc1e[i]*ccc1e[i] + ccc2e[i]*ccc2e[i] ) 
   cee2.append(arg) 
   arg = math.sqrt( cc1e[i]*cc1e[i] + cc3e[i]*cc3e[i] ) 
   ee3.append(arg) 
   arg = math.sqrt( ccc1e[i]*ccc1e[i] + ccc3e[i]*ccc3e[i] ) 
   cee3.append(arg) 

data_comp.insert(1,'calibCoeff_2Err',ee2)  
data_comp.insert(1,'calibCoeff_cor_2Err',cee2)  
data_comp.insert(1,'calibCoeff_3Err',ee3)  
data_comp.insert(1,'calibCoeff_cor_3Err',cee3)  

# take differences 
getDiff_2df("calibCoeff"        ,"calibCoeffErr"        ,"rh_cc_diff"  ,"rh_cce_diff"    ,data_df1,data_df2,data_diff)
getDiff_2df("calibCoeff_cor"    ,"calibCoeffErr_cor"    ,"rh_ccc_diff" ,"rh_ccce_diff"   ,data_df1,data_df2,data_diff)
getDiff_2df("calibCoeff"        ,"calibCoeffErr"        ,"bl_cc_diff"  ,"bl_cce_diff"    ,data_df1,data_df3,data_diff)
getDiff_2df("calibCoeff_cor"    ,"calibCoeffErr_cor"    ,"bl_ccc_diff" ,"bl_ccce_diff"   ,data_df1,data_df3,data_diff)

# getDiff_2df("calibCoeffFree"    ,"calibCoeffFreeErr"    ,"fcc_diff"  ,"fcce_diff"  ,data_df1,data_df2,data_diff)
# getDiff_2df("calibCoeffFree_cor","calibCoeffFreeErr_cor","fccc_diff" ,"fccce_diff" ,data_df1,data_df2,data_diff)

print data_diff

# plot data

# general setup  
NCOL = 1
NROW = 2

tickSize      = 16
xAxisFontSize = 16
yAxisFontSize = 16

# calib coeffs (NO FREE PROTON)
fig = plt.figure(1) 
plt.subplot(NROW,NCOL,1)

currentAxis = plt.gca() # grab current axis 
axis    = "calibCoeff_cor"
axisErr = "calibCoeffErr_cor"
data_df1.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[0], s=mSize, color=color[0], ax=currentAxis)
data_df2.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
data_df3.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.legend(["DF","RH","BL"])
currentAxis.set_xlabel("Probe"           , fontsize=xAxisFontSize) 
currentAxis.set_ylabel("Calib Coeff [Cor] (Hz)", fontsize=yAxisFontSize) 
currentAxis.tick_params(labelsize=tickSize)

plt.subplot(NROW,NCOL,2)
currentAxis = plt.gca() # grab current axis 
axis    = "rh_ccc_diff"
axisErr = "rh_ccce_diff"
data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
axis    = "bl_ccc_diff"
axisErr = "bl_ccce_diff"
data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.set_xlabel("Probe"                , fontsize=xAxisFontSize) 
currentAxis.set_ylabel("Calib Coeff [Cor] Diff (Hz)", fontsize=yAxisFontSize) 
currentAxis.tick_params(labelsize=tickSize)

for ax in fig.get_axes():
    ax.label_outer()

# plot calib coeffs against one another 
fig = plt.figure(2) 

label_1 = "Calibration Coefficients (DF)"
label_2 = "Calibration Coefficients (RH or BL)"

currentAxis = plt.gca() # grab current axis 
data_comp.plot(kind="scatter", x="calibCoeff_cor_1", y="calibCoeff_cor_2", yerr="calibCoeff_cor_2Err", marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
data_comp.plot(kind="scatter", x="calibCoeff_cor_1", y="calibCoeff_cor_3", yerr="calibCoeff_cor_3Err", marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.set_xlabel(label_1, fontsize = xAxisFontSize) 
currentAxis.set_ylabel(label_2, fontsize = yAxisFontSize)
currentAxis.tick_params(labelsize=tickSize)

for ax in fig.get_axes():
    ax.label_outer()

# # calib coeffs (WITH FREE PROTON)
# fig = plt.figure(2) 
# plt.subplot(NROW,NCOL,1)
# 
# currentAxis = plt.gca() # grab current axis 
# axis    = "calibCoeffFree"
# axisErr = "calibCoeffFreeErr"
# data_df.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[0], s=mSize, color=color[0], ax=currentAxis)
# axis    = "calibCoeffFree_cor"
# axisErr = "calibCoeffFreeErr_cor"
# data_df.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
# currentAxis.legend(["Without Misalignment","With Misalignment"])
# currentAxis.set_xlabel("Probe") 
# currentAxis.set_ylabel("Calib Coeff (Hz)") 
# 
# plt.subplot(NROW,NCOL,2)
# currentAxis = plt.gca() # grab current axis 
# axis    = "fcc_diff"
# axisErr = "fcce_diff"
# data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
# currentAxis.set_xlabel("Probe") 
# currentAxis.set_ylabel("Calib Coeff Diff (Hz)") 
# 
# for ax in fig.get_axes():
#     ax.label_outer()


# show all plots
plt.show()

