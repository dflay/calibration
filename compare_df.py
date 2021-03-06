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
import seaborn as sns
import math

#_______________________________________________________________________________
def printToScreen_listOfTypes(axis,df):
   # input axis = list of axes
   #       df   = existing dataframe 

   df_new = pd.DataFrame() 

   i=0
   for ax in axis:
      x = df[ax].tolist()
      df_new.insert(i,ax,x)      
      i = i + 1 
      x[:] = []  

   print df_new

   return
#_______________________________________________________________________________
def printToScreen(axis1,axis2,df):
   print("Data for {0} and {1}".format(axis1,axis2) ) 
   x = df[axis1].tolist()
   y = df[axis2].tolist()

   N = len(x) 
   for i in xrange(0,N): 
      print("probe {0:2d}: {1:.3f}, {2:.3f}".format(i+1,x[i],y[i]) )

   return
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

plotBingzhi = True 
plotRun2    = False 
removeBlind = False # using UNBLINDED files now! 

print("Loading data...")

pd.set_option("display.precision", 2) 

# create a pandas dataframe, reading in the csv file  
csv_path = "./output/unblinded/05-04-20/run-1/calibData_05-04-20.csv"
print("Reading data from: {0}".format(csv_path)) 
data_df1 = pd.read_csv(csv_path,index_col=False) # index_col = False when you don't have an index column
print data_df1.columns

csv_path = "./output/unblinded/05-04-20/run-2/calibData_05-04-20.csv"
print("Reading data from: {0}".format(csv_path)) 
data_df1_r2 = pd.read_csv(csv_path,index_col=False) # index_col = False when you don't have an index column
# print data_df1_r2[['deltaB_pp_z','deltaB_tr_z']]

csv_path = "./input/ran-hong/run-1_unblind_05-04-20.csv"
print("Reading data from: {0}".format(csv_path)) 
data_df2 = pd.read_csv(csv_path,index_col=False) # index_col = False when you don't have an index column
print data_df2.columns

csv_path = "./input/bingzhi-li/run-1_unblind_05-04-20.csv"
print("Reading data from: {0}".format(csv_path)) 
data_df3 = pd.read_csv(csv_path,index_col=False) # index_col = False when you don't have an index column
print data_df3.columns

print("--> Done.") 

blindDF = 2.54
blindRH = 8.60 
blindBL = 3.14

if(removeBlind):
   print("WARNING! REMOVING BLINDS!") 
   # data_df1['calibCoeff']     -= blindDF 
   # data_df1['calibCoeff_cor'] -= blindDF 
   # data_df2['calibCoeff']     -= blindRH 
   # data_df2['calibCoeff_cor'] -= blindRH 
   data_df3['calibCoeff']     -= blindBL 
   data_df3['calibCoeff_cor'] -= blindBL 

print("Making plots...") 

# marker parameters  
color  = ["blue","red","#20B010","magenta"]
mStyle = ["s","o","v","s"]
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

cc4  = data_df1_r2['calibCoeff'].tolist()
cc4e = data_df1_r2['calibCoeffErr'].tolist()
data_comp.insert(1,'calibCoeff_4',cc4)

ccc4  = data_df1_r2['calibCoeff_cor'].tolist()
ccc4e = data_df1_r2['calibCoeffErr_cor'].tolist()
data_comp.insert(1,'calibCoeff_cor_4',ccc4)

# estimate of uncertainty on DF vs RH 
NP  = len(cc1) 
ee2  = [] 
cee2 = [] 
ee3  = [] 
cee3 = [] 
ee4  = [] 
cee4 = [] 
for i in xrange(0,NP):
   arg = math.sqrt( cc1e[i]*cc1e[i] + cc2e[i]*cc2e[i] ) 
   ee2.append(arg) 
   arg = math.sqrt( ccc1e[i]*ccc1e[i] + ccc2e[i]*ccc2e[i] ) 
   cee2.append(arg) 
   arg = math.sqrt( cc1e[i]*cc1e[i] + cc3e[i]*cc3e[i] ) 
   ee3.append(arg) 
   arg = math.sqrt( ccc1e[i]*ccc1e[i] + ccc3e[i]*ccc3e[i] ) 
   cee3.append(arg) 
   arg = math.sqrt( cc1e[i]*cc1e[i] + cc4e[i]*cc4e[i] ) 
   ee4.append(arg) 
   arg = math.sqrt( ccc1e[i]*ccc1e[i] + ccc4e[i]*ccc4e[i] ) 
   cee4.append(arg) 

data_comp.insert(1,'calibCoeff_2Err',ee2)  
data_comp.insert(1,'calibCoeff_cor_2Err',cee2)  
data_comp.insert(1,'calibCoeff_3Err',ee3)  
data_comp.insert(1,'calibCoeff_cor_3Err',cee3)  
data_comp.insert(1,'calibCoeff_4Err',ee4)  
data_comp.insert(1,'calibCoeff_cor_4Err',cee4)  

# take differences 
getDiff_2df("calibCoeff"        ,"calibCoeffErr"        ,"rh_cc_diff"  ,"rh_cce_diff"    ,data_df1,data_df2,data_diff)
getDiff_2df("calibCoeff_cor"    ,"calibCoeffErr_cor"    ,"rh_ccc_diff" ,"rh_ccce_diff"   ,data_df1,data_df2,data_diff)
getDiff_2df("calibCoeff"        ,"calibCoeffErr"        ,"bl_cc_diff"  ,"bl_cce_diff"    ,data_df1,data_df3,data_diff)
getDiff_2df("calibCoeff_cor"    ,"calibCoeffErr_cor"    ,"bl_ccc_diff" ,"bl_ccce_diff"   ,data_df1,data_df3,data_diff)
getDiff_2df("calibCoeff"        ,"calibCoeffErr"        ,"r2_cc_diff"  ,"r2_cce_diff"    ,data_df1,data_df1_r2,data_diff)
getDiff_2df("calibCoeff_cor"    ,"calibCoeffErr_cor"    ,"r2_ccc_diff" ,"r2_ccce_diff"   ,data_df1,data_df1_r2,data_diff)

# getDiff_2df("calibCoeffFree"    ,"calibCoeffFreeErr"    ,"fcc_diff"  ,"fcce_diff"  ,data_df1,data_df2,data_diff)
# getDiff_2df("calibCoeffFree_cor","calibCoeffFreeErr_cor","fccc_diff" ,"fccce_diff" ,data_df1,data_df2,data_diff)

axisList = ["calibCoeff_cor","calibCoeffErr_cor","TotalSystematicError"]

print("DF")
printToScreen_listOfTypes(axisList,data_df1) 

print("DF (Run 2)")
printToScreen_listOfTypes(axisList,data_df1_r2) 

print("RH")
printToScreen_listOfTypes(axisList,data_df2) 

print("BL")
printToScreen_listOfTypes(axisList,data_df3)
# print data_df3.columns 

# print("DF") 
# printToScreen("calibCoeff","calibCoeffErr",data_df1) 
# print("RH") 
# printToScreen("calibCoeff","calibCoeffErr",data_df2) 

# print("DF, dB/dx (shim)")
# printToScreen("dBdx_shim","dBdx_shimErr",data_df1) 
# print("DF, dB/dy (shim)")
# printToScreen("dBdy_shim","dBdy_shimErr",data_df1) 
# print("DF, dB/dz (shim)")
# printToScreen("dBdz_shim","dBdz_shimErr",data_df1) 

# some stats
mean=0
stdev=0
mean,stdev = getStats("rh_ccc_diff",data_diff)
print("ccc_rhdf: mean = {0:.3f}, stdev = {1:.3f}".format(mean,stdev) )
mean,stdev = getStats("bl_ccc_diff",data_diff)
print("ccc_bldf: mean = {0:.3f}, stdev = {1:.3f}".format(mean,stdev) )

# general setup  
NCOL = 1
NROW = 2

tickSize      = 16
xAxisFontSize = 16
yAxisFontSize = 16

myFont = 'Helvetica'

legend = ['DF','RH']

if plotBingzhi:
   legend.append('BL')

if plotRun2: 
   legend.append("DF (Run 2)")

# sns.set_style("ticks")

# calib coeffs (NO FREE PROTON)
fig = plt.figure(1) 
plt.subplot(NROW,NCOL,1)

currentAxis = plt.gca() # grab current axis 
axis    = "calibCoeff_cor"
axisErr = "calibCoeffErr_cor"
data_df1.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[0], s=mSize, color=color[0], ax=currentAxis)
data_df2.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
if plotBingzhi:
   data_df3.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)

if plotRun2:
   data_df1_r2.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[3], s=mSize, color=color[3], ax=currentAxis)

# currentAxis.legend(legend)
currentAxis.set_xlabel("Probe"                   , fontsize=xAxisFontSize,fontname=myFont) 
currentAxis.set_ylabel("Calib Coeff (Hz)", fontsize=yAxisFontSize,fontname=myFont) 
currentAxis.set_ylim(bottom=100, top=300)
currentAxis.tick_params(labelsize=tickSize)

plt.subplot(NROW,NCOL,2)
currentAxis = plt.gca() # grab current axis 
axis    = "rh_ccc_diff"
axisErr = "rh_ccce_diff"
data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
if plotBingzhi:
   axis    = "bl_ccc_diff"
   axisErr = "bl_ccce_diff"
   data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)

if plotRun2:
   axis    = "r2_ccc_diff"
   axisErr = "r2_ccce_diff"
   data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[3], s=mSize, color=color[3], ax=currentAxis)

currentAxis.set_xlabel("Probe"                  , fontsize=xAxisFontSize,fontname=myFont) 
currentAxis.set_ylabel("Difference (Hz)", fontsize=yAxisFontSize,fontname=myFont) 
currentAxis.set_ylim(bottom=-4, top=12)
currentAxis.tick_params(labelsize=tickSize)

for ax in fig.get_axes():
    ax.label_outer()

# calib coeffs (NO FREE PROTON)
fig = plt.figure(2) 

currentAxis = plt.gca() # grab current axis 

diffRH = data_diff['rh_ccc_diff'].tolist()
plt.hist(diffRH,bins=50)
plt.title("Histogram of Differences (RH-DF)") 
plt.xlabel("RH-DF (Hz)",fontsize=xAxisFontSize,fontname=myFont) 
plt.ylabel("Counts"    ,fontsize=yAxisFontSize,fontname=myFont) 
plt.tick_params(labelsize=tickSize)

fig = plt.figure(3) 

currentAxis = plt.gca() # grab current axis 

diffBL = data_diff['bl_ccc_diff'].tolist()
plt.hist(diffBL,bins=50)
plt.title("Histogram of Differences (BL-DF)") 
plt.xlabel("BL-DF (Hz)",fontsize=xAxisFontSize,fontname=myFont) 
plt.ylabel("Counts"    ,fontsize=yAxisFontSize,fontname=myFont) 
plt.tick_params(labelsize=tickSize)

fig = plt.figure(4) 

currentAxis = plt.gca() # grab current axis 

diffBL = data_diff['r2_ccc_diff'].tolist()
plt.hist(diffBL,bins=50)
plt.title("Histogram of Differences (DF2-DF)") 
plt.xlabel("DF2-DF (Hz)",fontsize=xAxisFontSize,fontname=myFont) 
plt.ylabel("Counts"    ,fontsize=yAxisFontSize,fontname=myFont) 
plt.tick_params(labelsize=tickSize)

# data_diff.plot.hist(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
# if plotBingzhi:
#    axis    = "bl_cc_diff"
#    axisErr = "bl_cce_diff"
#    data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
# currentAxis.set_xlabel("Probe"                  , fontsize=xAxisFontSize,fontname=myFont) 
# currentAxis.set_ylabel("Blinded Difference (Hz)", fontsize=yAxisFontSize,fontname=myFont) 
# currentAxis.set_ylim(bottom=-4, top=12)
# currentAxis.tick_params(labelsize=tickSize)
# 
# for ax in fig.get_axes():
#     ax.label_outer()
 
 
# # plot calib coeffs against one another 
# fig = plt.figure(2) 
# 
# label_1 = "Calibration Coefficients (DF)"
# label_2 = "Calibration Coefficients (RH or BL)"
# 
# currentAxis = plt.gca() # grab current axis 
# data_comp.plot(kind="scatter", x="calibCoeff_cor_1", y="calibCoeff_cor_2", yerr="calibCoeff_cor_2Err", marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
# data_comp.plot(kind="scatter", x="calibCoeff_cor_1", y="calibCoeff_cor_3", yerr="calibCoeff_cor_3Err", marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
# currentAxis.set_xlabel(label_1, fontsize = xAxisFontSize) 
# currentAxis.set_ylabel(label_2, fontsize = yAxisFontSize)
# currentAxis.tick_params(labelsize=tickSize)
# 
# for ax in fig.get_axes():
#     ax.label_outer()

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
plt.tight_layout()
plt.show()

