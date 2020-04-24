# read in calibration coeffs from other analyzers and create a json object 

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
def getDiff_lists(x,xe,y,ye,z,ze):
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

plotBingzhi = False 

legend = ["DF","RH"]

if(plotBingzhi):
   legend.append("BL")  

# create file paths
csv_path_rh = "./input/ran-hong/run-1_04-18-20.csv" 
csv_path_bl = "./input/bingzhi-li/run-1_04-14-20.csv" 
csv_path_df = "./output/blinded/flay/04-22-20/run-1/calibData_04-22-20.csv"

# create a pandas dataframe, reading in the csv file  
print("Reading data from: {0}".format(csv_path_rh)) 
data_rh = pd.read_csv(csv_path_rh,index_col=False) # index_col = False when you don't have an index column
print("Reading data from: {0}".format(csv_path_bl)) 
data_bl = pd.read_csv(csv_path_bl,index_col=False) # index_col = False when you don't have an index column
print("Reading data from: {0}".format(csv_path_df)) 
data_df = pd.read_csv(csv_path_df,index_col=False) # index_col = False when you don't have an index column

# adjust bingzhi coefficients (he plots on absolute scale) 
data_bl['shim_x_c'] -= 61.78E+6
data_bl['shim_y_c'] -= 61.78E+6
data_bl['shim_z_c'] -= 61.78E+6

# marker parameters  
color  = ["blue","red","#20B010"]
mStyle = ["s","o","v"]
mSize  = 80

# axis details 
tickSize      = 16 
xAxisFontSize = 16 
yAxisFontSize = 16  

# calculate differences 
probeList = data_df['Probe'].tolist() # a list of the probe numbers to add to the new data frame 
data_diff = pd.DataFrame() # new dataframe of differences 
data_diff.insert(0,"Probe",probeList)

# take differences relative to DF
# for Ran
getDiff("calibCoeff"    ,"calibCoeffErr"    ,"cc_rhdf"  ,"cce_rhdf"  ,data_df,data_rh,data_diff)
getDiff("calibCoeff_cor","calibCoeffErr_cor","ccc_rhdf" ,"ccce_rhdf" ,data_df,data_rh,data_diff)
getDiff("misCor_x"      ,"misCor_xErr"      ,"mcx_rhdf" ,"mcxe_rhdf" ,data_df,data_rh,data_diff)
getDiff("misCor_y"      ,"misCor_yErr"      ,"mcy_rhdf" ,"mcye_rhdf" ,data_df,data_rh,data_diff)
getDiff("misCor_z"      ,"misCor_zErr"      ,"mcz_rhdf" ,"mcze_rhdf" ,data_df,data_rh,data_diff)
getDiff("mis_x"         ,"mis_xErr"         ,"mx_rhdf"  ,"mxe_rhdf"  ,data_df,data_rh,data_diff)
getDiff("mis_y"         ,"mis_yErr"         ,"my_rhdf"  ,"mye_rhdf"  ,data_df,data_rh,data_diff)
getDiff("mis_z"         ,"mis_zErr"         ,"mz_rhdf"  ,"mze_rhdf"  ,data_df,data_rh,data_diff)
getDiff("shim_x_a"      ,"shim_x_aErr"      ,"sxa_rhdf" ,"sxae_rhdf" ,data_df,data_rh,data_diff)
getDiff("shim_x_b"      ,"shim_x_bErr"      ,"sxb_rhdf" ,"sxbe_rhdf" ,data_df,data_rh,data_diff)
getDiff("shim_x_c"      ,"shim_x_cErr"      ,"sxc_rhdf" ,"sxce_rhdf" ,data_df,data_rh,data_diff)
getDiff("shim_y_a"      ,"shim_y_aErr"      ,"sya_rhdf" ,"syae_rhdf" ,data_df,data_rh,data_diff)
getDiff("shim_y_b"      ,"shim_y_bErr"      ,"syb_rhdf" ,"sybe_rhdf" ,data_df,data_rh,data_diff)
getDiff("shim_y_c"      ,"shim_y_cErr"      ,"syc_rhdf" ,"syce_rhdf" ,data_df,data_rh,data_diff)
getDiff("shim_z_a"      ,"shim_z_aErr"      ,"sza_rhdf" ,"szae_rhdf" ,data_df,data_rh,data_diff)
getDiff("shim_z_b"      ,"shim_z_bErr"      ,"szb_rhdf" ,"szbe_rhdf" ,data_df,data_rh,data_diff)
getDiff("shim_z_c"      ,"shim_z_cErr"      ,"szc_rhdf" ,"szce_rhdf" ,data_df,data_rh,data_diff)
getDiff("deltaB_pp_x"   ,"deltaB_pp_xErr"   ,"dbppx_rhdf" ,"dbppxe_rhdf" ,data_df,data_rh,data_diff)
getDiff("deltaB_pp_y"   ,"deltaB_pp_yErr"   ,"dbppy_rhdf" ,"dbppye_rhdf" ,data_df,data_rh,data_diff)
getDiff("deltaB_pp_z"   ,"deltaB_pp_zErr"   ,"dbppz_rhdf" ,"dbppze_rhdf" ,data_df,data_rh,data_diff)
getDiff("deltaB_tr_x"   ,"deltaB_tr_xErr"   ,"dbtrx_rhdf" ,"dbtrxe_rhdf" ,data_df,data_rh,data_diff)
getDiff("deltaB_tr_y"   ,"deltaB_tr_yErr"   ,"dbtry_rhdf" ,"dbtrye_rhdf" ,data_df,data_rh,data_diff)
getDiff("deltaB_tr_z"   ,"deltaB_tr_zErr"   ,"dbtrz_rhdf" ,"dbtrze_rhdf" ,data_df,data_rh,data_diff)
getDiff("dBdx_imp"      ,"dBdx_impErr"   ,"impx_rhdf" ,"impxe_rhdf" ,data_df,data_rh,data_diff)
getDiff("dBdy_imp"      ,"dBdy_impErr"   ,"impy_rhdf" ,"impye_rhdf" ,data_df,data_rh,data_diff)
getDiff("dBdz_imp"      ,"dBdz_impErr"   ,"impz_rhdf" ,"impze_rhdf" ,data_df,data_rh,data_diff)

# for Bingzhi 
getDiff("calibCoeff"    ,"calibCoeffErr"    ,"cc_bldf"  ,"cce_bldf"  ,data_df,data_bl,data_diff)
getDiff("calibCoeff_cor","calibCoeffErr_cor","ccc_bldf" ,"ccce_bldf" ,data_df,data_bl,data_diff)
getDiff("misCor_x"      ,"misCor_xErr"      ,"mcx_bldf" ,"mcxe_bldf" ,data_df,data_bl,data_diff)
getDiff("misCor_y"      ,"misCor_yErr"      ,"mcy_bldf" ,"mcye_bldf" ,data_df,data_bl,data_diff)
getDiff("misCor_z"      ,"misCor_zErr"      ,"mcz_bldf" ,"mcze_bldf" ,data_df,data_bl,data_diff)
getDiff("mis_x"         ,"mis_xErr"         ,"mx_bldf"  ,"mxe_bldf"  ,data_df,data_bl,data_diff)
getDiff("mis_y"         ,"mis_yErr"         ,"my_bldf"  ,"mye_bldf"  ,data_df,data_bl,data_diff)
getDiff("mis_z"         ,"mis_zErr"         ,"mz_bldf"  ,"mze_bldf"  ,data_df,data_bl,data_diff)
getDiff("shim_x_a"      ,"shim_x_aErr"      ,"sxa_bldf" ,"sxae_bldf" ,data_df,data_bl,data_diff)
getDiff("shim_x_b"      ,"shim_x_bErr"      ,"sxb_bldf" ,"sxbe_bldf" ,data_df,data_bl,data_diff)
getDiff("shim_x_c"      ,"shim_x_cErr"      ,"sxc_bldf" ,"sxce_bldf" ,data_df,data_bl,data_diff)
getDiff("shim_y_a"      ,"shim_y_aErr"      ,"sya_bldf" ,"syae_bldf" ,data_df,data_bl,data_diff)
getDiff("shim_y_b"      ,"shim_y_bErr"      ,"syb_bldf" ,"sybe_bldf" ,data_df,data_bl,data_diff)
getDiff("shim_y_c"      ,"shim_y_cErr"      ,"syc_bldf" ,"syce_bldf" ,data_df,data_bl,data_diff)
getDiff("shim_z_a"      ,"shim_z_aErr"      ,"sza_bldf" ,"szae_bldf" ,data_df,data_bl,data_diff)
getDiff("shim_z_b"      ,"shim_z_bErr"      ,"szb_bldf" ,"szbe_bldf" ,data_df,data_bl,data_diff)
getDiff("shim_z_c"      ,"shim_z_cErr"      ,"szc_bldf" ,"szce_bldf" ,data_df,data_bl,data_diff)
getDiff("deltaB_pp_x"   ,"deltaB_pp_xErr"   ,"dbppx_bldf" ,"dbppxe_bldf" ,data_df,data_bl,data_diff)
getDiff("deltaB_pp_y"   ,"deltaB_pp_yErr"   ,"dbppy_bldf" ,"dbppye_bldf" ,data_df,data_bl,data_diff)
getDiff("deltaB_pp_z"   ,"deltaB_pp_zErr"   ,"dbppz_bldf" ,"dbppze_bldf" ,data_df,data_bl,data_diff)
getDiff("deltaB_tr_x"   ,"deltaB_tr_xErr"   ,"dbtrx_bldf" ,"dbtrxe_bldf" ,data_df,data_bl,data_diff)
getDiff("deltaB_tr_y"   ,"deltaB_tr_yErr"   ,"dbtry_bldf" ,"dbtrye_bldf" ,data_df,data_bl,data_diff)
getDiff("deltaB_tr_z"   ,"deltaB_tr_zErr"   ,"dbtrz_bldf" ,"dbtrze_bldf" ,data_df,data_bl,data_diff)
getDiff("dBdx_imp"      ,"dBdx_impErr"      ,"impx_bldf"  ,"impxe_bldf"  ,data_df,data_bl,data_diff)
getDiff("dBdy_imp"      ,"dBdy_impErr"      ,"impy_bldf"  ,"impye_bldf"  ,data_df,data_bl,data_diff)
getDiff("dBdz_imp"      ,"dBdz_impErr"      ,"impz_bldf"  ,"impze_bldf"  ,data_df,data_bl,data_diff)

# print data_diff

# some stats
mean=0
stdev=0
mean,stdev = getStats("cc_rhdf",data_diff)
print("cc_rhdf: mean = {0:.3f}, stdev = {1:.3f}".format(mean,stdev) )
mean,stdev = getStats("cc_bldf",data_diff)
print("cc_bldf: mean = {0:.3f}, stdev = {1:.3f}".format(mean,stdev) )
mean,stdev = getStats("ccc_rhdf",data_diff)
print("ccc_rhdf: mean = {0:.3f}, stdev = {1:.3f}".format(mean,stdev) )
mean,stdev = getStats("ccc_bldf",data_diff)
print("ccc_bldf: mean = {0:.3f}, stdev = {1:.3f}".format(mean,stdev) )

# plot data

# calib coeffs (NO MISALIGNMENT)
NCOL = 1
NROW = 2
fig = plt.figure(1) 
plt.subplot(NROW,NCOL,1)

currentAxis = plt.gca() # grab current axis 
axis    = "calibCoeff"
axisErr = "calibCoeffErr"
data_df.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[0], s=mSize, color=color[0], ax=currentAxis)
data_rh.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
if plotBingzhi: 
   data_bl.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
# currentAxis.legend(legend)
currentAxis.set_xlabel("Probe"           ,fontsize=xAxisFontSize) 
currentAxis.set_ylabel("Calib Coeff (Hz)",fontsize=yAxisFontSize)
currentAxis.set_ylim(bottom=100, top=300)
currentAxis.tick_params(labelsize=tickSize) 

plt.subplot(NROW,NCOL,2)
currentAxis = plt.gca() # grab current axis 
axis    = "cc_rhdf"
axisErr = "cce_rhdf"
data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
axis    = "cc_bldf"                                                        
axisErr = "cce_bldf"                                                       
if plotBingzhi: 
   data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.set_xlabel("Probe"          ,fontsize=xAxisFontSize) 
currentAxis.set_ylabel("Difference (Hz)",fontsize=yAxisFontSize) 
currentAxis.set_ylim(bottom=-4, top=12)
currentAxis.tick_params(labelsize=tickSize) 

for ax in fig.get_axes():
    ax.label_outer()

# calib coeffs (WITH MISALIGNMENT)
fig = plt.figure(2) 
plt.subplot(NROW,NCOL,1)
currentAxis = plt.gca() # grab current axis 
axis    = "calibCoeff_cor"
axisErr = "calibCoeffErr_cor"
data_df.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[0], s=mSize, color=color[0], ax=currentAxis)
data_rh.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
if plotBingzhi: 
   data_bl.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
# currentAxis.legend(["DF","RH","BL"])
currentAxis.set_xlabel("Probe"                 ,fontsize=xAxisFontSize) 
currentAxis.set_ylabel("Calib Coeff [Cor] (Hz)",fontsize=yAxisFontSize)
currentAxis.set_ylim(bottom=100, top=300)
currentAxis.tick_params(labelsize=tickSize) 
 
plt.subplot(NROW,NCOL,2)
currentAxis = plt.gca() # grab current axis 
axis    = "ccc_rhdf"
axisErr = "ccce_rhdf"
data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
axis    = "ccc_bldf"
axisErr = "ccce_bldf"
if plotBingzhi: 
   data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.set_xlabel("Probe"                      ,fontsize=xAxisFontSize) 
currentAxis.set_ylabel("Calib Coeff [Cor] Diff (Hz)",fontsize=yAxisFontSize) 
currentAxis.set_ylim(bottom=-4, top=12)
currentAxis.tick_params(labelsize=tickSize) 

for ax in fig.get_axes():
    ax.label_outer()

# misalignment corrections: x
title = ["misCor_x","misCor_y","misCor_z"]
fig = plt.figure(3)
NROW = 2
NCOL = 3
plt.subplot(NROW,NCOL,1)
currentAxis = plt.gca() # grab current axis 
axis    = "misCor_x"
axisErr = "misCor_xErr"
data_df.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[0], marker=mStyle[0], s=mSize, color=color[0], ax=currentAxis)
data_rh.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[0], marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
if plotBingzhi: 
   data_bl.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[0], marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.legend(legend)
currentAxis.set_xlabel("") 
currentAxis.set_ylabel("misCor (Hz)") 

plt.subplot(NROW,NCOL,4)
currentAxis = plt.gca() # grab current axis 
axis    = "mcx_rhdf"
axisErr = "mcxe_rhdf"
data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
axis    = "mcx_bldf"
axisErr = "mcxe_bldf"
if plotBingzhi: 
   data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.set_xlabel("Probe") 
currentAxis.set_ylabel("misCor Diff (Hz)") 

# misalignment corrections: y 
plt.subplot(NROW,NCOL,2)
currentAxis = plt.gca() # grab current axis 
axis    = "misCor_y"
axisErr = "misCor_yErr"
data_df.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[1], marker=mStyle[0], s=mSize, color=color[0], ax=currentAxis)
data_rh.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[1], marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
if plotBingzhi: 
   data_bl.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[1], marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.legend(legend)
currentAxis.set_xlabel("") 
currentAxis.set_ylabel("") 

plt.subplot(NROW,NCOL,5)
currentAxis = plt.gca() # grab current axis 
axis    = "mcy_rhdf"
axisErr = "mcye_rhdf"
data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
axis    = "mcy_bldf"
axisErr = "mcye_bldf"
if plotBingzhi: 
   data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.set_xlabel("Probe") 
currentAxis.set_ylabel("") 

# misalignment corrections: z 
plt.subplot(NROW,NCOL,3)
currentAxis = plt.gca() # grab current axis 
axis    = "misCor_z"
axisErr = "misCor_zErr"
data_df.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[2], marker=mStyle[0], s=mSize, color=color[0], ax=currentAxis)
data_rh.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[2], marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
if plotBingzhi: 
   data_bl.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[2], marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.legend(legend)
currentAxis.set_xlabel("") 
currentAxis.set_ylabel("") 

plt.subplot(NROW,NCOL,6)
currentAxis = plt.gca() # grab current axis 
axis    = "mcz_rhdf"
axisErr = "mcze_rhdf"
data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
axis    = "mcz_bldf"
axisErr = "mcze_bldf"
if plotBingzhi: 
   data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.set_xlabel("Probe") 
currentAxis.set_ylabel("") 

# shimmed gradient fits, a coeff
fig = plt.figure(4)
NROW = 2
NCOL = 3
title = ["shim_x_a","shim_y_a","shim_z_a"]
plt.subplot(NROW,NCOL,1)
currentAxis = plt.gca() # grab current axis 
axis    = "shim_x_a"
axisErr = "shim_x_aErr"
data_df.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[0], marker=mStyle[0], s=mSize, color=color[0], ax=currentAxis)
data_rh.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[0], marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
if plotBingzhi: 
   data_bl.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[0], marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.legend(legend)
currentAxis.set_xlabel("") 
currentAxis.set_ylabel("shim a Par (Hz/mm^2)") 

plt.subplot(NROW,NCOL,4)
currentAxis = plt.gca() # grab current axis 
axis    = "sxa_rhdf"
axisErr = "sxae_rhdf"
data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
axis    = "sxa_bldf"
axisErr = "sxae_bldf"
if plotBingzhi: 
   data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.set_xlabel("Probe") 
currentAxis.set_ylabel("shim a Par Diff (Hz/mm^2)") 

# shimmed field gradient a: y 
plt.subplot(NROW,NCOL,2)
currentAxis = plt.gca() # grab current axis 
axis    = "shim_y_a"
axisErr = "shim_y_aErr"
data_df.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[1], marker=mStyle[0], s=mSize, color=color[0], ax=currentAxis)
data_rh.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[1], marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
if plotBingzhi: 
   data_bl.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[1], marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.legend(legend)
currentAxis.set_xlabel("") 
currentAxis.set_ylabel("") 

plt.subplot(NROW,NCOL,5)
currentAxis = plt.gca() # grab current axis 
axis    = "sya_rhdf"
axisErr = "syae_rhdf"
data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
axis    = "sya_bldf"
axisErr = "syae_bldf"
if plotBingzhi: 
   data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.set_xlabel("Probe") 
currentAxis.set_ylabel("") 

# shimmed field gradient a: z  
plt.subplot(NROW,NCOL,3)
currentAxis = plt.gca() # grab current axis 
axis    = "shim_z_a"
axisErr = "shim_z_aErr"
data_df.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[2], marker=mStyle[0], s=mSize, color=color[0], ax=currentAxis)
data_rh.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[2], marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
if plotBingzhi: 
   data_bl.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[2], marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.legend(legend)
currentAxis.set_xlabel("") 
currentAxis.set_ylabel("") 

plt.subplot(NROW,NCOL,6)
currentAxis = plt.gca() # grab current axis 
axis    = "sza_rhdf"
axisErr = "szae_rhdf"
data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
axis    = "sza_bldf"
axisErr = "szae_bldf"
if plotBingzhi: 
   data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.set_xlabel("Probe") 
currentAxis.set_ylabel("") 

# shimmed field gradient b   
title = ["shim_x_b","shim_y_b","shim_z_b"]
fig = plt.figure(5)
NROW = 2
NCOL = 3
plt.subplot(NROW,NCOL,1)
currentAxis = plt.gca() # grab current axis 
axis    = "shim_x_b"
axisErr = "shim_x_bErr"
data_df.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[0], marker=mStyle[0], s=mSize, color=color[0], ax=currentAxis)
data_rh.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[0], marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
if plotBingzhi: 
   data_bl.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[0], marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.legend(legend)
currentAxis.set_xlabel("") 
currentAxis.set_ylabel("shim b Par (Hz/mm)") 

plt.subplot(NROW,NCOL,4)
currentAxis = plt.gca() # grab current axis 
axis    = "sxb_rhdf"
axisErr = "sxbe_rhdf"
data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
axis    = "sxb_bldf"
axisErr = "sxbe_bldf"
if plotBingzhi: 
   data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.set_xlabel("Probe") 
currentAxis.set_ylabel("shim b Par Diff (Hz/mm)") 

# shimmed field gradient b: y  
plt.subplot(NROW,NCOL,2)
currentAxis = plt.gca() # grab current axis 
axis    = "shim_y_b"
axisErr = "shim_y_bErr"
data_df.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[1], marker=mStyle[0], s=mSize, color=color[0], ax=currentAxis)
data_rh.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[1], marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
if plotBingzhi: 
   data_bl.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[1], marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.legend(legend)
currentAxis.set_xlabel("") 
currentAxis.set_ylabel("") 

plt.subplot(NROW,NCOL,5)
currentAxis = plt.gca() # grab current axis 
axis    = "syb_rhdf"
axisErr = "sybe_rhdf"
data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
axis    = "syb_bldf"
axisErr = "sybe_bldf"
if plotBingzhi: 
  data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.set_xlabel("Probe") 
currentAxis.set_ylabel("") 

# shimmed field gradient b: z  
plt.subplot(NROW,NCOL,3)
currentAxis = plt.gca() # grab current axis 
axis    = "shim_z_b"
axisErr = "shim_z_bErr"
data_df.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[2], marker=mStyle[0], s=mSize, color=color[0], ax=currentAxis)
data_rh.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[2], marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
if plotBingzhi: 
   data_bl.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[2], marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.legend(legend)
currentAxis.set_xlabel("") 
currentAxis.set_ylabel("") 

plt.subplot(NROW,NCOL,6)
currentAxis = plt.gca() # grab current axis 
axis    = "szb_rhdf"
axisErr = "szbe_rhdf"
data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
axis    = "szb_bldf"
axisErr = "szbe_bldf"
if plotBingzhi: 
   data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.set_xlabel("Probe") 
currentAxis.set_ylabel("") 

# shimmed gradient fits c: x
title = ["shim_x_c","shim_y_c","shim_z_c"]
fig = plt.figure(6)
NROW = 2
NCOL = 3
plt.subplot(NROW,NCOL,1)
currentAxis = plt.gca() # grab current axis 
axis    = "shim_x_c"
axisErr = "shim_x_cErr"
data_df.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[0], marker=mStyle[0], s=mSize, color=color[0], ax=currentAxis)
data_rh.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[0], marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
if plotBingzhi: 
   data_bl.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[0], marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.legend(legend)
currentAxis.set_xlabel("") 
currentAxis.set_ylabel("shim c Par (Hz)") 

plt.subplot(NROW,NCOL,4)
currentAxis = plt.gca() # grab current axis 
axis    = "sxc_rhdf"
axisErr = "sxce_rhdf"
data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
axis    = "sxc_bldf"
axisErr = "sxce_bldf"
if plotBingzhi: 
   data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.set_xlabel("Probe") 
currentAxis.set_ylabel("shim c Par Diff (Hz)") 

# shimmed field gradient c: y 
plt.subplot(NROW,NCOL,2)
currentAxis = plt.gca() # grab current axis 
axis    = "shim_y_c"
axisErr = "shim_y_cErr"
data_df.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[1], marker=mStyle[0], s=mSize, color=color[0], ax=currentAxis)
data_rh.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[1], marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
if plotBingzhi: 
   data_bl.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[1], marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.legend(legend)
currentAxis.set_xlabel("") 
currentAxis.set_ylabel("") 

plt.subplot(NROW,NCOL,5)
currentAxis = plt.gca() # grab current axis 
axis    = "syc_rhdf"
axisErr = "syce_rhdf"
data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
axis    = "syc_bldf"
axisErr = "syce_bldf"
if plotBingzhi: 
   data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.set_xlabel("Probe") 
currentAxis.set_ylabel("") 

# shimmed field gradient c: z
plt.subplot(NROW,NCOL,3)
currentAxis = plt.gca() # grab current axis 
axis    = "shim_z_c"
axisErr = "shim_z_cErr"
data_df.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[2], marker=mStyle[0], s=mSize, color=color[0], ax=currentAxis)
data_rh.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[2], marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
if plotBingzhi: 
   data_bl.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[2], marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.legend(legend)
currentAxis.set_xlabel("") 
currentAxis.set_ylabel("") 

plt.subplot(NROW,NCOL,6)
currentAxis = plt.gca() # grab current axis 
axis    = "szc_rhdf"
axisErr = "szce_rhdf"
data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
axis    = "szc_bldf"
axisErr = "szce_bldf"
if plotBingzhi: 
   data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.set_xlabel("Probe") 
currentAxis.set_ylabel("") 

# deltaB PP  
title = ["deltaB_pp_x","deltaB_pp_y","deltaB_pp_z"]
fig = plt.figure(7)
NROW = 2
NCOL = 3
plt.subplot(NROW,NCOL,1)
currentAxis = plt.gca() # grab current axis 
axis    = "deltaB_pp_x"
axisErr = "deltaB_pp_xErr"
data_df.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[0], marker=mStyle[0], s=mSize, color=color[0], ax=currentAxis)
data_rh.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[0], marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
if plotBingzhi: 
   data_bl.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[0], marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.legend(legend)
currentAxis.set_xlabel("") 
currentAxis.set_ylabel("DeltaB (Hz)") 

plt.subplot(NROW,NCOL,4)
currentAxis = plt.gca() # grab current axis 
axis    = "dbppx_rhdf"
axisErr = "dbppxe_rhdf"
data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
axis    = "dbppx_bldf"
axisErr = "dbppxe_bldf"
if plotBingzhi: 
   data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.set_xlabel("Probe") 
currentAxis.set_ylabel("DeltaB Diff (Hz)") 

# dBy 
plt.subplot(NROW,NCOL,2)
currentAxis = plt.gca() # grab current axis 
axis    = "deltaB_pp_y"
axisErr = "deltaB_pp_yErr"
data_df.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[1], marker=mStyle[0], s=mSize, color=color[0], ax=currentAxis)
data_rh.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[1], marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
if plotBingzhi: 
   data_bl.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[1], marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.legend(legend)
currentAxis.set_xlabel("") 
currentAxis.set_ylabel("") 

plt.subplot(NROW,NCOL,5)
currentAxis = plt.gca() # grab current axis 
axis    = "dbppy_rhdf"
axisErr = "dbppye_rhdf"
data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
axis    = "dbppy_bldf"
axisErr = "dbppye_bldf"
if plotBingzhi: 
   data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.set_xlabel("Probe") 
currentAxis.set_ylabel("") 

# dBz 
plt.subplot(NROW,NCOL,3)
currentAxis = plt.gca() # grab current axis 
axis    = "deltaB_pp_z"
axisErr = "deltaB_pp_zErr"
data_df.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[2], marker=mStyle[0], s=mSize, color=color[0], ax=currentAxis)
data_rh.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[2], marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
if plotBingzhi: 
   data_bl.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[2], marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.legend(legend)
currentAxis.set_xlabel("") 
currentAxis.set_ylabel("") 

plt.subplot(NROW,NCOL,6)
currentAxis = plt.gca() # grab current axis 
axis    = "dbppz_rhdf"
axisErr = "dbppze_rhdf"
data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
axis    = "dbppz_bldf"
axisErr = "dbppze_bldf"
if plotBingzhi: 
   data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.set_xlabel("Probe") 
currentAxis.set_ylabel("") 

# deltaB TR  
title = ["deltaB_tr_x","deltaB_tr_y","deltaB_tr_z"]
fig = plt.figure(8)
NROW = 2
NCOL = 3
plt.subplot(NROW,NCOL,1)
currentAxis = plt.gca() # grab current axis 
axis    = "deltaB_tr_x"
axisErr = "deltaB_tr_xErr"
data_df.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[0], marker=mStyle[0], s=mSize, color=color[0], ax=currentAxis)
data_rh.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[0], marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
if plotBingzhi: 
   data_bl.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[0], marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.legend(legend)
currentAxis.set_xlabel("") 
currentAxis.set_ylabel("DeltaB (Hz)") 

plt.subplot(NROW,NCOL,4)
currentAxis = plt.gca() # grab current axis 
axis    = "dbtrx_rhdf"
axisErr = "dbtrxe_rhdf"
data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
axis    = "dbtrx_bldf"
axisErr = "dbtrxe_bldf"
if plotBingzhi: 
   data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.set_xlabel("Probe") 
currentAxis.set_ylabel("DeltaB Diff (Hz)") 

# dBy 
plt.subplot(NROW,NCOL,2)
currentAxis = plt.gca() # grab current axis 
axis    = "deltaB_tr_y"
axisErr = "deltaB_tr_yErr"
data_df.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[1], marker=mStyle[0], s=mSize, color=color[0], ax=currentAxis)
data_rh.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[1], marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
if plotBingzhi: 
   data_bl.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[1], marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.legend(legend)
currentAxis.set_xlabel("") 
currentAxis.set_ylabel("") 

plt.subplot(NROW,NCOL,5)
currentAxis = plt.gca() # grab current axis 
axis    = "dbtry_rhdf"
axisErr = "dbtrye_rhdf"
data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
axis    = "dbtry_bldf"
axisErr = "dbtrye_bldf"
if plotBingzhi: 
   data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.set_xlabel("Probe") 
currentAxis.set_ylabel("") 

# dBz 
plt.subplot(NROW,NCOL,3)
currentAxis = plt.gca() # grab current axis 
axis    = "deltaB_tr_z"
axisErr = "deltaB_tr_zErr"
data_df.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[2], marker=mStyle[0], s=mSize, color=color[0], ax=currentAxis)
data_rh.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[2], marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
if plotBingzhi: 
   data_bl.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[2], marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.legend(legend)
currentAxis.set_xlabel("") 
currentAxis.set_ylabel("") 

plt.subplot(NROW,NCOL,6)
currentAxis = plt.gca() # grab current axis 
axis    = "dbtrz_rhdf"
axisErr = "dbtrze_rhdf"
data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
axis    = "dbtrz_bldf"
axisErr = "dbtrze_bldf"
if plotBingzhi: 
   data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.set_xlabel("Probe") 
currentAxis.set_ylabel("") 

# imposed gradients  
title = ["dBdx_imp","dBdy_imp","dBdz_imp"]
fig = plt.figure(9)
NROW = 2
NCOL = 3
plt.subplot(NROW,NCOL,1)
currentAxis = plt.gca() # grab current axis 
axis    = "dBdx_imp"
axisErr = "dBdx_impErr"
data_df.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[0], marker=mStyle[0], s=mSize, color=color[0], ax=currentAxis)
data_rh.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[0], marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
if plotBingzhi: 
   data_bl.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[0], marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.legend(legend)
currentAxis.set_xlabel("") 
currentAxis.set_ylabel("dB/dq [imp] (Hz/mm)") 

plt.subplot(NROW,NCOL,4)
currentAxis = plt.gca() # grab current axis 
axis    = "impx_rhdf"
axisErr = "impxe_rhdf"
data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
axis    = "impx_bldf"
axisErr = "impxe_bldf"
if plotBingzhi: 
   data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.set_xlabel("Probe") 
currentAxis.set_ylabel("dB/dq Diff (Hz/mm)") 

# y 
plt.subplot(NROW,NCOL,2)
currentAxis = plt.gca() # grab current axis 
axis    = "dBdy_imp"
axisErr = "dBdy_impErr"
data_df.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[1], marker=mStyle[0], s=mSize, color=color[0], ax=currentAxis)
data_rh.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[1], marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
if plotBingzhi: 
   data_bl.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[1], marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.legend(legend)
currentAxis.set_xlabel("") 
currentAxis.set_ylabel("") 

plt.subplot(NROW,NCOL,5)
currentAxis = plt.gca() # grab current axis 
axis    = "impy_rhdf"
axisErr = "impye_rhdf"
data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
axis    = "impy_bldf"
axisErr = "impye_bldf"
if plotBingzhi: 
   data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.set_xlabel("Probe") 
currentAxis.set_ylabel("") 

# dBz 
plt.subplot(NROW,NCOL,3)
currentAxis = plt.gca() # grab current axis 
axis    = "dBdz_imp"
axisErr = "dBdz_impErr"
data_df.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[2], marker=mStyle[0], s=mSize, color=color[0], ax=currentAxis)
data_rh.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[2], marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
if plotBingzhi: 
   data_bl.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, title=title[2], marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.legend(legend)
currentAxis.set_xlabel("") 
currentAxis.set_ylabel("") 

plt.subplot(NROW,NCOL,6)
currentAxis = plt.gca() # grab current axis 
axis    = "impz_rhdf"
axisErr = "impze_rhdf"
data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[1], s=mSize, color=color[1], ax=currentAxis)
axis    = "impz_bldf"
axisErr = "impze_bldf"
if plotBingzhi: 
   data_diff.plot(kind="scatter", x="Probe", y=axis, yerr=axisErr, marker=mStyle[2], s=mSize, color=color[2], ax=currentAxis)
currentAxis.set_xlabel("Probe") 
currentAxis.set_ylabel("") 

# show all plots
plt.show()

