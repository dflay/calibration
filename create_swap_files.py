#! /usr/bin/python 

# Create swap data files for Ran
# Extract trolley swap positional data from analysis and print to file
# line format: probe,swap-1,...,swap-N 

import os
import sys
import pandas as pd 


prefix = "./output/blinded/flay/04-25-20/run-1/csv"

header = ["timestamp","freq","freqErr","temp","tempErr","x","xErr","y","yErr","phi","phiErr"]

header_out = ["probe"]
for i in xrange(0,9):
   header_out.append( "swap-{0:02d}".format(i+1) )
print header_out

pList   = []

# output dataframe
df_out = pd.DataFrame()

for i in xrange(0,17):
   # probe number 
   probe = i+1
   pList.append( "{0:02d}".format(probe) )
   # read in data
   inpath = "{0}/trly-swap-data_pr-{1:02d}.csv".format(prefix,probe) 
   df = pd.read_csv(inpath,names=header,index_col=False)
   # get phi data
   phi = df["phi"].tolist() 
   for entry in phi:
      swapStr = "{0:.3f}".format(entry)
      pList.append(swapStr)
   # add empty spots if necessary 
   N = len(pList)
   if(N<10):
      M = 10-N
      for j in xrange(0,M):
         pList.append("") 
   # check the row 
   print pList 
   # now add the row to a data frame
   row = pd.Series(pList) 
   # add to output dataframe 
   df_out = df_out.append(row,ignore_index=True) 
   # cleanup 
   pList[:] = [] 
   phi[:]   = []  

df_out.columns = header_out

print("Output data: ") 
print df_out 
