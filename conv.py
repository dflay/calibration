# read in calibration coeffs from other analyzers and create a json object 

import csv
import json
import pandas as pd  

# whose data are we looking at 
usr       = "ran-hong" 
dataFile  = "run-1_04-03-20"

# test DF output 
# csv_path  = "/home/dflay/Workplace/dflay/calibration/output/blinded/flay/04-04-20/calibData_04-04-20.csv"

# create file paths
csv_path  = "./input/"          + usr + "/" + dataFile + ".csv" 
json_path = "./output/blinded/" + usr + "/" + dataFile + ".json"

# create a pandas dataframe, reading in the csv file  
print("Reading data from: {0}".format(csv_path)) 
df = pd.read_csv(csv_path,index_col=False) # index_col = False when you don't have an index column
print df

# create json data; start with python dictionary  
myDict = {} 

# convert columns to python lists
# then this list is the value of a given key (labeled as col here)   
for col in df.columns: 
   theList     = df[col].tolist()
   myDict[col] = theList 

# print the python dictionary to screen 
# print myDict

# print to a JSON file 
jsonFile = open(json_path,'w')
json.dump(myDict,jsonFile,indent=4)

print("JSON data printed to: {0}".format(json_path)) 
