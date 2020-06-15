# read in calibration coeffs from other analyzers and create a json object 

import csv
import json
import pandas as pd  

usr      = "bingzhi-li" 
dataFile = "run-2_03-06-20"

csv_path  = "./input/"          + usr + "/" + dataFile + ".csv" 
json_path = "./output/blinded/" + usr + "/" + dataFile + ".json"

# create a pandas dataframe, reading in the csv file  
print("Reading data from: {0}".format(csv_path)) 
df = pd.read_csv(csv_path)
print df.columns

# create json data; start with python dictionary  
myDict = {} 

# convert columns to python lists
# then this list is the value of a given key (labeled as col here)   
for col in df.columns: 
   theList     = df[col].tolist()
   myDict[col] = theList 

# print the python dictionary to screen 
print myDict

# print to a JSON file 
jsonFile = open(json_path,'w')
json.dump(myDict,jsonFile,indent=4)

print("Data printed to: {0}".format(json_path)) 
