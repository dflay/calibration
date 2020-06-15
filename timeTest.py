
import datetime 
import time 

from myfuncs import utcNow 

theTime_0 = utcNow()

print("Waiting...")
time.sleep(10) 

theTime_1 = utcNow() 
diff = float(theTime_1 - theTime_0) 

print("Elapsed time: {0:.3f} sec".format(diff)) 

