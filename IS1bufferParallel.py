### 
#!python

import subprocess
import multiprocessing

# set the corners of your ROI in geographic coordinates - this example is for Maryland
#lon_min = -180 * 1000
#lon_max = 180 * 1000
#lat_min = -52 * 1000
#lat_max = 52 * 1000

# list to store subsetting commands
cmd_list = []
#####Step 1: reduce grids number first
##### read from IS1 gridded at 0.125 degree csv data
import pandas as pd

# Read CSV file into a DataFrame
df = pd.read_csv('/gpfs/data1/vclgp/xiongl/HeightComparisonGEDI_IS2_IS1/Out/is1_global_0dot125degree_gridded.csv')

# Display the DataFrame
#print(df)

# Loop through every row
for index, row in df.iterrows():
        lon = row['x']
        lat = row['y']
        # pass necessary info to subsetting command
        #lon = -180 + i *0.125
        #lat = -52 + j * 0.125
        cmd = 'python IS1bufferFunction.py --lon %f  --lat %f ' % (lon, lat)
        cmd_list.append(cmd)

def os_quiet(cmd):
    try:
        #os.system(cmd)
        # Execute the command using subprocess
        subprocess.call(cmd, shell=True)
    except:
        pass

import time
import os
# set up parallel processing
#print(cmd_list)
nCPU = len(cmd_list)
if nCPU > 10: 
   nCPU = 10 # number of cores to use
   
pool = multiprocessing.Pool(nCPU) # Set up multi-processing
t0 = time.time() # start time
## Apply the process_input function to each input value in parallel

pool.map(os_quiet, cmd_list) # run multiple commands at once
elapsed_time = time.time() - t0 # end time
print(time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))






