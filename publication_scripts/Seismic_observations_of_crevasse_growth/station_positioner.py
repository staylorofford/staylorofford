#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Take hourly GPS position from GPS clocks and generate average
station position over a given period of time.
"""

import glob
import numpy as np
from scipy import stats

# Set start end dates (inclusive)

start_doy = 122
end_doy = 138

# Set file location

location_file_dir = '/home/sam/RISSIN_CSV/'
location_files = glob.glob(location_file_dir + '*NZTM.csv')

# Parse relevant data

Eastings = [[] for i in range(len(location_files))]
Northings = [[] for i in range(len(location_files))]

stations = []
station_position = [[] for i in range(len(location_files))]

l = -1
for location_file in location_files:
    l += 1
    
    stations.append(location_file.split('/')[-1][:5])
    
    with open(location_file, 'r') as openfile:
        
        for row in openfile:
            
            cols = row.split(',')
            
            doy = int(cols[0])
            
            if start_doy <= doy <= end_doy:
            
                Eastings[l].append(int(round(float(cols[-1][:cols[-1].index('\n')]))))
                Northings[l].append(int(round(float(cols[-2]))))
                
# Calculate average position (mode)    
            
for l in range(len(location_files)):
    
    station_position[l] = [stats.mode(Eastings[l])[0][0],
                            stats.mode(Northings[l])[0][0]]
    
#    station_position[l] = [np.mean(Eastings[l]), np.mean(Northings[l])]
    
#    station_position[l] = [np.median(Eastings[l]), np.median(Northings[l])]
    
# Return result to the operator
            
for l in range(len(location_files)):
    
    print(stations[l] + ',' + str(station_position[l][0]) + ',' + str(station_position[l][1]))
    
    
            
            

