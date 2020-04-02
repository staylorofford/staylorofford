#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Take x timestep positions from geodetic quality GPS locations
(within a few metres of seismic GPS clocks) and calculate mean
horizontal and vertical position over a given period.
"""

import datetime
import glob
import math
import numpy as np
from scipy import stats

# Set start end dates (inclusive)

start_doy = 122
end_doy = 138

# Set file location

location_file_dir = '/home/sam/VEL4SAM/'
location_files = glob.glob(location_file_dir + 'tas_position*.txt')

# Parse relevant data

Eastings = [[] for i in range(len(location_files))]
Northings = [[] for i in range(len(location_files))]
Elevation = [[] for i in range(len(location_files))]
E_err = [[] for i in range(len(location_files))]
N_err = [[] for i in range(len(location_files))]
Elev_err = [[] for i in range(len(location_files))]

stations = []
station_position = [[] for i in range(len(location_files))]

l = -1
for location_file in location_files:
    l += 1
    
    stations.append(location_file.split('/')[-1][-8:-4])
    
    with open(location_file, 'r') as openfile:
        
        for row in openfile:
            
            init_cols = row.split(' ')
            
            # Remove empty "columns"
            
            cols = []
            for col in init_cols:
                
                if len(col) > 1:
                    cols.append(col)
            
            # Ignore no data periods
            
            if str(cols[1]) == 'NaN':
                
                continue
            
            # Extract data
            
            doy = datetime.datetime.fromordinal(int(math.floor(float(cols[0]) / 86400)) - 366).timetuple().tm_yday
            
            if start_doy <= doy <= end_doy:

                Eastings[l].append(float(cols[1]))
                Northings[l].append(float(cols[2]))
                Elevation[l].append(float(cols[3]))
                E_err[l].append(1 / (float(cols[4])**2))
                N_err[l].append(1 / (float(cols[5])**2))
                Elev_err[l].append(1 / (float(cols[6])**2))
                
# Calculate average position (weighted mean)   

for l in range(len(location_files)):
    
    if len(Eastings[l]) == 0:
        
        continue
    
    else:
    
        station_position[l] = [np.average(Eastings[l], weights = E_err[l]),
                                 np.average(Northings[l], weights = N_err[l]),
                                 np.average(Elevation[l], weights = Elev_err[l])]
    
# Return result to the operator
            
for l in range(len(location_files)):
    
    print(stations[l] + ',' + str(station_position[l][0]) + ',' + str(station_position[l][1]) + ',' + str(station_position[l][2]))