#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Grid search algorithm that takes arrival time picks from S-files
and calculates 2D epicentres as the grid cell with lowest L2 residual.
"""

import datetime
import glob
import math
import numpy as np




def get_arrival_times(S_file):
    
    '''
    Parse S-file for P wave arrival times
    at recording stations.
    '''
    
    arrival_time_data = [[], []]
    with open(S_file) as openfile:
        
        r = -1
        for row in openfile:
            r += 1
            
            # Split row by spaces
            
            row_split = row[:-2].split(' ')
            
            # Get non-empty columns in the row
            
            cols = []
            for col in row_split:
                
                if len(col) != 0:
                    
                    cols.append(col)
                    
            # Collect data
                    
            if r == 0:
                
                # Get date
                
                year = cols[0]
                
                if len(cols[1]) > 2:
                    
                    month = cols[1][:-2]
                    day = cols[1][-2:]
                    
                else:
                    
                    
                    month= cols[1]
                    day = cols[2]
                
                if len(month) < 2:
                    month = '0' + month
                    
                if len(day) < 2:
                    day = '0' + day
                
                date = datetime.datetime.strptime(year + '-' + month + '-' + day,
                                                  '%Y-%m-%d')
                
            elif r > 3:
                
                if len(cols) == 0:
                    
                    continue
                
                # Get arrival time data (assume no day increments)
                
                station = cols[0][:-2] # Ignore component details
                
                seconds = cols[-1]
                
                if len(cols[-2]) > 2:
                    
                    minutes = cols[-2][-2:]
                    hours = cols[-2][:-2]
                    
                else:
                    
                    minutes = cols[-2]
                    hours = cols[-3]
                    
                total_seconds = (int(hours) * 3600 + int(minutes) * 60 + float(seconds))
                
                arrival_time_data[0].append(station)
                arrival_time_data[1].append(total_seconds)
                
    return arrival_time_data




def relative_timing(arrival_time_data):
    
    '''
    Make arrival time data relative for all stations.
    '''
    
    reference_station = arrival_time_data[0][0]
    reference_arrival_time = arrival_time_data[1][0] 
    
    relative_arrival_time_data = [[], []]
    for i in range(len(arrival_time_data[0])):
        
        relative_arrival_time_data[0].append(arrival_time_data[0][i])
        relative_arrival_time_data[1].append(arrival_time_data[1][i] - reference_arrival_time)
        
    return relative_arrival_time_data, reference_station




def generate_tt(xmin, xmax, ymin, ymax, xstep, ystep, velocity, station_positions):
    
    '''
    Generate 2D travel time grid.
    '''

    gridx, gridy = np.meshgrid(np.linspace(xmin, xmax, (xmax - xmin) / xstep + 1), 
                               np.linspace(ymin, ymax, (ymax - ymin) / ystep + 1))
    
    # Generate travel times from each station to each grid point
    
    travel_times = [[[] for i in range(len(gridx))] for j in range(len(gridx[0]))]
    
    for i in range(len(gridx)):
        for j in range(len(gridx[i])):
            
            x = gridx[i][j]
            y = gridy[i][j]
            
            for station_position in station_positions:
                
                distance = calculate_distance(x, y, station_position[0], station_position[1])
                
                travel_times[j][i].append(distance / velocity)
                
    return gridx, gridy, travel_times
    



def calculate_distance(x1, y1, x2, y2):
    
    '''
    Calculate distance between two points in a two-dimensional
    cartesian coordinate system.
    '''
    
    distance = math.sqrt((x2 - x1)**2 + (y2 - y1)**2)
    
    return distance




def relative_travel_time(travel_times, stations):
    
    '''
    Make travel times relative to a reference station.
    '''

    relative_travel_times = [[[[] for i in range(len(gridx))] for j in range(len(gridx[0]))] for k in range(len(stations))]

    for k in range(len(stations)):
        for j in range(len(travel_times)):
            for i in range(len(travel_times[j])):
                for l in range(len(travel_times[j][i])):
                    
                    relative_travel_times[k][j][i].append(travel_times[j][i][l] - travel_times[j][i][k])
            
    return relative_travel_times




def L2_residual(relative_arrival_times, relative_travel_times, stations):
    
    '''
    Calculate L2 residual for data and model output.
    '''
    
    L2_residual = 0
    
    for i in range(len(relative_arrival_times[0])):
        
        data_station = relative_arrival_times[0][i]
        
        for j in range(len(relative_travel_times)):
            
            model_station = stations[j]
            
            if data_station == model_station:
                
                # Scale by 1000 before squaring, to avoid square scaling reversals
                # with > 1, < 1 residuals.
                
                L2_residual += (1000000 * (relative_arrival_times[1][i] - relative_travel_times[j])**2)
        
    return math.sqrt(L2_residual)




# Get arrival time data

S_file_dir = '/home/sam/PAPERV2_SEISMIC_3COMP_EVENTS/'
S_files = glob.glob(S_file_dir + '*.S*')

events = []
for S_file in S_files:
    
    events.append(get_arrival_times(S_file))

# Define station names and positions (XYZ)

stations = ['TSNC1', 'TSNL2', 'TSNC3', 'TSNR3', 'TSNR2', 'TSNL3']

station_positions = [[0.00, 0.00],
                     [-0.53, 0.34],
                     [2.12, -0.11],
                     [1.01, -0.66],
                     [-0.77, -0.35],
                     [1.19, 0.43]]

# Set grid size and scale (km)

xmin = -4
xmax = 4
ymin = -3
ymax = 3
xstep = 0.05
ystep = 0.05

# Set grid velocity (km / s)

velocity = 1.65

# Calculate travel time grid and travel times

gridx, gridy, travel_times = generate_tt(xmin, xmax, ymin, ymax, xstep, ystep,
                                         velocity, station_positions)

# Calculate relative travel times for each station

relative_tt = relative_travel_time(travel_times, stations)

# Perform grid search

for e in range(len(events)):
    
    event = events[e]
    epicentre_idx = [0, 0]
    min_residual = 999999
    L2_residuals = [[[] for i in range(len(gridx))] for j in range(len(gridx[0]))]
    
    # Make arrival time relative to a given station

    relative_arrival_time_data, reference_station = relative_timing(event)
    
    # Find index of reference station in station list
    
    reference_station_idx = stations.index(reference_station)
    
    # Look through all grid cells
    
    for i in range(len(gridx)):
        
        for j in range(len(gridx[i])):
    
            # Make travel times relative to the reference station for the grid cell
    
            relative_travel_times = relative_tt[reference_station_idx][j][i]
            
            # Calculate the L2 residual for the grid cell
    
            L2_residuals[j][i] = L2_residual(relative_arrival_time_data, relative_travel_times, stations)
            
            if L2_residuals[j][i] < min_residual:
                
                min_residual = L2_residuals[j][i]
                epicentre_idx = [j, i] 
                epicentre = [gridx[i][j], gridy[i][j]]
                
    print(S_files[e].split('/')[-1] + ',' + str(min_residual) + ',' + str(epicentre_idx) + ',' + str(epicentre))
    np.save(S_files[e].split('/')[-1] + '_residuals.npy', L2_residuals)
            
            
    

