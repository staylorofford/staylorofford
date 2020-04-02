#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Stack GRID2D results and plot them over a given time window.
"""

import datetime
import glob
import numpy as np
import math
import matplotlib.pyplot as plt
from bisect import bisect_left

def takeClosest(myList, myNumber):
    """
    Assumes myList is sorted. Returns closest value to myNumber.

    If two numbers are equally close, return the smallest number.
    From Lauritz V. Thaulow on stackoverflow
    """
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return myList[0], 0
    if pos == len(myList):
        return myList[-1], -1
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
       return after, pos
    else:
       return before, pos - 1




def distance(p1, p2):
    
    '''
    Calculate 2D distance between two points
    '''
    
    distance = math.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2)
    
    return distance

# Give parameters

## Numpy grid output location

event_output_directory = '/home/sam/EVENTS_IT3/TYPE_A/4/'

## Start and end dates for processing of day-long seismic streams

start_year = '2016'
start_month = '01'
start_day = '04'
start_hour = '00'

end_year = '2016'
end_month = '08'
end_day = '04'
end_hour = '23'

min_date = datetime.datetime.strptime(start_year + '-' + start_month + '-' + start_day + 'T' + start_hour,
                                      '%Y-%m-%dT%H')

max_date = datetime.datetime.strptime(end_year + '-' + end_month + '-' + end_day + 'T' + end_hour,
                                      '%Y-%m-%dT%H')

## Station positions

station_positions = [[1374.062, 5164.009, 0.821],
                     [1374.498, 5166.080, 0.904],
                     [1374.261, 5163.408, 0.816],
                     [1374.783, 5165.049, 0.857],
                     [1373.533, 5163.353, 0.801],
                     [1373.689, 5165.150, 0.857]]

## Grid shape

xmin = 1373.0
xmax = 1375.0
ymin = 5163.0
ymax = 5166.5
xstep = 0.05
ystep = 0.05

# Amount to shift position along bearing line every iteration
# Note: finer amounts are better, 
    
step = 0.01

## Location and uncertainty details

long_axis_direction = [0.1995, 0.9799]
bearing = math.atan(long_axis_direction[0] / long_axis_direction[1])
delta_E = math.sin(bearing) * step 
delta_N = math.cos(bearing) * step

grid_files = glob.glob(event_output_directory + '*.xcorrvaluegrid.npy')
grid_files.sort()

gridx, gridy = np.meshgrid(np.linspace(xmin, xmax, (xmax - xmin + 1) / xstep), np.linspace(ymin, ymax, (ymax - ymin + 1) / ystep))

# Create coordinate axis lists for matching values to
    
xvals = gridx[0].tolist()
yvals  = []
for x in gridy:
    yvals.append(x[0])

grid_values = np.zeros(gridx.shape)

for d in range(int(math.ceil((max_date - min_date).total_seconds()/86400))):
    
    locations = []
    last_date = 0
    
    new_min_date = min_date + datetime.timedelta(days = d)
    new_max_date = min_date + datetime.timedelta(days = d + 1)

    for grid_file in grid_files:
        
        # Apply date controls
        
        date = datetime.datetime.strptime(grid_file.split('/')[-1][:19], '%Y-%m-%dT%H:%M:%S')

        if new_min_date <= date <= new_max_date:
            pass
        elif date > new_max_date:
            break
        else:
            continue
        
#        last_date = date
        
        xcorr_value_grid = np.load(grid_file)
        
        max_value = np.max(xcorr_value_grid)
        max_value_indices = np.where(xcorr_value_grid == xcorr_value_grid.max())
        
        # Find location of max value (by convention, if two values are the same
        # take the one farthest from the origin)
        
        E = gridx[max_value_indices][-1]
        N = gridy[max_value_indices][-1]
        
        E_1 = E
        E_2 = E
        N_1 = N
        N_2 = N
        
        transect_data = [[max_value, E, N]]
        
        c = -1  
        while (xmin <= E_1 <= xmax) and (ymin <= N_1 <= ymax):
            c += 1
            
            # shift position forward and round to nearest grid cell
            
            E_1 = round(float(E + c * delta_E) / step) * step
            N_1 = round(float(N + c * delta_N) / step) * step
        
            # get the value for new coords from the xcorr value grid
            
            a = takeClosest(xvals, E_1)[-1]
            b = takeClosest(yvals, N_1)[-1]
            
            try:
                
                transect_data.append([xcorr_value_grid[b][a], xvals[a], yvals[b]])
                
            except:
                
                break
            
        # The same, but backwards
            
        c = 0
        while (xmin <= E_2 <= xmax) and (ymin <= N_2 <= ymax):
            c += 1
            
            E_2 = round(float(E - c * delta_E) / step) * step
            N_2 = round(float(N - c * delta_N) / step) * step
            a = takeClosest(xvals, E_2)[-1]
            b = takeClosest(yvals, N_2)[-1]
            
            try:
                
                transect_data.insert(0, [xcorr_value_grid[b][a], xvals[a], yvals[b]])
                
            except:
                
                break
            
        # Remove any duplicate values
        
        unique_transect_data = []
        
        for i in range(len(transect_data)):
            try:
                unique_transect_data.index(transect_data[i][1])
                unique_transect_data.index(transect_data[i][2])
            except:
                unique_transect_data.append(transect_data[i])
        
        transect_data = unique_transect_data
        
        # Calculate the second derivative of the transect values
        
        transect_second_derivative = []
        
        if len(transect_data) > 1:
            
            for i in range(1, len(transect_data) - 1):
                    
                # Calculate the second derivative using the second order
                # central difference and taking the separation as the mean
                # separation between the three points on the transect
                
                transect_second_derivative.append((abs(transect_data[i + 1][0] - 2 * 
                               transect_data[i][0] + transect_data[i - 1][0])
                / (np.mean([distance(transect_data[i][1:], transect_data[i - 1][1:]),
                            distance(transect_data[i][1:], transect_data[i + 1][1:])]))))
                
        # Handle edge points
                
        transect_second_derivative.insert(0, 0)
        transect_second_derivative.append(0)
        
        # Find plateau limits
        
        if len(transect_second_derivative) >= 7:
            
        # NOTE: this requirement means events with max values near the
        # grid edge do not produce a location. If the grid edge are outside
        # the seismic network, this is OK as events that locate outside the
        # network are generally network-external and thus any derived location
        # will be incorrect due to limitations in the cross-correlation method.
        
        # First find the points of max slope closest to the plateau centre        
                
            plateau_edge_indices = []
            first_slope_max = None
            second_slope_max = None
            for i in range(len(transect_data)):
                if (transect_data[i][1] == E) and (transect_data[i][2] == N):
                    transect_centre_index = i
            
            for i in range(len(transect_second_derivative)):
                
                # Move outward from the transect centre
                
                a = transect_centre_index - i
                b = transect_centre_index + i
                
                try:
                    if ((transect_second_derivative[a + 1] < transect_second_derivative[a] > transect_second_derivative[a - 1]) and
                        (a >= 0)):
                        if not first_slope_max:
                            first_slope_max = True
                            plateau_edge_indices.append(a)
                except:
                    pass
    
                try:                    
                    if ((transect_second_derivative[b - 1] < transect_second_derivative[b] > transect_second_derivative[b + 1]) and
                        (b >= 0)):
                        if not second_slope_max:
                            second_slope_max = True
                            plateau_edge_indices.append(b)
                except:
                    pass
        
            # Define the plateau as all transect_data within the plateau edge indices
            
    #        if len(plateau_edge_indices) == 1: 
    #            
    #            # If there is only one edge, assume the other is at the same distance
    #            
    #            plateau_edge_indices = [transect_centre_index - abs(transect_centre_index - plateau_edge_indices[0]),
    #                                    transect_centre_index + abs(transect_centre_index - plateau_edge_indices[0])]
    #            
    #            plateau_edge_indices.sort()
    #            
    #            if (plateau_edge_indices[0]) < 0:
    #                plateau_edge_indices[0] = 0
    #            
    #            if (plateau_edge_indices[1] > len(transect_data) - 1):
    #                plateau_edge_indices[1] = len(transect_data) - 1
                
            try:
                
                plateau_data = transect_data[min(plateau_edge_indices): max(plateau_edge_indices)]
                
            except:
                
                # This fails when no plateau edge can be defined
                
                continue
                
            # Best location is taken as the centre of the plateau
            # with uncertainty as the half width of the plateau
            
            try:
                # Let all single-edged (or no-edged) plateaus fail
                
                if plateau_edge_indices[0] != plateau_edge_indices[1]:
    
                    location = transect_data[int(round(np.mean(plateau_edge_indices)))][1:]
                    uncertainty = distance(transect_data[min(plateau_edge_indices)][1:],
                                        transect_data[max(plateau_edge_indices)][1:]) / 2
        
                else:
                
                    # Uncertainty is over the grid cell itself
                
                    location = [E, N]
                
                    uncertainty = distance(transect_data[min(plateau_edge_indices) - 1][1:],
                                         transect_data[max(plateau_edge_indices) + 1][1:]) / 4
                
                locations.append([location[0], location[1], uncertainty, max_value])
                
            except:
                
                continue
                                                     
    #        print(E, N)
    #        print(location, uncertainty, max_value)
    
#            if max_value < 0.7: continue
                        
            plt.imshow(xcorr_value_grid, extent = [xmin, xmax, ymin, ymax], origin = 'lower')
            plt.colorbar()
            for station_position in station_positions:
                plt.scatter(station_position[0], station_position[1], color = 'k')
            plt.xlim(xmin, xmax)
            plt.ylim(ymin, ymax)
            for data in transect_data[min(plateau_edge_indices) : max(plateau_edge_indices)]:
                plt.scatter(data[1], data[2], color ='green')
            plt.scatter(E, N, color = 'purple')
            plt.scatter(location[0], location[1], color = 'red')
            plt.savefig(grid_file[:-4] + '.png', format='png')
            plt.close()  
        
        
        #    date = datetime.datetime.strptime(grid_file.split('/')[-1][:13], '%Y-%m-%dT%H')
        #    if min_date <= date <= max_date:
        #
        #        xcorr_value_grid = np.load(grid_file)
        #        
        #        if np.max(xcorr_value_grid) >= 0:
        #    
        #            for i in range(len(xcorr_value_grid)):
        #                for j in range(len(xcorr_value_grid[i])):
        #    
        #                        grid_values[i][j] += xcorr_value_grid[i][j]
        
        #plt.contour(grid_values)            
        #plt.imshow(grid_values, extent = [xmin, xmax, ymin, ymax])
        #plt.colorbar()
        #for station_position in station_positions:
        #    plt.scatter(station_position[0], station_position[1], color = 'k')
            
        #plt.xlim(xmin, xmax)
        #plt.ylim(ymin, ymax)
        #plt.show()  
                
    # Save data in csv form
            
    outfile = open('location_file.csv', 'w')      
    outfile = open('location_file.csv', 'a')
    for location in locations:
        outfile.write(str(location[0]) + ',' + str(location[1]) + ',' + str(location[2]) + ',' + str(location[3]) + '\n')
        
    # Accumulate events in cells
        
    total_num = 0
    location_grid = [[0 for x in range(len(xvals))] for y in range(len(yvals))]
    for y in yvals:
        for x in xvals:
            for location in locations:
                if (location[0] == x) and (location[1] == y):
                    # X-corr value filtering
#                    if location[3] < 0.7: continue
                    total_num += 1
                    location_grid[yvals.index(y)][xvals.index(x)] += 1   
                    
#    print(last_date, total_num)
    
#    if last_date == 0: continue
     
    outfile = open(str(last_date) + '_location_grid.csv', 'w')      
    outfile = open(str(last_date) + '_location_grid.csv', 'a')               
    for y in yvals:
        for x in xvals:
            outfile.write(str(x) + ',' + str(y) +',' + str(location_grid[yvals.index(y)][xvals.index(x)]) + '\n')