#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sum epicentre results from grid_search.py.
"""

import numpy as np

xmin = -3
xmax = 3
ymin = -4
ymax = 4
numx = 601
numy = 801

# Create coordinate axis lists for matching values to
    
xvals = np.linspace(xmin, xmax, numx).tolist()
yvals = np.linspace(ymin, ymax, numy).tolist()

locations = []
with open('/home/sam/GRIDSEARCH/gridsearch.csv', 'r') as openfile:
    for row in openfile:
        cols = row.split('[')
        locations.append([float(cols[-1].split(',')[0]), float(cols[-1].split(',')[1][:-2])])
                

# Accumulate events in cells
    
total_num = 0
location_grid = [[0 for x in range(len(xvals))] for y in range(len(yvals))]
for y in yvals:
    for x in xvals:
        for location in locations:
            if ((abs(location[0] - x) < 0.001) and (abs(location[1] - y) < 0.001)): 
                total_num += 1
                location_grid[yvals.index(y)][xvals.index(x)] += 1   
                
outfile = open('location_grid.csv', 'w')      
outfile = open('location_grid.csv', 'a')               
for y in yvals:
    for x in xvals:
        if location_grid[yvals.index(y)][xvals.index(x)] > 0:
            outfile.write(str(x) + ',' + str(y) +',' + str(location_grid[yvals.index(y)][xvals.index(x)]) + '\n')