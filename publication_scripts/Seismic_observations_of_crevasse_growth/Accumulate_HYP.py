#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sum hypocentre results from NLLoc (code pulled from Plot_GRID2D.py)
"""

import numpy as np

xmin = 1373.0
xmax = 1375.0
ymin = 5163.0
ymax = 5166.5
numx = 41
numy = 71

# Create coordinate axis lists for matching values to
    
xvals = np.linspace(xmin, xmax, numx).tolist()
yvals = np.linspace(ymin, ymax, numy).tolist()

locations = []
with open('/home/sam/NLLoc/nlloc_TGSN/HYP.csv', 'r') as openfile:
    for row in openfile:
        locations.append([])
        for col in row.split(','):
            locations[-1].append(float(col))

    # Accumulate events in cells
        
    total_num = 0
    location_grid = [[0 for x in range(len(xvals))] for y in range(len(yvals))]
    for y in yvals:
        for x in xvals:
            for location in locations:
                if (location[0] == x) and (location[1] == y):
                    total_num += 1
                    location_grid[yvals.index(y)][xvals.index(x)] += 1   
     
    outfile = open('location_grid.csv', 'w')      
    outfile = open('location_grid.csv', 'a')               
    for y in yvals:
        for x in xvals:
            outfile.write(str(x) + ',' + str(y) +',' + str(location_grid[yvals.index(y)][xvals.index(x)]) + '\n')