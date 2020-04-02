#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot residual grids from grid_search.py
"""

import glob
import numpy as np
import matplotlib.pyplot as plt

residual_grid_dir = '/home/sam/GRIDSEARCH/'
residual_grids = glob.glob(residual_grid_dir + '*npy')

station_positions = [[0.00, 0.00],
                     [-0.53, 0.34],
                     [2.12, -0.11],
                     [1.01, -0.66],
                     [-0.77, -0.35],
                     [1.19, 0.43]]

xmin = -3
xmax = 3
ymin = -4
ymax = 4

for residual_grid in residual_grids:
    
#    grid = np.load(residual_grid)
    
#    plt.imshow(grid, origin = 'lower', extent = [xmin, xmax, ymin, ymax])
#    plt.colorbar()
    
    with open(residual_grid_dir + 'gridsearch.csv', 'r') as openfile:
        
        for row in openfile:
            
            if row[:19] == residual_grid.split('/')[-1][:19]:
                
                cols = row.split('[')
                
                y = float(cols[-1].split(',')[0])
                x = float(cols[-1].split(',')[1][:-2])
                
                plt.scatter(x, y, color = 'red')
    
    for station in station_positions:
        
        plt.scatter(station[0], station[1], color = 'k')
    
plt.show()
