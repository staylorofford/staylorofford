#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fit linear plane to point cloud using least squares,
then translate horizontal station positions vertically onto plane,
then transform 3D coordinates into an in-plane 2D coordinate system.

Plane fitting adjusted from https://gist.github.com/amroamroamro/1db8d69b4b65e8bc66a6.
3D -> 2D translation from https://stackoverflow.com/questions/23814234/convert-3d-plane-to-2d.
"""

import math
import numpy as np
import scipy.linalg
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

# Load horizontal station positions from csv file

station_position_file = '/home/sam/RISSIN_CSV/mode_clock_122_138_pos.csv'
stations_to_plot = ['TSNC1', 'TSNC3', 'TSNL2', 'TSNL3', 'TSNR2', 'TSNR3']

sta, sta_x, sta_y, sta_z = [[], [], [], []] 
with open(station_position_file, 'r') as openfile:
    
    for row in openfile:
        
        cols = row.split(',')
        
        if cols[0] in stations_to_plot:
        
            sta.append(cols[0])
            sta_x.append(float(cols[1]))
            sta_y.append(float(cols[2]))
            sta_z.append(float(cols[-1]))
            
# Load geodetic station positions from csv file

g_station_position_file = '/home/sam/RISSIN_CSV/geodetic_122_138_pos.csv'
g_stations_to_plot = ['TSNC1', 'TSNC3']

g_sta, g_sta_x, g_sta_y, g_sta_z = [[], [], [], []] 
with open(g_station_position_file, 'r') as openfile:
    
    for row in openfile:
        
        cols = row.split(',')
        
        if cols[0] in g_stations_to_plot:
        
            g_sta.append(cols[0])
            g_sta_x.append(float(cols[1]))
            g_sta_y.append(float(cols[2]))
            g_sta_z.append(float(cols[3]))

# Load glacier surface points from ASCII xyz file

XYZ_file = '/home/sam/RISSIN_CSV/lower_tasman_elevation.xyz'

x, y, z = [[], [], []]
with open(XYZ_file, 'r') as openfile:
    
    for row in openfile:
        
        cols = row.split(' ')
        
        # Conditional to ignore no data points
        # and points from Ball Glacier
        
        if float(cols[2][:4]) > 1000: 
            
            continue
    
        else:
            
            x.append(float(cols[0]))
            y.append(float(cols[1]))
            z.append(float(cols[2]))
        
data = np.c_[x,y,z]

# Best-fit linear plane to data

A = np.c_[data[:,0], data[:,1], np.ones(data.shape[0])]
C,_,_,_ = scipy.linalg.lstsq(A, data[:,2])    # coefficients

# Generate plotting dataset of plane over grid

mn = np.min(data, axis=0)
mx = np.max(data, axis=0)
X,Y = np.meshgrid(np.linspace(mn[0], mx[0], 20), np.linspace(mn[1], mx[1], 20))

Z = C[0]*X + C[1]*Y + C[2]

print('Plane parameters are:')
print('a = ' + str(C[0]) + ', b = ' + str(C[1]) + ', c = ' + str(C[2]))
print('\n')

# Generate station positions on plane by translating them vertically

plane_sta_z = []
for i in range(len(sta)):
    
    plane_sta_z.append(C[0] * sta_x[i] + C[1] * sta_y[i] + C[2])
    
# Calculate 3D station distances to test later 3D -> 2D translation of points

sta_distances_3D = [[] for i in range(len(sta))]
for i in range(len(sta)):
    
    for j in range(len(sta)):
        
        sta_distances_3D[i].append(math.sqrt((sta_x[j] - sta_x[i])**2 +
                        (sta_y[j] - sta_y[i])**2 + (plane_sta_z[j] - plane_sta_z[i])**2))

# Calculate 2D axes X, Y in plane

z_plane = [C[0], C[1], C[2]]
z_plane /= np.linalg.norm(z_plane)

U = [sta_x[sta.index('TSNC1')], sta_y[sta.index('TSNC1')], plane_sta_z[sta.index('TSNC1')]]
O = U
U /= np.linalg.norm(U)

x_plane = [U[i] - np.dot(U, z_plane) * z_plane[i] for i in range(len(z_plane))]
x_plane /= np.linalg.norm(x_plane)

y_plane = np.cross(z_plane, x_plane)

print('Plane origin (arbitrarily set) is:')
print('O = [' + str(O[0]) + ',' + str(O[1]) + ',' + str(O[2]) + ']')
print('\n')

print('In-plane coordinate vectors are:')
print('x = [' + str(x_plane[0]) + ',' + str(x_plane[1]) + ',' + str(x_plane[2]) + ']')
print('y = [' + str(y_plane[0]) + ',' + str(y_plane[1]) + ',' + str(y_plane[2]) + ']')
print('z = [' + str(z_plane[0]) + ',' + str(z_plane[1]) + ',' + str(z_plane[2]) + ']')
print('\n')

print('In-plane coordinate vector dot products are:')
print('x cross z = ' + str(np.dot(x_plane, z_plane)))
print('y cross z = ' + str(np.dot(y_plane, z_plane)))
print('x cross y = ' + str(np.dot(x_plane, y_plane)))
print('\n')

# Translate DEM points onto plane for use as backdrop

Po = [np.dot(x_plane, O), np.dot(y_plane, O)]

x_2D, y_2D = [], []
for i in range(len(x)):
    
    Pi = [x[i], y[i], C[0] * x[i] + C[1] * y[i] + C[2]]
    
    x_2D.append(np.dot(x_plane, Pi) - Po[0])
    y_2D.append(np.dot(y_plane, Pi) - Po[1])
    
#    print(str(x_2D[i] / 1000) + ',' + str(y_2D[i] / 1000)) 

# Translate on-plane station positions into 2D (in-plane) coordinate system
# and show station position distribution.

print('Station x,y positions in 2D (in-plane) coordinate system:')

plt.figure()
    
sta_x_2D = []
sta_y_2D = []
for i in range(len(sta)):
    
    Pi = [sta_x[i], sta_y[i], plane_sta_z[i]]
    
    sta_x_2D.append(np.dot(x_plane, Pi) - Po[0])
    sta_y_2D.append(np.dot(y_plane, Pi) - Po[1])
    
    print(sta[i] + ',' + str(sta_x_2D[i]) + ',' + str(sta_y_2D[i]))
    plt.text(sta_x_2D[i], sta_y_2D[i], sta[i])

plt.scatter(x_2D, y_2D, alpha = 0.2 , color = 'black')
plt.scatter(sta_x_2D, sta_y_2D)
print('\n')

# Calculate station distances in-plane and compare to 3D distances

print('Interstation distance disparities between 3D coordinate system and 2D (in-plane) coordinate system:')

for i in range(len(sta)):

    for j in range(len(sta)):
        
        distance_2D = math.sqrt((sta_x_2D[j] - sta_x_2D[i])**2 + (sta_y_2D[j] - sta_y_2D[i])**2)
        
        distance_3D = sta_distances_3D[i][j]
        
        print(sta[i] + ' to ' + sta[j] + ' 3D - 2D distance = ' + str(distance_2D - distance_3D))
        
    print('\n')

# Plot points and fitted surface for comparison

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(X, Y, Z, rstride = 1, cstride = 1, alpha = 0.2, color = 'black')
ax.scatter(data[:,0], data[:,1], data[:,2], c = 'green', s = 10, alpha = 0.1)
ax.scatter(sta_x, sta_y, sta_z, c = 'red', s = 50)
ax.scatter(g_sta_x, g_sta_y, g_sta_z, c = 'purple', s = 50)
ax.scatter(sta_x, sta_y, plane_sta_z, c = 'blue', s = 50)
plt.xlabel('X')
plt.ylabel('Y')
ax.set_zlabel('Z')
ax.axis('equal')
ax.axis('tight')
plt.show()
    