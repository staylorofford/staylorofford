"""
Calculate, for a given network configuration and a given detection system, which events (earthquakes) would be
detected independent of event size or attenuation.
"""

import math
import numpy as np
from obspy.taup import TauPyModel
import pandas as pd
import matplotlib.pyplot as plt


# Define velocity model

spherical_velocity_model = TauPyModel(model="iasp91")

# Define the grid to calculate detection capability over

minlatitude, maxlatitude = -49, -33  # minimum and maximum latitude for event
minlongitude, maxlongitude = 163, -175  # western and eastern longitude for events
mindepth, maxdepth = 0, 40  # minimum and maximum depth (+ve down) for events (integer only)
dlat, dlon, dh = 0.1, 0.1, 10  # grid spacing values: difference in latitude, longitude, and height
# between grid points in each respective direction. Note: methods assumes grid spacing
# is sufficiently small to capture all variations in the detectability. As grid is
# defined in spherical coordinates, the operator is reminded to consider the relationship
# between 1 degree of longitude and distance at varying degrees of latitude, i.e. that
# on regional or national scales the true spacing of grid points will be irregular.

# Give the path to the station file: a csv containing in each row the site code, latitude, longitude, and elevation
# (+ve up, in metres) of each station to use in the detectability calculations. The file should contain a header
# denoting these column titles in the form of "Station", "Latitude", "Longitude", "Elevation" in the corresponding
# positions of the data with each row containing only this data for one station. The code will find these columns and
# only parse the data in these columns. For example, a header with columns:
# Station,Network,Name,Latitude,Longitude,Elevation,Datum,Start Date,End Date
# will lead to only data in the 1st, 4th, 5th, and 6th columns being parsed.
station_file = '/home/samto/PROCESSING/station_file.csv'

# Define the detection system parameters

detections_for_event = 10  # number of detections required to define an event
detection_grouping_time = 15  # seconds between any detections for them to be grouped at a grid point

# Build the grid

gridx = list(range(minlongitude, maxlongitude + dlon, dlon))
gridy = list(range(minlatitude, maxlatitude + dlat, dlat))
gridz = list(range(mindepth, maxdepth + dh, dh))

# Parse the station metadata from file

station_metadata = pd.read_csv(station_file, header=0, index_col=0)

# Calculate distances and travel times from each grid point to each station

p_wave_travel_times = [[[([0] * len(station_metadata))
                         for k in range(len(gridz))]
                        for j in range(len(gridy))]
                       for i in range(len(gridx))]
for i in range(len(gridx)):
    for j in range(len(gridy)):
        for k in range(len(gridz)):
            for m in range(len(station_metadata)):
                hypocentre = [gridx[i], gridy[j], gridz[j]]
                station_location = [station_metadata.iloc[m]['Latitude'],
                                    station_metadata.iloc[m]['Longitude'],
                                    -1/1000.0 * station_metadata.iloc[m]['Elevation']]
                # Calculate epicentral distance in degrees
                delta = math.degrees(2 * (math.asin(((math.sin(1 / 2 * math.radians(abs(hypocentre[0] -
                                                                                        station_location[0])))) ** 2 +
                                                     math.cos(math.radians(hypocentre[0])) *
                                                     math.cos(math.radians(station_location[0])) *
                                                     (math.sin(1 / 2 * math.radians(abs(hypocentre[1] -
                                                                                        station_location[1])))) ** 2) **
                                                    (1 / 2))))
                # Calculate arrival times of P wave
                arrivals = spherical_velocity_model.get_travel_times(source_depth_in_km=hypocentre[2],
                                                                     receiver_depth_in_km=min(0, station_location[2]),
                                                                     distance_in_degree=delta,
                                                                     phase_list=['p', 'P'])
                p_tt = None
                for arrival in arrivals:
                    if arrival.name == 'p' or arrival.name == 'P':
                        p_tt = arrival.time
                        break
                p_wave_travel_times[i][j][k][m] = p_tt

# For each grid point, see if a detection would be made.
# This follows the assumption that if any N picks occur within T time for a given hypocentre a detection would be made.
# So the detection system assumed to be in use is one that takes all detections occuring in a similar time period and
# backprojects them to grid points. If N backprojected detections occur at a grid point within T time, then it considers
# that a verification that an earthquake produces those detections and hence produces an event detection.
grid_point_detection_values = [[([0] * len(gridz))
                                for j in range(len(gridy))]
                               for i in range(len(gridx))]
for i in range(len(p_wave_travel_times)):
    for j in range(len(p_wave_travel_times[i])):
        for k in range(len(p_wave_travel_times[i][j])):
            # We have the set M arrival times, where M is the number of stations in the network
            grid_point_p_wave_travel_times = p_wave_travel_times[i][j][k]
            # Sort the travel times
            grid_point_p_wave_travel_times.sort()
            # Run through the sorted list and see if any N travel times occur within T time
            for m in range(len(grid_point_p_wave_travel_times) - detections_for_event):
                dt = (grid_point_p_wave_travel_times[m + detections_for_event] -
                      grid_point_p_wave_travel_times[m]).total_seconds()
                if dt <= detection_grouping_time:
                    # A detection will occur, set the detection value to True (1) and move on to the next grid point
                    grid_point_detection_values[i][j][k] = 1
                    break

# Plot the data in contour plots descending through the grid

for k in range(len(gridz)):
    # Extract the data from the grid at this depth
    depth_slice = [([False] * len(gridy))
                   for i in range(len(gridx))]
    for i in range(len(gridx)):
        for j in range(len(gridy)):
            depth_slice[i][j] = grid_point_detection_values[i][j][k]
    # Plot the data in a countour map
    plt.contour(gridx, gridy, depth_slice, levels=[0, 1])
    plt.savefig('network_detectability_' + str(gridz[k]) + 'km.png')
