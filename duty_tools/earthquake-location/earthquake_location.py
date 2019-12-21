#!/usr/bin/env python3

"""
Produce earthquake locations
"""

import argparse
import datetime

# Parse arguments

parser = argparse.ArgumentParser()
parser.add_argument('--arrival-time-file', type=str, help='Comma-separated file containing event arrival times details.'
                                                          'First row is header with columns: '
                                                          'Event #, site1_P, site1_S, site2_P,...,siteN_S.'
                                                          'All other rows contain arrival time data as:'
                                                          'event number (e.g. 1),'
                                                          'P wave arrival time at given site in form '
                                                          '"YYYY-MM-DD HH:MM:SS" (UTC),'
                                                          'S wave arrival time at given site in form '
                                                          '"YYYY-MM-DD HH:MM:SS" (UTC).'
                                                          'If no data exist then a date prior '
                                                          'to 1900 should be entered.')
# Note, here the arrival time information is limited by the formatting in Excel. In the finished version of the code
# all time input should be ISO8601 or nan if no data exists.
parser.add_argument('--network-file', type=str, help='Comma-separated file containing network details.'
                                                     'First row is header with columns: '
                                                     'site, latitude, longitude, depth.'
                                                     'All other rows contain site details as: '
                                                     'site code,'
                                                     'latitude (decimal degrees),'
                                                     'longitude (decimal degrees),'
                                                     'depth (metres, +ve direction is down')
parser.add_argument('--velocity-file', type=str, help='Comma-separated file containing velocity model.'
                                                      'First row is header with columns: '
                                                      'upper_depth, lower_depth, velocity.'
                                                      'All other rows contain site details as: '
                                                      'site code,'
                                                      'latitude (decimal degrees),'
                                                      'longitude (decimal degrees),'
                                                      'depth (metres, +ve direction is down')
parser.add_argument('--method', type=str, help='Earthquake location method to use, options are: grid_search,'
                                               'RVT_semblance')
args = parser.parse_args()

# Parse files

arrival_time_data = []
with open(args.arrival_time_file, 'r') as openfile:
    header = -1
    for row in openfile:
        if header == -1:
            header = 0
            arrival_time_data_header = row[:-1].split(',')
        else:
            cols = row[:-1].split(',')
            arrival_time_data.append([])
            arrival_time_data[-1].append(cols[0])
            for col in cols[1:]:
                arrival_time_data[-1].append(datetime.datetime.strptime(col, '%Y-%m-%d %H:%M:%S'))

network_data = []
with open(args.network_file, 'r') as openfile:
    header = -1
    for row in openfile:
        if header == -1:
            header = 0
        else:
            network_data.append(row[:-1].split(','))

velocity_model = []
with open(args.velocity_file, 'r') as openfile:
    header = -1
    for row in openfile:
        if header == -1:
            header = 0
        else:
            velocity_model.append(row[:-1].split(','))

# Build travel time grid to each site

#

print(arrival_time_data)
print(network_data)
print(velocity_model)
