#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2-dimensional grid search using cross-correlation to determine most-likely
epicentre. Designed to use the outputs of spectrogram_event_detection.py.
"""

import concurrent.futures
import datetime
import obspy
import os
import numpy as np
import math
import matplotlib.pyplot as plt

def calculate_distance(x1, y1, x2, y2):
    
    '''
    Calculate distance between two points in a two-dimensional
    cartesian coordinate system.
    '''
    
    distance = math.sqrt((x2 - x1)**2 + (y2 - y1)**2)
    
    return distance




def xcorr_location(event):
    
    '''
    Functionsied version of main body code to do event location.
    Saves xcorr grid to a .npy file.
    '''
    
    print('Processing event ' + event)
    
    # Load events into python as streams
    # Note: this will fail if streams are corrupted
    try:
        
        stream = obspy.read(event_directory + event)
        
    except:
        
        pass
    
    # Parse pre and post trigger times into reference lists
    
    stations = []
    trigger_on = []
    trigger_off = []
    with open(event_directory + event[:-6] + '.csv', 'r') as openfile:
        for row in openfile:
            stations.append(row.split(',')[0])
            trigger_on.append(datetime.datetime.strptime(row.split(',')[1], '%Y-%m-%dT%H:%M:%S.%fZ'))
            trigger_off.append(datetime.datetime.strptime(row.split(',')[2][:-1], '%Y-%m-%dT%H:%M:%S.%fZ'))
    
    # Filter the waveform to remove noise outside of the signal spectral band
    
    stream.filter('bandpass', freqmin = filter_band[0], freqmax = filter_band[-1])
    stream.detrend(type = 'demean')
    stream.detrend(type = 'simple')
    
    start_times = []
    end_times = []
    waveforms = []
    
    # Smooth the data to cross correlate better
    
    for trace in stream:
        
        # Get pre and post trigger times
        
        station = trace.stats.station
        start_times.append(trace.stats.starttime)
        end_times.append(trace.stats.endtime)
        
        pre_trigger_time = (trigger_on[stations.index(station)] - start_times[-1].datetime).total_seconds()
        post_trigger_time = (end_times[-1].datetime - trigger_off[stations.index(station)]).total_seconds()
        
        pre_trigger_samples = int(sampling_rate * pre_trigger_time)
        post_trigger_samples = int(sampling_rate * post_trigger_time) 
#        plt.plot(trace.data, color = 'red')
        # Cut the trace to only the triggered data
    
        trace_data = [0] * (pre_trigger_samples) + \
                      trace.data[pre_trigger_samples : - post_trigger_samples].tolist() + \
                      [0] * (post_trigger_samples)
        
        for i in range(len(trace_data)):
           
           trace_data[i] = abs(trace_data[i])
    
        smoothed_trace_data = trace_data
    
        # Smooth the trace data with a running mean (applied to absolute data)
        
#        plt.plot(trace_data, color = 'blue', alpha = 0.5)
        
#        smoothed_trace_data = [0] * len(trace.data)
#
#        for i in range(pre_trigger_samples, 
#                       len(trace_data) - post_trigger_samples):
#            
#            smoothed_trace_data[i] = sum(trace_data[i - smooth_length : i + smooth_length]) / (2 * smooth_length + 1)
        
        # Normalise the data so cross-correlations aren't bias
#        plt.plot(smoothed_trace_data, color = 'green'); plt.show()
        maxval = max(max(smoothed_trace_data), abs(min(smoothed_trace_data)))
        
        normalised_STD = []
        for k in range(len(smoothed_trace_data)):
                    
            normalised_STD.append(smoothed_trace_data[k] / maxval)
        
        trace_data = normalised_STD
        
        waveforms.append(trace_data)
    
    # Extend the trace_data so their relative sample timing is kept
    
    start_time = min(start_times)
    end_time = max(end_times)

    for w in range(len(waveforms)):
#        plt.plot(waveforms[w], color = 'red')
        waveform_start_time = start_times[w]
        waveform_end_time = end_times[w]

        start_padding = int(sampling_rate * (waveform_start_time - start_time))
        end_padding = int(sampling_rate * (end_time - waveform_end_time))
        
        if start_padding > 0:
            
            for p in range(start_padding):
            
                waveforms[w].insert(0, 0)
            
        if end_padding > 0:
            
            for p in range( end_padding):
            
                waveforms[w].append(0)
#        plt.plot(waveforms[w], color = 'blue', alpha = 0.5); plt.show()
    # Shift the position of each waveform by the negative travel time between
    # the waveform's recording station and a position in the grid, then calculate
    # the cross-correlation coefficient between the shifted waveforms.
    
    xcorr_value_grid = np.zeros(gridx.shape)
    
    for i in range(len(gridx)):
        for j in range(len(gridx[i])):
            
            shifted_waveforms = []
            
            # Shift waveforms
            
            for w in range(len(waveforms)):
                
                travel_time = [travel_times[j][i][stations.index(stream[w].stats.station)]]
                
                shifted_waveform = [0] * len(waveforms[w])
#                plt.plot(waveforms[w], color = 'red')                
                for k in range(len(waveforms[w])):
                    
                    # This will fail if the shift goes outside the original data window
                    # Fix it by expanding the original data window!

                    try:
                        
                        shifted_waveform[k] = waveforms[w][k - int(round(travel_time[0] * sampling_rate))]
                        
                    except:
                        
                        pass
#                plt.plot(shifted_waveform, color = 'blue', alpha = 0.5); plt.show()
                shifted_waveforms.append(np.asarray(shifted_waveform))
                
            # Calculate cross-correlation coefficients
            
            shifted_waveforms = np.asarray(shifted_waveforms)
            
            correlation_coefficients = np.corrcoef(shifted_waveforms)
            
            # Calculate mean cross-correlation coefficient
            
            mean_cc = 0
            num_cc = 0
            for row in correlation_coefficients:
                for col in row:
                    if float(col) != 1:
                        mean_cc += col
                        num_cc += 1
                        
            mean_cc /= num_cc
            
            xcorr_value_grid[i][j] = mean_cc
            
    np.save(event_directory + event[:-6] + '.xcorrvaluegrid.npy', xcorr_value_grid)




# Define station names and positions (XYZ)

stations = ['TSNC1', 'TSNC3', 'TSNL2', 'TSNL3', 'TSNR2', 'TSNR3']

station_positions = [[1374.062, 5164.009, 0.821],
                     [1374.498, 5166.080, 0.904],
                     [1374.261, 5163.408, 0.816],
                     [1374.783, 5165.049, 0.857],
                     [1373.533, 5163.353, 0.801],
                     [1373.689, 5165.150, 0.857]]

# Set grid size and scale (km)

xmin = 1373
xmax = 1375
ymin = 5163
ymax = 5166.5
xstep = 0.05
ystep = 0.05

# Set grid velocity (km / s)

velocity = 1.65

sampling_rate = 250

event_directory = '/home/sam/EVENTS_IT3/TYPE_A/4/'

# Set smoothing length (in samples) to use in running mean smoothing of waveforms

smooth_length = 25

## Set frequency band for bandpass filtering (generally same as detection band)

filter_band = [1, 25]

# Parse event file paths

files = os.listdir(event_directory)
events = []
for event in files:
    if (event[-6:] == '.MSEED') and (event[:-6] + '.xcorrvaluegrid.npy' not in files):
        events.append(event)
    
# Make grid
    
gridx, gridy = np.meshgrid(np.linspace(xmin, xmax, (xmax - xmin + 1) / xstep), np.linspace(ymin, ymax, (ymax - ymin + 1) / ystep))

# Generate travel times from each station to each grid point

travel_times = [[[] for i in range(len(gridx))] for j in range(len(gridx[0]))]

for i in range(len(gridx)):
    for j in range(len(gridx[i])):
        
        x = gridx[i][j]
        y = gridy[i][j]
        
        for station_position in station_positions:
            
            distance = calculate_distance(x, y, station_position[0], station_position[1])
            
            travel_times[j][i].append(distance / velocity)
            
        # Make travel times relative
        
        relative_tt = [[] for k in range(len(travel_times[j][i]))]
        for k in range(len(travel_times[j][i])):
            
            relative_tt[k] = travel_times[j][i][k] - travel_times[j][i][0]
            
        travel_times[j][i] = relative_tt
            
    
# Locate events
    
with concurrent.futures.ProcessPoolExecutor() as executor:
        
    executor.map(xcorr_location, events)
        
#for event in events:
    
#    xcorr_location(event)
                    

        
        
        
    