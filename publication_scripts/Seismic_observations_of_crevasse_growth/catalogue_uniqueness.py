#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Separate events from catalogues so that all events are only in one catalogue.
Works by events having second-decimated start times.
NOTE: this version assumes the lowest directory is the threshold value dir.
"""

def forward_difference(data):
    
    '''
    Calculate forward difference
    '''
    
    diff = []
    for i in range(1, len(data)):
        diff.append(float(data[i]) - float(data[i - 1]))
    return diff




import os
import datetime
import matplotlib.pyplot as plt

#catalogue_root_dir = '/home/samto/EVENTS_IT1/'
#catalogues = ['TYPE_A', 'TYPE_B', 'TYPE_C', 'TYPE_D', 'TYPE_E', 'TYPE_F']
#threshold = '5'
#catalogues_to_purify = ['TYPE_A', 'TYPE_B', 'TYPE_C', 'TYPE_D', 'TYPE_E', 'TYPE_F']

catalogue_root_dir = '/home/sam/EVENTS_IT3/'
catalogues = ['TYPE_A', 'TYPE_B', 'TYPE_C', 'TYPE_D']
#catalogues = ['TYPE_A']#,'TYPE_C','TYPE_D']
thresholds = ['4', '4', '3', '3']
catalogues_to_purify = catalogues

import numpy as np
import scipy.optimize as linefit

def func(x,m,c):
    # function for line fit optimisation
    return m*x+c

# Get all events in the catalogues

events = [[] for i in range(len(catalogues))]
for i in range(len(catalogues)):
    
    events[i] = os.listdir(catalogue_root_dir + catalogues[i] + '/' + thresholds[i] + '/')
    
for event_series in events:
    
    cumulative_event_sum = 0
    cumulative_event_sums = []
    event_times = []
    
    # Generate cumulative sum list
    # and list of event times
    
    sorted_PE = sorted(event_series)
    
    for event in sorted_PE:
        
        if event[-6:] != '.MSEED': continue
        
        cumulative_event_sum += 1
        cumulative_event_sums.append(cumulative_event_sum)
    
        event_time = datetime.datetime.strptime(event[:event.index('.')], '%Y-%m-%dT%H:%M:%S')
        event_times.append(event_time)
        
    plt.plot(event_times, cumulative_event_sums)
    plt.xlim(datetime.datetime.strptime('2016-05-01', '%Y-%m-%d'),
             datetime.datetime.strptime('2016-05-19', '%Y-%m-%d'))
    plt.show()
    
#    event_reltimes = []
#    for i in range(len(event_times)):
#        event_reltimes.append((event_times[i] - event_times[0]).total_seconds())
#        
#    x=10
#    slopes=[]
#    for i in range(len(event_times)-x):
#        line=linefit.curve_fit(func, event_reltimes[i:i+x], cumulative_event_sums[i:i+x])
#        slopes.append(line[0][0])
#    
#    for slope in slopes:
#        slopes[slopes.index(slope)] /= max(slopes)
#    
#    plt.plot(event_times[:-x], slopes)
#    plt.xlim(datetime.datetime.strptime('2016-05-01', '%Y-%m-%d'),
#             datetime.datetime.strptime('2016-05-19', '%Y-%m-%d'))
#plt.show()
    
#
## Remove all shared events from some catalogues
#    
#purified_events = [[] for i in range(len(catalogues))]
#refuse_events = [[] for i in range(len(catalogues))]
#events_copy = events
#if len(catalogues) > 1:
#    for i in range(len(catalogues)):
#        for j in range(len(catalogues_to_purify)):
#            k = catalogues.index(catalogues_to_purify[j])
#            if k != i:
#                for event in events[i]:
#                    for test_event in events[k]:
#                        if event == test_event:
#                            refuse_events[i].append(event)
#                            break
#                    purified_events[i].append(event)
#else:
#    purified_events = events
#    
## Plot it
#            
#for j in range(len(purified_events)):
#                
#    cumulative_event_sum = 0
#    cumulative_event_sums = []
#    event_times = []
#    
#    # Generate cumulative sum list
#    # and list of event times
#    
#    sorted_PE = sorted(purified_events[j])
#    
#    for event in sorted_PE:
#        
#        cumulative_event_sum += 1
#        cumulative_event_sums.append(cumulative_event_sum)
#    
#        event_time = datetime.datetime.strptime(event[:event.index('.')], '%Y-%m-%dT%H:%M:%S')
#        event_times.append(event_time)
#
#    # Plot cumulative event sum over time
#        
#    plt.plot(event_times, cumulative_event_sums)
#    plt.show()
#
#for j in range(len(refuse_events)):
#                
#    refuse_cumulative_event_sum = 0
#    refuse_cumulative_event_sums = []
#    refuse_event_times = []
#    
#    # Generate cumulative sum list
#    # and list of event times
#    
#    sorted_RE = sorted(refuse_events[j])
#    
#    for event in sorted_RE:
#        
#        refuse_cumulative_event_sum += 1
#        refuse_cumulative_event_sums.append(refuse_cumulative_event_sum)
#    
#        event_time = datetime.datetime.strptime(event[:event.index('.')], '%Y-%m-%dT%H:%M:%S')
#        refuse_event_times.append(event_time)
#        
#    # Plot cumulative event sum over time
#        
#    plt.plot(refuse_event_times, refuse_cumulative_event_sums)
#    plt.show()