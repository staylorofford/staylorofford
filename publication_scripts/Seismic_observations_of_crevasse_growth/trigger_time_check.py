#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot trigger times over seismograms to check spectrogram detection functionality.
"""

import glob
import obspy
import matplotlib.pyplot as plt
import datetime

event_directory = '/home/sam/EVENTS_IT3/TYPE_A/4/'

seismic_event_files = glob.glob(event_directory + '*MSEED')
trigger_time_files = glob.glob(event_directory + '*csv')

seismic_event_files.sort()
trigger_time_files.sort()

for e in range(len(seismic_event_files)):
    
    event = obspy.read(seismic_event_files[e])
    
    with open(seismic_event_files[e][:-6] + '.csv', 'r') as openfile:
        
        rc = 0  
        for row in openfile:
            rc += 1
    
    if rc < 3:
        print('Less than 3 re-triggers for event ' + seismic_event_files[e])
        continue

    mintime = obspy.UTCDateTime('9999-01-01T00:00:00').datetime
    
    for trace in event:
        
        mintime = min(trace.stats.starttime.datetime, mintime)
    
    c = 0
    for trace in event:
        c += 1  
        
        plotted = False
        
        data = trace.data
        start = trace.stats.starttime.datetime
        timediff = (start - mintime).total_seconds()
#        print(start, mintime, timediff)
        station = trace.stats.station

        with open(seismic_event_files[e][:-6] + '.csv', 'r') as openfile:
            
            for row in openfile:
    
                if row.split(',')[0] == station:
                    sample = 250 * (obspy.UTCDateTime(row.split(',')[1], iso8601=True).datetime - start).total_seconds() + timediff
                    off_sample = 250 * (obspy.UTCDateTime(row.split(',')[2], iso8601=True).datetime - start).total_seconds() + timediff
                    plt.subplot(int(str(len(event)) + '1' + str(c)))
                    plt.plot(data, color = 'black')
                    plt.axvline(sample, color = 'red')
                    plt.axvline(off_sample, color = 'red')
                    plotted = True
                    
                    break
        
    if plotted == True:
        
        plt.show(block = False)
        event.plot()
        plt.close('all')
        
        