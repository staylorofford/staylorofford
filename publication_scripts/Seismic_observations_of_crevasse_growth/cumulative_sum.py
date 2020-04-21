#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate cumulative sum plot for events
"""

import os
import datetime
import matplotlib.pyplot as plt

# Get list of events

event_dir = '/home/samto/PERSONAL_SCIENCE/EVENTS/TYPE_D/'
events = os.listdir(event_dir)
events.sort()

cumulative_event_sum = 0
cumulative_event_sums = []
event_times = []

# Generate cumulative sum list
# and list of event times

for event in events:
    
    cumulative_event_sum += 1
    cumulative_event_sums.append(cumulative_event_sum)

    event_time = datetime.datetime.strptime(event[:event.index('.')], '%Y-%m-%dT%H:%M:%S')
    event_times.append(event_time)
    
# Plot cumulative event sum over time
    
plt.plot(event_times, cumulative_event_sums)
plt.show()