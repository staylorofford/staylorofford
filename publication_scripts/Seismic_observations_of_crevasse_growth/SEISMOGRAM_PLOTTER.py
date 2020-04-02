#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May  9 13:24:50 2017

@author: taylorsa
"""

import obspy
import os
import matplotlib.pyplot as plt

# look at all seismograms in a folder


opendir = '/home/samto/PERSONAL_SCIENCE/EVENTS/TYPE_A/'
openfiles=os.listdir(opendir)
    
# only allow .MSEED files to be in the event list    
    
i=0
for e in range(len(openfiles)):
    try:
        if openfiles[i][-6:]=='.MSEED':
            i+=1
        else:
            openfiles.pop(openfiles.index(openfiles[i]))
    except:
        openfiles.pop(openfiles.index(openfiles[i]))

openfiles.sort()

last_event=len(openfiles)
valid_number=0
while valid_number==0:
    print('What event would you like to begin at?')
    start_event = input('Event number: ')
    try:
        int(start_event)
        valid_number=1
        
        # reset if outside bounds
        
        if int(start_event)<0:
            valid_number=0              
        elif int(start_event)>last_event:
            valid_number=0
    except:
        print('That is not a valid number. Try again.')

for o in range(len(openfiles)):
    
    # load in the file
    
    openfile=openfiles[o+int(start_event)]
    
#    if openfiles not in [1, 2]: continue # ONLY ALLOW EVENTS 2,3
    
    st=obspy.read(opendir+openfile)
    stog=obspy.Stream()

    for tr in st:
#         if tr.stats.station not in ['TSNM1','TSNM2','TSNM3','TSNC4']:
        if tr.stats.station=='TSNC1':
             stog.append(tr)
              
     # z component and normalise
             
    st=stog   
    st.detrend(type='demean')
    st=st.select(channel='*Z') # ONLY USE Z COMP
    
#    # trim to data
# 
#    if o+int(start_event)==1:
#        st.trim(starttime=obspy.UTCDateTime("2016-04-21T06:36:01.1"), endtime=obspy.UTCDateTime("2016-04-21T06:36:02.81"))
#    elif o+int(start_event)==2:
#        st.trim(starttime=obspy.UTCDateTime("2016-05-15T06:43:38.00"))#, endtime=obspy.UTCDateTime("2016-05-15T06:43:45.8"))

#    st.trim(starttime=obspy.UTCDateTime("2016-05-15T06:43:43.00"), endtime=obspy.UTCDateTime("2016-05-15T06:43:44.5"))
    
    st.normalize()
    fig=plt.figure(figsize=[2.75591, 2.75591])
    Zfig=plt.plot(fig=fig, show=False)
    plt.close()
    st.plot(fig=Zfig,show=False,type="relative")
    
    ax=plt.gca()
#    ax.set_xticks(np.arange(0,len(st[0])/250.0,step[o+int(start_event)]))
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.tick_params(which = 'both', direction = 'out')
    ax.tick_params(which = 'both', length = 3)
    ax.tick_params(which = 'both', width = 1.1)
    ax.tick_params(axis='both',which='major',labelsize=10)
    ax.set_xlabel('Time (s)', fontsize=10, labelpad=5)
#    ax.set_ylabel('Normalised amplitude', fontsize=10, labelpad=5) # cant set this position well enough
    ax.get_yaxis().set_label_coords(-0.10,0.5)
    ax.set_title("", fontsize=0)
    ax.legend()
    
#    mngr = plt.get_current_fig_manager()
#    mngr.window.setGeometry(0,0,1920,1000) # take up the left hand screen
    
    x=-1
    y=0
    for trace in st:
        x+=1
        if x>5:
            x-=6
            y+=1
        print('Next spectrogram is for station '+str(trace.stats.station))
        trace.spectrogram(show=False, wlen=0.1, per_lap=0.95, cmap='jet') # optimised parameters for clarity and speed
        
        ax=plt.gca()
#        ax.set_xticks(np.arange(0,len(st[0])/250.0,step[o+int(start_event)]))
        ax.spines['top'].set_visible(True)
        ax.spines['right'].set_visible(True)
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
        ax.tick_params(which = 'both', direction = 'out')
        ax.tick_params(which = 'both', length = 3)
        ax.tick_params(which = 'both', width = 1.1)   
        ax.tick_params(axis='both',which='major',labelsize=10)
        ax.set_xlabel('Time (s)', fontsize=10, labelpad=5)
        ax.set_ylabel('Frequency (Hz)', fontsize=10, labelpad=5)
#        ax.get_yaxis().set_label_coords(-0.15,0.5)
        ax.set_title("")
        
#        mngr = plt.get_current_fig_manager()
#        mngr.window.setGeometry(1920+x*640,0+y*540,640,540) # tesselate these windows across the two remaining screens
        
    plt.show(block=False)
    
    # query the user to move on / print time of event
    
    print(openfile)
    print(o+int(start_event))
    print('Template?')
    
    ans=input('[Y/n]: ')
    if ans=='Y':
        for tr in st:
            print('Potential template at time '+str(tr.stats.starttime))
            break
    plt.close('all')