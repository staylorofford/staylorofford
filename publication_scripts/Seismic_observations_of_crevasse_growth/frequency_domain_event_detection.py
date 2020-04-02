#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Frequency domain seismic detection algorithm,
uses spectrum files saved from spectrum_generation.
DEPRECIATED!!! USE spectrogram_event_detection.py instead!
"""

# Import packages

import datetime
import glob
import numpy as np
import obspy
import matplotlib.pyplot as plt

# Set parameters

## Directory to save event files to

event_output_directory = '/media/sam/61D05F6577F6DB39/SCIENCE/FDD/'

## Directory to load spectrum files from

spectrum_directory = '/home/samto/spectrums/'

## Seismic component to use in processing

stream_component = 'Z'

## Stations to process

stream_stations = ['TSNC1', 'TSNC3', 'TSNL2', 'TSNL3', 'TSNR2', 'TSNR3']

## Start and end dates for processing of day-long seismic streams

start_year = '2016'
start_month = '04'
start_day = '21'

end_year = '2017'
end_month = '08'
end_day = '01'

## Set signal and noise frequency bands

signal_bands = [[1, 20]]
noise_bands = [[21, 125]]

## Set weighting for each band

signal_band_weights = [1]
noise_band_weights = [1]

## Set FFT window length (seconds)

FFT_window_length = 1

## Set FFT window overlap (fractional percentage)

FFT_window_overlap = 0.9

## Set threshold required to trigger detection consideration

trigger_threshold = 10

## Set minimum trigger length for a trigger to be registered

minimum_trigger_length = 0.2

## Set pre and post event time (s) for event files

pre_event_time = 10
post_event_time = 10

## For network detection: set maximum interstation detection delay time (seconds)

delay_time = 60

## For network detection: set minimum time between successive events

interevent_time = 0

## For network detection: set number of detections within successive delay times
## required for a detection (default = 1, only enables network detection if 
## station_threshold > 1)

station_threshold = 3

# Convert start and end dates into datetime objects, and get them as julian days in their respective years
    
start_date = datetime.datetime.strptime(start_year + '-' + start_month + '-'+ start_day, '%Y-%m-%d')
end_date = datetime.datetime.strptime(end_year + '-' + end_month + '-'+ end_day, '%Y-%m-%d')

start_date_doy = start_date.timetuple().tm_yday
end_date_doy = end_date.timetuple().tm_yday

# Look through data in all streams within the processing window

years = range(int(start_year), int(start_year) + int(end_year) - int(start_year))

times = []

for year in years:
    
    for doy in range(366):
    
        if (year == int(start_year)) and (doy < int(start_date_doy)): continue                        
        elif (year == int(end_year)) and (doy > int(end_date_doy)): continue            
        else:
            
            spectrum_files = glob.glob(spectrum_directory + '*' + str(year) + '*' + str(doy) + '*spectrums.npy')
                        
            # Apply component and station filtering
            
            for spectrum_file in spectrum_files:
                
                stream_file_metadata = spectrum_file.split('/')[-1].split('.')
                component = stream_file_metadata[3][-1]
                station = stream_file_metadata[0]
                
                if component != stream_component: continue
                if station not in stream_stations: continue
                
                # Load the spectrum file
            
                spectrum_list = np.load(spectrum_file).tolist()
                
                # Join spectra lists together
                
                times = []
                spectrums = []
                for i in range(len(spectrum_list)):
                    for j in range(len(spectrum_list[i])):
                        if len(spectrum_list[i][j]) > 0:
                            times.append([])
                            spectrums.append([])
                            times[-1] = (spectrum_list[i][j][0])
                            # Miss the first frequency entry as it is not representative of
                            # particular harmonics.
                            spectrums[-1] = np.absolute(spectrum_list[i][j][1][1:].real)
                
                times = np.array(times)
                spectrums = np.array(spectrums)
                
                # Calculate triggers
                
                print('Generating triggers')
                
                trigger_on = False
                triggers = []
                
                for i in range(len(spectrums)):
                    
                    midtime = times[i]
                    spectrum = spectrums[i]
                
                    signal = 0
                    noise = 0      
                    
                    # Calculate band ratios (normalise by length of band)
                    
                    for signal_band in signal_bands:
                        
                        signal += sum(abs(spectrum[signal_band[0] : signal_band[1]])) * signal_band_weights[signal_bands.index(signal_band)] \
                                    / (signal_band[1] - signal_band[0])
                        
                    for noise_band in noise_bands:
                        
                        noise += sum(abs(spectrum[noise_band[0] : noise_band[1]])) * noise_band_weights[noise_bands.index(noise_band)] \
                                    / (noise_band[1] - noise_band[0])
                    
                    try:
                        
                        band_ratio = signal / float(noise)
                        
                    except:
                        
                        continue
                    
                    # Determine and save triggers
                    
                    if (trigger_on == False) and (band_ratio >= trigger_threshold):
                        
                        trigger_on = True    
                        trigger_on_index = i
                        
                    elif (trigger_on == True) and (band_ratio < trigger_threshold):
        
                        trigger_on = False
                        trigger_off_index = i
                        
                        if times[trigger_off_index] - times[trigger_on_index] >= minimum_trigger_length:
                        
                            triggers.append([times[trigger_on_index],
                                             times[trigger_off_index],
                                             station])
                                             
                            print(triggers[-1])
                        
            # Sort all the day's triggers chronologically                    
            
            triggers.sort()
            
            # Find coincident triggers and write event files
            
            print('Locating coincident triggers')
        
            last_off = 0 # Set the end time of the previous trigger
            coincidence_triggers = []
            
            count = -1
            while count < (len(triggers) - 1):
                
                count += 1
                
                event = False

                # Get the first trigger
                 
                on, off, station = triggers[count]

                # Prepare the event output
                
                event_on_time = on
                event_stations = [station]
                coincidence_sum = 1
                              
                # Avoid overlapping events
                              
                if on < last_off + interevent_time: continue
                
                # Compare current trigger with temporally-nearby triggers
                    
                for later_trigger in triggers:
                    
                    lt_on, lt_off, lt_station = later_trigger
                    
                    # Check triggers are appropriately temporally clustered
                    
                    if lt_on > off + delay_time: break
                    
                    # Only allow one trigger per station
        
                    if lt_station in event_stations: continue
                                         
                    # Update event data for new trigger
                    
                    event_stations.append(lt_station)
                    coincidence_sum += 1
                    event_off_time = max([off, lt_off])        
                    event_length = event_off_time - event_on_time                
                    last_off = event_off_time
                    
                    # If the coincidence sum is high enough, save event data
                    
                    if coincidence_sum >= station_threshold:
                        
                        event = True
                        
                if event == True:
                    
                        coincidence_triggers.append(event_on_time,
                                event_off_time,
                                event_stations,
                                event_length)
                                
                        print(coincidence_triggers[-1])
                        
            # Save event files
                        
            print('Saving event files')
        
            for coincidence_trigger in coincidence_triggers:
                
                print('Saving ' + str(coincidence_trigger))

                # Build seismic stream
                # NOTE: contains all components and stations available in stream files
                
                stream = obspy.Stream()
                
                for stream_file in stream_files:
                    
                    stream += stream_file.trim(starttime = coincidence_trigger[0] - pre_event_time,
                                               endtime = coincidence_trigger[1] + post_event_time)
 
#                stream.write(event_output_directory + coincidence_trigger[0] + '.MSEED',
#                             format = 'MSEED', cencoding = 'STEIM2') 