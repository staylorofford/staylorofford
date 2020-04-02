#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Frequency domain seismic detection algorithm
using temporal as well as frequency data,
uses spectrum files saved from spectrum_generation.
"""

# Import packages

import datetime
import glob
import numpy as np
import obspy
import matplotlib.pyplot as plt

# Set parameters

## Directory to save event files to

#event_output_directory = '/media/sam/61D05F6577F6DB39/SCIENCE/FDD/'
#event_output_directory = '/home/samto/PERSONAL_SCIENCE/EVENTS/TYPE_D/'
event_output_directory = '/home/sam/EVENTS_IT3/TYPE_A/5/'

## Directory to load spectrum files from (all spectrum files in one folder)

spectrum_directory = '/media/sam/61D05F6577F6DB39/SCIENCE/Spectrums/'

## Stream root directory contains individual day-long streams of each station's components
## within julian day directories (same data as spectrums were generated from!)

stream_root_directory = '/media/sam/61D05F6577F6DB39/SCIENCE/day_volumes_S/'
#stream_root_directory = '/SEISMIC_DATA/day_volumes_S/'
#stream_root_directory = '/home/samto/day_volumes_S/'

## Seismic component to use in processing

stream_component = 'Z'

## Stations to process

stream_stations = ['TSNC1', 'TSNC3', 'TSNL2', 'TSNL3', 'TSNR2', 'TSNR3']#, 'TSNM1', 'TSNM3']

## Start and end dates for processing of day-long seismic streams

start_year = '2016'
start_month = '04'
start_day = '21'

end_year = '2017'
end_month = '08'
end_day = '01'

## Set signal and noise frequency bands

signal_bands = [[1, 20]]
noise_bands = signal_bands

## Set weighting for each band

signal_band_weights = [1]
noise_band_weights = [1]

## Set threshold required to trigger detection consideratio

trigger_threshold = 5
trigger_overload = 999999

## Set trigger length range for a trigger to be register
## NOTE: this in in FFT windows

trigger_min_len = 0
trigger_max_len = 2

## Set minimum and maximum event durations

event_min_len = 0
event_max_len = 60

## Set pre and post event time (s) for event files

#pre_event_time = 5
#post_event_time = 5
pre_post_time = 10

## For network detection: set maximum interstation detection delay time in FFT windows

delay_time = 2

## For network detection: set minimum time between successive events (in FFT windows)

interevent_time = 2

## For network detection: set number of detections within successive delay times
## required for a detection (default = 1, only enables network detection if 
## station_threshold > 1)

station_threshold = 3

## Enter FFT window length (seconds) for mode II event file writing

FFT_window_len = 1

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
            stream_files = glob.glob(stream_root_directory + 'Y' + str(year) + '/R' + str(doy) + '.01/*')
                        
            # Apply component and station filtering

            triggers = []
            band_ratios = []
            trigger_durations = []
            trigger_SNRs = []
            event_durations = []
            event_trigger_SNRs = []
            all_event_triggers = []
            
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
                            
                # Ensure times are chronological
                # Unsorting may result from multithreading during spectrum generation
                times, spectrums = zip(*sorted(zip(times, spectrums)))
                
                times = np.array(times)
                spectrums = np.array(spectrums)
                
                # Calculate trigger section:
                
                print('Generating triggers for station ' + station + ' on day ' + str(doy) + ' in ' + str(year))
                
                trigger_on = False
                
                # Generate first noise value
                
                noise = 0
                
                for noise_band in noise_bands:
                        
                    noise += sum(abs(spectrums[0][noise_band[0] : noise_band[1]])) * noise_band_weights[noise_bands.index(noise_band)] \
                            / (noise_band[1] - noise_band[0])
                            
                # Compare signal values to noise values and generate triggers
                # Each signal value is compared to the previous "signal" value
                # (noise value) until the signal/noise threshold is passed,
                # at this point each successive signal value is compared to
                # the noise value at the time of triggering. This allows long
                # triggers to occur which is useful for separating different
                # spectral signatures.
                
                for i in range(1, len(spectrums)):
                    
                    midtime = times[i]
                    spectrum = spectrums[i]
                
                    signal = 0
                    
                    # Calculate band ratios (normalise by length of band)
                    
                    for signal_band in signal_bands:
                        
                        signal += sum(abs(spectrum[signal_band[0] : signal_band[1]])) * signal_band_weights[signal_bands.index(signal_band)] \
                                    / (signal_band[1] - signal_band[0])
                    
                    try:
                        
                        band_ratio = signal / noise
                        band_ratios.append(100 * band_ratio)
                        
                    except:
                        
                        print('Signal:noise calculation failed at index ' + str(i) + ' !')
                        continue
                    
                    # Determine and save triggers
                    
                    if (trigger_on == False) and (trigger_overload > band_ratio >= trigger_threshold):
                        
                        trigger_on = True    
                        trigger_on_index = i
                        
                    elif (trigger_on == True) and (band_ratio < trigger_threshold):
        
                        trigger_on = False
                        noise = signal
                        trigger_off_index = i
                        trigger_durations.append(trigger_off_index - trigger_on_index)
#                        print(trigger_on_index, trigger_off_index)
#                        print(times[trigger_on_index], times[trigger_off_index])
                        triggers.append([times[trigger_on_index],
                                         times[trigger_off_index],
                                         station])
                                         
#                        print(triggers[-1])        
                        
                    elif (trigger_on == False):
                        
                        noise = signal
                        
                    if band_ratio >= trigger_threshold:
                        
                        trigger_SNRs.append(100 * band_ratio)
                        
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
                trigger_SNR = trigger_SNRs[count]
                
                # Avoid overlapping events
                              
                if on < last_off + interevent_time: continue
            
                # Check trigger duration is desired
                
                if not (trigger_min_len <= (off - on) <= trigger_max_len): continue
            
                # Prepare the event output
                
                event_on_time = on
                event_stations = [station]
                event_SNRs = [trigger_SNR]
                event_triggers = [triggers[count]]
                coincidence_sum = 1
                
                # Compare current trigger with temporally-nearby triggers
                
                for later_trigger in triggers[count + 1:]:
                    
                    lt_on, lt_off, lt_station = later_trigger
                    lt_trigger_SNR = trigger_SNRs[triggers.index(later_trigger)]

                    # Check trigger duration is desired

                    if not (trigger_min_len <= (lt_off - lt_on) <= trigger_max_len): continue

                    # Check triggers are appropriately temporally clustered
                    
                    if lt_on > off + delay_time: break

                    # Only allow one trigger per station
        
                    if lt_station in event_stations: continue

                    # Update event data for new trigger
                    
                    event_stations.append(lt_station)
                    event_triggers.append(later_trigger)
                    event_SNRs.append(lt_trigger_SNR)
                    coincidence_sum += 1
                    event_off_time = max([off, lt_off])        
                    event_length = event_off_time - event_on_time                
                    last_off = event_off_time
                    
                    # If the coincidence sum is high enough, save event data
                    
                    if coincidence_sum >= station_threshold:

                        event = True
                        
                if event == True:
                    
                        if (event_min_len < event_length < event_max_len):
                    
                            coincidence_triggers.append([event_on_time,
                                    event_off_time,
                                    event_stations,
                                    event_length])
        
                            event_durations.append(event_length)
                            event_trigger_SNRs.extend(event_SNRs)
                            
                            # Save trigger details for alternative uses
                            # such as trimming streams before xcorr
                            
                            all_event_triggers.append(event_triggers)
                                    
                            print(coincidence_triggers[-1])
                        
            # Save event files:
            
#            plt.figure()
#            plt.subplot(221)
#            plt.hist(band_ratios, range(0, 1000))   
#            plt.title('Band ratio distribution', fontsize = 10)
#            plt.subplot(222)
#            plt.hist(trigger_durations, range(0, 20))
#            plt.title('Trigger duration distribution', fontsize = 10)
#            plt.subplot(223)
#            plt.hist(event_trigger_SNRs, range(0, 1000))
#            plt.title('Event trigger band ratio distribution', fontsize = 10)
#            plt.subplot(224)
#            plt.hist(event_durations, range(0, 20))
#            plt.title('Event duration distribution', fontsize = 10)
#            plt.show()
            
            # Build seismic stream
            # NOTE: contains all components and stations available in stream files

            print('Saving event files')

            # Mode I: save events with pre- and post-event time
            # with data from all triggered stations
            
            # suffers from same issue as mode 2 ! seemingly random
            # stations with very poor time cropping!!
                
            stream = obspy.Stream()
                
#            for stream_file in stream_files:
#                    
#                stream += obspy.read(stream_file)
#        
#            for coincidence_trigger in coincidence_triggers:
#               
#                event_stream = obspy.Stream()
#                
#                print('Saving ' + str(coincidence_trigger))
#                
#                for trace in stream:
#                        
#                    if (trace.stats.station in coincidence_trigger[2]) and (trace.stats.channel[-1] == stream_component):
#
#                        event_stream += trace.copy().trim(starttime = coincidence_trigger[0] - pre_event_time,
#                                           endtime = coincidence_trigger[1] + post_event_time)
#
#                event_stream.write(event_output_directory + str(coincidence_trigger[0]) + '.MSEED',
#                             format = 'MSEED', cencoding = 'STEIM2')

#            # Mode II: save events with each stream trimmed to its trigger duration

            for event_triggers in all_event_triggers:
                
                print('Saving ' + str(coincidence_triggers[all_event_triggers.index(event_triggers)]))
                
                event_stream = obspy.Stream()
                
                for event_trigger in event_triggers:
                    
                    event_stream += obspy.read(glob.glob(stream_root_directory + 'Y' + str(year) + '/R' + str(doy) + '.01/' + event_trigger[2] + '*' + stream_component + '*')[0],
                                               starttime = event_trigger[0] - pre_post_time * FFT_window_len, endtime = event_trigger[1] + pre_post_time * FFT_window_len)
    
                event_stream.write(event_output_directory + str(coincidence_triggers[all_event_triggers.index(event_triggers)][0]) + \
                                   '.MSEED', format = 'MSEED', cencoding = 'STEIM2')