# -*- coding: utf-8 -*-
"""
Qualatative seismic event inspection code. This code uses continuous 3-component
seismic streams and seismic event files derived from an STA:LTA detection algorithm
to create combined component displacement timeseries, spectrograms, and particle
motion plots with sonified combined displacement timeseries overlaid as audio.

Dependencies: non-default python modules: matplotlib, numpy, obspy, scipy, and
ffmpeg installed on a linux operating system.
"""

import datetime
import math
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
import obspy
from obspy.signal.polarization import polarization_analysis
from obspy.imaging.spectrogram import spectrogram
import os
import scipy
import scipy.io.wavfile as wav
import subprocess

# Set parameters

## Catalogue root directory contains folders each containing event files

catalogue_directory = '/media/sam/61D05F6577F6DB39/SCIENCE/paper_seismic_catalogue/'
catalogues = os.listdir(catalogue_directory)

event_file_suffix = '.MSEED' # set suffix for event files (allows directory to contain non-event files)

### Catalogues and stations to ignore

ignore_catalogues = []
ignore_stations = ['TSNM1', 'TSNM2', 'TSNM3']

## Set only some files to process (catalogue independent)

filelist = []

## Set a limit to the number of events to process per catalogue

process_limit = 3

## Stream root directory contains individual day-long streams of each station's components
## within julian day directories

stream_directory = '/media/sam/61D05F6577F6DB39/SCIENCE/day_volumes_S/'
stream_years = os.walk(stream_directory).next()[1]
year_doys = [[] for year in stream_years]
for year in stream_years:   year_doys[stream_years.index(year)] = os.walk(stream_directory+year+'/').next()[1]
    
earliest_time = '2016-04-21T00:00:00.000Z' # earliest time in streams
    
## The response function file in RESP format. Ensure this is generalised to the 'temp' sensor metadata used later in this code.
    
respf = '/home/sam/Downloads/RESP.XX.NS361..SHZ.L28LB.395.2490.4_5.27'
pre_filt = (1.0, 10.0, 120.0, 125.0) # set filter to apply before response removal to reduce noise introduction
    
## Parameters for event output    
    
pre_event_time = 2 # desired pre-event time
post_event_time = 4 # desired post-event time

old_poet = 5 # event file post-event time from STA:LTA detection algorithm
    
#! CODE BEGINS !#

# Run through each catalogue
    
for catalogue in catalogues:
    
    if catalogue in ignore_catalogues: continue # ignore certain catalogues
    
    events = os.listdir(catalogue_directory + catalogue) # get all event files in catalogue directory
   
    p = -1   
   
    for event in events: # process all files in catalogue directory
    
            if filelist: # only engage when a filelist is given
                
                if event not in filelist: continue # only process some events

            # Check file is an event file

            if event[-(len(event_file_suffix)):] != event_file_suffix: continue # only process event files
            
            p += 1
            
            if p == process_limit: break # stop processing catalogue events once the event limit is reached
            
            print '\n Processing event ' + event + ' in catalogue ' + catalogue + '\n'

            # Get steam for day of event
            
            year = event[:4]
            time = datetime.datetime.strptime(event[:-7], '%Y-%m-%dT%H:%M:%S.%f')
            # NOTE: assumes an ISO8601 format to event onset times contained in the event file name
            doy = str(time.timetuple().tm_yday)
             
            stations = [] 
            raw_stream = obspy.Stream()
            raw_stream = raw_stream.clear()
            for component_stream in os.listdir(stream_directory + 'Y' + year + '/' + 'R' + doy + '.01/'):
                station = component_stream.split('.')[0]
                
                ## Ignore stations with certain codes                
                
                if station in ignore_stations: continue
                
                stations.append(station)
                raw_stream += obspy.read(stream_directory + 'Y' + year + '/' + 'R' + doy + '.01/' + component_stream)
                
            stations = list(set(stations)) # remove repetition of stations
                
            # Remake event from stream with desired pre- and post-event times
                
            raw_stream.trim(starttime = (obspy.UTCDateTime(time) - pre_event_time) - 3, endtime = (obspy.UTCDateTime(time) + post_event_time) + 1) 
            # NOTE: -3 and +1 are hardcoded for particle motion calculation trimming the stream

            # Remove instrument response

            ## Prior to instrument response removal store all trace metadata
            ## for reapplication after using a general RESP file
            
            tr_stat = []
            tr_net = []
            tr_loc = []
            tr_chan = []
            tr_star = []
            tr_end = []                
            
            for trace in raw_stream:
                
                tr_stat.append(trace.stats.station)
                tr_net.append(trace.stats.network)
                tr_loc.append(trace.stats.location)
                tr_chan.append(trace.stats.channel)
                tr_star.append(trace.stats.starttime)
                tr_end.append(trace.stats.endtime)
                
                trace.stats.station = 'temp'
                trace.stats.network = 'temp'
                trace.stats.location = 'temp'
                trace.stats.channel = 'temp'
                trace.stats.starttime = earliest_time  
            
            # the actual response removal

            seedresp = {'filename': respf, 'units': 'DIS'}                
            stream = raw_stream.simulate(paz_remove=None, pre_filt=pre_filt, seedresp=seedresp)
            
            ## Restore trace metadata
            
            t = 0
            for trace in stream:
                
                trace.stats.station = tr_stat[t]
                trace.stats.network = tr_net[t]
                trace.stats.location = tr_loc[t]
                trace.stats.channel = tr_chan[t]
                trace.stats.starttime = tr_star[t]
                
                t += 1                  
            
            stream.detrend(type = 'simple')
                            
            # Combine components into one stream
            
            global sampling_rate
            combined_streams = [[] for i in range(len(stations))]
            for station in stations:
                combined_stream = [0] * len(stream[0])
                for trace in stream:
                    if trace.stats.station == station:
                        for i in range(len(trace.data)):
                            sampling_rate = float(trace.stats.sampling_rate)
                            combined_stream[i] += trace.data[i]
                combined_streams[stations.index(station)] = combined_stream        
            
            # Calculate stream SNR from event file and continue processing
            # with combined stream for the station with highest event file SNR
            
            eventstream = obspy.read(catalogue_directory + catalogue + '/' + event)
            endtime = eventstream[0].stats.endtime
            
            pre_event_samples=int(pre_event_time * sampling_rate)
            event_samples=int(((endtime - old_poet) - obspy.UTCDateTime(time)) * sampling_rate)
            
            SNRs = []
            for combined_stream in combined_streams:
            
                noise=np.abs(combined_stream[:pre_event_samples - int(0.5 * sampling_rate)]) # take noise as all data up until 0.5 second before the event
                signal=np.abs(combined_stream[pre_event_samples : pre_event_samples + event_samples]) # take signal as all data during the event
                noise_mean=sum(noise) / float(len(noise))
                signal_mean=sum(signal) / float(len(signal))       
                SNRs.append(signal_mean / noise_mean)
                
            seismogram = combined_streams[SNRs.index(max(SNRs))]  / max([abs(min(combined_streams[SNRs.index(max(SNRs))])), abs(max(combined_streams[SNRs.index(max(SNRs))]))]) # normalization applied
            station = stations[SNRs.index(max(SNRs))]
            
            # Keep 3 component data separately for use in particle motion calculation
            
            seismogram_components = obspy.Stream()
            seismogram_components = seismogram_components.clear()
            seismogram_components += stream[tr_stat.index(station) + 2]
            seismogram_components += stream[tr_stat.index(station) + 1]
            seismogram_components += stream[tr_stat.index(station) + 0]
            
            # Generate figure for output
            
            print 'Plotting to figure - if necessary please change plot axes labels and ticks for your data in the code body.'            
            
            ## Add combined seismogram to figure
            
            rel_time = np.linspace(0, obspy.UTCDateTime(tr_end[tr_stat.index(station)]) - obspy.UTCDateTime(tr_star[tr_stat.index(station)]), len(seismogram)).tolist()
            rel_ticks = range(0, int(math.floor((obspy.UTCDateTime(tr_end[tr_stat.index(station)]) - obspy.UTCDateTime(tr_star[tr_stat.index(station)])))))
        
            plt.figure(figsize = [7.55906, 9.448825])
            # NOTE: figure size is hardcoded to be a nice size for desktop computer screens
            ax1 = plt.subplot(411)
            ax1.plot(rel_time, 0.99 * seismogram, color = 'k')
            ax1.set_xlim(3, rel_time[-1] - 1.5)                 
            ax1.set_yticks([-1, 0, 1])
            ax1.set_ylim([-1.1, 1.1])
            ax1.spines['top'].set_visible(True)
            ax1.spines['right'].set_visible(True)
            ax1.yaxis.set_ticks_position('both')
            ax1.xaxis.set_ticks_position('both')
            ax1.tick_params(which = 'both', direction = 'out')
            ax1.tick_params(which = 'both', length = 5)
            ax1.tick_params(which = 'both', width = 1.1)
            ax1.tick_params(axis='x',which='both',labelsize=0)
            ax1.set_ylabel('normalised\namplitude\n', fontsize = 12)
            ax1.get_yaxis().set_label_coords(-0.08,0.5)        
            
            plt.title(str.upper(catalogue) + ' : ' + event[:-6], fontsize = 12, y = 1.1)
            
            ## Generate spectrogram for combined seismogram and add to figure         
            
            print 'Generating spectrogram plots for station ' + station

            ax2 = plt.subplot(412)          
            spectrogram(np.array(seismogram), samp_rate = sampling_rate, per_lap = 0.99, wlen = 1.0, cmap='jet', axes = ax2)
            ax2.set_xlim(3, rel_time[-1] - 1.5)
            ax2.set_yticks([0, 40, 80, 120])
# HARDCODE !!!
            ax2.spines['top'].set_visible(True)
            ax2.spines['right'].set_visible(True)
            ax2.yaxis.set_ticks_position('both')
            ax2.xaxis.set_ticks_position('both')
            ax2.tick_params(which = 'both', direction = 'out')
            ax2.tick_params(which = 'both', length = 5)
            ax2.tick_params(which = 'both', width = 1.1)
            ax2.tick_params(axis='x',which='both',labelsize=0)
            ax2.set_ylabel('frequency\n(Hz)\n', fontsize = 12)
            ax2.get_yaxis().set_label_coords(-0.08,0.5)

            ## Generate polarization plots and add to figure

            print 'Generating polarization plots for station ' + station
            
            particle_motion = polarization_analysis(seismogram_components,\
            win_len = 1.0, win_frac = 0.01, frqlow = 1.0, frqhigh = 120.0, method = 'vidale',\
# HARDCODE !!!
            stime = seismogram_components[0].stats.starttime, etime = seismogram_components[0].stats.endtime)
            
            pm_rel_time = np.linspace((obspy.UTCDateTime(particle_motion['timestamp'][0]) - obspy.UTCDateTime(tr_star[tr_stat.index(station)])), \
            ((obspy.UTCDateTime(particle_motion['timestamp'][0]) - obspy.UTCDateTime(tr_star[tr_stat.index(station)])) +\
            (obspy.UTCDateTime(particle_motion['timestamp'][-1]) - obspy.UTCDateTime(particle_motion['timestamp'][0]))), len(particle_motion['timestamp']))          
            
            ax3 = plt.subplot(413)
            ax3.plot(pm_rel_time, particle_motion['azimuth'], color = 'k')
            ax3.set_xlim(3, rel_time[-1] - 1.5)  
            ax3.set_yticks([0, 60, 120, 180])
# HARDCODE !!!
            ax3.set_ylim([-10, 190])
            ax3.set_xlim(3, rel_time[-1] - 1.5)        
            ax3.spines['top'].set_visible(True)
            ax3.spines['right'].set_visible(True)
            ax3.yaxis.set_ticks_position('both')
            ax3.xaxis.set_ticks_position('both')
            ax3.tick_params(which = 'both', direction = 'out')
            ax3.tick_params(which = 'both', length = 5)
            ax3.tick_params(which = 'both', width = 1.1)
            ax3.tick_params(axis='x',which='both',labelsize=0)
            ax3.set_ylabel('polarization\nazimuth\n(deg)', fontsize = 12)
            ax3.get_yaxis().set_label_coords(-0.08,0.5)
            
            ax4 = plt.subplot(414)
            ax4.plot(pm_rel_time, particle_motion['incidence'], color = 'k')
            ax4.set_xlim(3, rel_time[-1] - 1.5)  
            ax4.set_xticklabels(range(len(rel_ticks[rel_ticks.index(3) : rel_ticks.index(int(math.ceil(rel_time[-1] - 1.5)))])))
            ax4.set_yticks([0, 30, 60, 90])
# HARDCODE !!!
            ax4.set_ylim([-5, 95])
            ax4.set_xlim(3, rel_time[-1] - 1.5)
            ax4.spines['top'].set_visible(True)
            ax4.spines['right'].set_visible(True)
            ax4.yaxis.set_ticks_position('both')
            ax4.xaxis.set_ticks_position('both')
            ax4.tick_params(which = 'both', direction = 'out')
            ax4.tick_params(which = 'both', length = 5)
            ax4.tick_params(which = 'both', width = 1.1)
            ax4.set_xlabel('time (s)', fontsize = 12)
            ax4.get_xaxis().set_label_coords(0.5,-0.2)
            ax4.set_ylabel('polarization\ninclination\n(deg)', fontsize = 12)
            ax4.get_yaxis().set_label_coords(-0.08,0.5)
            
            plt.tight_layout()

            print 'Saving plots as .mp4 video'

            # Reduce combined seismogram data to the length shown in plot
            
            seismogram = seismogram[int(sampling_rate * 3) : -int(sampling_rate * 1.5)]

            # animate plot and save as .mp4
            
            fig = plt.gcf()            
            
            axlines = []
            for t in range(len(seismogram)):
                t /= sampling_rate
                t += 3
                ax1lines, = ax1.plot([t] * 4000, range(-2000, 2000), color = 'red', linewidth = 1)
                ax2lines, = ax2.plot([t] * 4000, range(-2000, 2000), color = 'red', linewidth = 1)
                ax3lines, = ax3.plot([t] * 4000, range(-2000, 2000), color = 'red', linewidth = 1)
                ax4lines, = ax4.plot([t] * 4000, range(-2000, 2000), color = 'red', linewidth = 1)
                axlines.append([ax1lines, ax2lines, ax3lines, ax4lines])
            myanimation = animation.ArtistAnimation(fig, axlines, interval = 1000 * (1 / sampling_rate), repeat = False)
            myanimation.save(catalogue_directory + catalogue + '/' + station + '_' + event[:event.index('.MSEED')] + '.mp4', writer = 'ffmpeg', codec = 'libx264', bitrate = -1)
            
            print 'Audioifying all-component seismogram for station ' + station
            
            # Prepare to shift seismogram into widest audible range

            stretch_factor = 20000 / (sampling_rate)                
            grain_duration = 1 / sampling_rate 
                
            # Calculate grain width and number of grains in seismogram
                
            grain_width = int(math.floor(grain_duration * sampling_rate))
            num_grains = int(math.floor(len(seismogram) / float(grain_width)))
            # NOTE: data loss will occur if the above fraction does not evaluate to an integer
            
            # Generate grain spectra and store it for interpolation        
        
            all_grain_spectras = [[] for i in range(num_grains)]
            for i in range(1, num_grains):
                
                grain = seismogram[(i - 1) * grain_width : i * grain_width]
                grain_spectra = np.fft.fft(grain)            
                all_grain_spectras[i-1] = grain_spectra
                
            # Generate interpolated grain spectra
    
            interpolated_grain_spectras = [[] for i in range(len(all_grain_spectras[0]))]
            for j in range(len(all_grain_spectras[0])):
                
                # Create spectra timeseries across all grains            
                
                spectral_timeseries = []            
                for i in range(1, num_grains-1):    
                    
                    spectral_timeseries.append(all_grain_spectras[ i - 1 ][j])
                
                spectra_function = scipy.interpolate.interp1d(range(len(spectral_timeseries)), spectral_timeseries)
                interpolated_spectra_timeseries = spectra_function(np.linspace(0, len(spectral_timeseries)-1, (stretch_factor * len(spectral_timeseries))))             
                interpolated_grain_spectras[j] = interpolated_spectra_timeseries
            
            # Regenerate seismogram from grain spectra
            
            interpolated_seismogram = []
            reduced_IS = []
            for i in range(len(interpolated_grain_spectras[0])):
                
                # Gather frequency data for each grain
                
                interpolated_grain_spectra = []
                for j in range(len(interpolated_grain_spectras)):
                    
                    interpolated_grain_spectra.append(interpolated_grain_spectras[j][i])
                    
                interpolated_seismogram.extend(np.real(np.fft.ifft(interpolated_grain_spectra)))
                
            # Convert seismogram to 16bit integer format so it plays correctly
            # with .wav file handling audio software, then save as .wav file
            # with speed-up to audiofy the seismogram.            
            
            interpolated_seismogram = 0.99 * ( interpolated_seismogram / max(abs(min(interpolated_seismogram)), max(interpolated_seismogram)) )
            interpolated_seismogram_16bit = np.int16(np.array(interpolated_seismogram) * float(2**15))
            wav.write(catalogue_directory + catalogue + '/' + station + '_' + event[:event.index('.MSEED')] + '.wav', (stretch_factor * sampling_rate), interpolated_seismogram_16bit)

            # Combine audio and visual then remove intermediate files

            subprocess.call('ffmpeg -i '+ catalogue_directory + catalogue + '/' + station + '_' + event[:event.index('.MSEED')] + '.mp4 -i ' + \
            catalogue_directory + catalogue + '/' + station + '_' + event[:event.index('.MSEED')] + '.wav -strict experimental ' + \
            catalogue_directory + catalogue + '/temp_' + station + '_' + event[:event.index('.MSEED')] + '.mp4', shell = True)
            
            subprocess.call('mv ' + catalogue_directory + catalogue + '/temp_' + station + '_' + event[:event.index('.MSEED')] + '.mp4 ' + \
            catalogue_directory + catalogue + '/' + station + '_' + event[:event.index('.MSEED')] + '.mp4', shell = True)

            subprocess.call('rm -rf ' + catalogue_directory + catalogue + '/' + station + '_' + event[:event.index('.MSEED')] + '.wav', shell = True)