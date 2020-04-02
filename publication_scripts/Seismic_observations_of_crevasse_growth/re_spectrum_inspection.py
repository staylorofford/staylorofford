#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Frequency domain seismic detection algorithm,
uses spectrum files saved from spectrum_generation.
"""

# Import packages

import datetime
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.animation as animation

# Set parameters

## Directory to load spectrum files from

#spectrum_directory = '/home/sam/spectrums/'
#spectrum_directory = '/media/sam/Seagate Backup Plus Drive/SCIENCE/Spectrums/'
spectrum_directory = '/home/sam/EVENTS_IT3/TYPE_A/5/'

## Stations to process

stream_stations = ['TSNC1', 'TSNC3', 'TSNL2', 'TSNL3', 'TSNR2', 'TSNR3']#, 'TSNM1', 'TSNM3']

## Start and end dates for processing of day-long seismic streams

start_year = '2016'
start_month = '04'
start_day = '01'

end_year = '2017'
end_month = '08'
end_day = '01'
            
## window length to plot in FFT windows

window_length = 15

## window overlap as a fractional percentage

window_overlap = 0.99999

## sampling details (samples in samples, sampling rate in Hz)

samples = 50
sampling_rate = 250

## Index focus (removes some frequencies from plotting)

indices = [0, 100] 

files = os.listdir(spectrum_directory)
events = []
for afile in files:
    if afile[-6:] == '.MSEED':
        events.append(afile)
        
events.sort()
                
# Generate trigger times

for event in events:
            
            spectrum_files = glob.glob(spectrum_directory + event[:-6] + '*spectrums.npy')
            
            all_spectrums = []
            all_times = []
            plot_stations = []
            
            for spectrum_file in spectrum_files:
                
                stream_file_metadata = spectrum_file.split('/')[-1].split('_')
                station = stream_file_metadata[-2]
                
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
                            
                all_times.append(times)
                all_spectrums.append(np.array(spectrums))
                plot_stations.append(station)
               
            try:
                
                times = np.array(times)
                spectrums = np.array(spectrums)
            
            except:
                
                continue
            
            # Plot the spectrum as an "animation" using refreshes of an imshow plot
            
            num_windows = int(float(1 / (1 - window_overlap)) * \
              (len(spectrums)) / float(window_length))
            
            frequencies = np.fft.fftfreq(samples, 1 / sampling_rate)
            
            fig = plt.figure(figsize = (8,8))
            
            im_list = [[] for j in range(len(all_spectrums))]
            
            numplots = len(all_spectrums)
            
            if numplots / 6 > 1:
                
                rows = 3
                cols = 3
                
            elif numplots / 3 > 1:
                
                rows = 2
                cols = 3
                
            else:
                
                rows = 1
                cols = 3
            
            j = -1
            for spectrum in all_spectrums:
                j += 1
                
                plt.subplot(str(rows) + str(cols) + str(j))
                # Uncomment to alter data, so each window scales relative to
                # its highest data value.
                maxval = max(spectrum[: window_length].flatten())
                im_list[j] = plt.imshow(spectrum[: window_length].transpose() / maxval, extent = [0, window_length, \
                                0, len(frequencies) - 1], origin = 'lower', cmap=cm.jet, aspect = 0.5, \
                                vmin = 0, vmax = 1)     
                # Uncomment to leave data alone: plotting colorbar will be relative to the highest value in all windows
#                im_list[j] = plt.imshow(spectrum[: window_length].transpose(), extent = [0, window_length, \
#                                frequencies[0], frequencies[-1]], origin = 'lower', cmap=cm.jet, aspect = 0.5)
                
                plt.title(plot_stations[j])
                plt.xticks([0, window_length])
                plt.yticks(range(len(frequencies)), frequencies)
#                plt.xticks([0, window_length], [all_times[j][0], all_times[j][-1]], rotation = 90)  
                plt.xticks(range(window_length), all_times[j], rotation = 90)
                
            def run_animation():
                
                '''
                Function from woodenflute on stackoverflow.com
                to allow for animation pausing.
                '''
                
                anim_running = True
            
                def onClick(event):
                    nonlocal anim_running
                    if anim_running:
                        spectrogram_animation.event_source.stop()
                        anim_running = False
                    else:
                        spectrogram_animation.event_source.start()
                        anim_running = True
                                
                def update_spectrogram(w):
        
                    '''
                    Function to update the current figure with
                    each successive spectrogram.
                    '''
                    
                    i = int(w * (1 - window_overlap) * window_length)
                    for j in range(len(all_spectrums)):
                        maxval = max(all_spectrums[j][i : i + window_length].flatten())
                        im_list[j].set_array(all_spectrums[j][i : i + window_length].transpose() / maxval)
                    fig.suptitle(times[i + int(window_length / 2)])
                        
                    return im_list
            
                fig.canvas.mpl_connect('button_press_event', onClick)
            
                spectrogram_animation = animation.FuncAnimation(fig, update_spectrogram, frames = range(1, num_windows), 
                                        interval=200)
                
                plt.show()
            
            run_animation()