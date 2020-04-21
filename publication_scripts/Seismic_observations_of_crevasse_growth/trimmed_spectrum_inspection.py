#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Frequency domain seismic detection algorithm,
uses spectrum files saved from spectrum_generation.
"""

# Import packages

import datetime
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.animation as animation

# Set parameters

## Directory to save event files to

#event_output_directory = '/media/sam/61D05F6577F6DB39/SCIENCE/FDD/'

## Directory to load spectrum files from

spectrum_directory = '/home/samto/Spectrums_it4/'

## Seismic component to use in processing

stream_component = 'Z'

## Stations to process

stream_stations = ['TSNM1', 'TSNM2', 'TSNM3']

## Start and end dates for processing of day-long seismic streams

start_year = '2016'
start_month = '05'
start_day = '01'

end_year = '2017'
end_month = '05'
end_day = '19'
            
## window length to plot in FFT windows

window_length = 300

## window overlap as a fractional percentage

window_overlap = 0.99999

## Index focus (removes some frequencies from plotting)

indices = [0, -1]

## FFT window length

sampling_rate = 250
FFT_win_len = 250000

FFT_freqs = np.fft.fftfreq(FFT_win_len, 1/sampling_rate)[:int(sampling_rate / 2) + 1][1:][indices[0]:indices[1]]
print(FFT_freqs)
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
        
#         The logic around more than one day of animation hasn't been done yet
#         so use this crappy conditional to choose your day!
        
#        if doy not in [125, 132, 133, 134, 135, 136, 137]: continue
    
        if (year == int(start_year)) and (doy < int(start_date_doy)): continue                        
        elif (year == int(end_year)) and (doy > int(end_date_doy)): continue            
        else:
            
            spectrum_files = glob.glob(spectrum_directory + '*' + str(year) + '*' + str(doy) + '*spectrums.npy')

            all_spectrums = []
            plot_stations = []
            
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
                
                # trim spectrums to set indices
                
                trimmed_spectrums = []
                s = -1
                for spectrum in spectrums:
                    s += 1

                    trimmed_spectrums.append(spectrum[indices[0] : indices[1]])
                    
                all_spectrums.append(np.array(trimmed_spectrums))
                plot_stations.append(station)
               
            try:
                
                times = np.array(times)
                spectrums = np.array(spectrums)
            
            except:
                
                continue
            
            # Plot the spectrum as an "animation" using refreshes of an imshow plot
            
            num_windows = int(float(1 / (1 - window_overlap)) * \
              (len(spectrums)) / float(window_length))
            
#            frequencies = range(len(spectrums[0]))
## WARNING - frequencies can only be defined like this for 1 second FFT window lengths
            
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
                
                rows = 3
                cols = 1    
            
            j = -1
            for spectrum in all_spectrums:
                j += 1
                
                plt.subplot(str(rows) + str(cols) + str(j))
                # Uncomment to alter data, so each window scales relative to
                # its highest data value.
                maxval = max(spectrum[: window_length].flatten())
                im_list[j] = plt.imshow(spectrum[: window_length].transpose() / maxval, extent = [0, window_length, \
                                0, len(FFT_freqs)], origin = 'lower', cmap=cm.jet, aspect = 0.5, \
                                vmin = 0, vmax = 1)     
                # Uncomment to leave data alone: plotting colorbar will be relative to the highest value in all windows
#                im_list[j] = plt.imshow(spectrum[: window_length].transpose(), extent = [0, window_length, \
#                                frequencies[0], frequencies[-1]], origin = 'lower', cmap=cm.jet, aspect = 0.5)
                
                plt.title(plot_stations[j])
                plt.xticks([0, window_length])
                plt.yticks(range(0, len(FFT_freqs)))
                
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
                        try:
                            maxval = max(all_spectrums[j][i : i + window_length].flatten())
                            im_list[j].set_array(all_spectrums[j][i : i + window_length].transpose() / maxval)
                        except:
                            pass
                    try:
                        fig.suptitle(times[i + int(window_length / 2)])
                    except:
                        pass
                        
                    return im_list
            
                fig.canvas.mpl_connect('button_press_event', onClick)
            
                spectrogram_animation = animation.FuncAnimation(fig, update_spectrogram, frames = range(1, num_windows), 
                                        interval=200)
                
                plt.show()
            
            run_animation()