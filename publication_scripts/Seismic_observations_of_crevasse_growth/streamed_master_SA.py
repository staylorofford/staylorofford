# -*- coding: utf-8 -*-
"""
Master script to use audification of seismograms code library.
"""

##### This code was written in a miniconda virtual enivornment
##### with obspy, numpy, portaudio, pyaudio, scipy, and matplotlib installed

import glob
import math
import matplotlib.pyplot as plt
import numpy as np
import obspy
import os
import scipy
import scipy.io.wavfile as wav

# Set parameters

buffer_length = 10 # number of seconds of seismogram to process in buffer

#interest_frequency_range = [5, 20]
#stretch_factor = 20000 / interest_frequency_range[1]

# Parse in data structure: root directory, day of year (doy) folders,
# and seismic station names.

seismic_data_root = "/media/sam/61D05F6577F6DB39/SCIENCE/day_volumes_S/Y2016/"
seismic_data_doys = os.listdir(seismic_data_root)
streams0 = os.listdir(seismic_data_root+seismic_data_doys[0])
stations = []
for stream in streams0:
    stations.append(stream[ : stream.index(".")]) # station name ends at first .
stations = list(set(stations)) # keep only unique names
del streams0

# Tier 1 loop to run through all doy folders

for doy in seismic_data_doys:
    
    print 'Loading data from day '+doy+'\n'
    
    # Tier 2 loop to load in all Z component seismograms in each doy folder
    
    for station in stations:
        
        print 'Processing seismogram from '+station+'\n'
#
#         Load and normalize seismogram with 'station*Z*' filename format

        fullstream = obspy.read(glob.glob(seismic_data_root+doy+'/'+station+"*Z*")[0])
        fullstream.detrend(type = 'simple')
        
        for trace in fullstream: 
            start_time = trace.stats.starttime
            sampling_rate = float(trace.stats.sampling_rate)

            stretch_factor = 20000 / (sampling_rate)                
            grain_duration = 1 / sampling_rate 
        
        # Test loop for >60 s sound generation        
        
        for t in range(int(math.floor(86400.0 / buffer_length))):
            
            print 'Generating the audio file for the '+str(buffer_length)+' second-long window with the '+str(t * buffer_length)+'th starting second'
            
            stream = fullstream.slice(starttime = (start_time + t * buffer_length), endtime = (start_time + (t + 1) * buffer_length))
            for trace in stream: seismogram = trace.data
        
    #        stream = obspy.read('/home/sam/Scripts/2016-04-26T14_36_34.900000Z.MSEED')
                
            # Calculate grain width and number of grains in seismogram
                
            grain_width = int(math.floor(grain_duration * sampling_rate))
            num_grains = int(math.floor(len(seismogram) / float(grain_width)))
            # NOTE: data loss will occur if the above fraction does not evaluate to an integer
            
            # Tier 3 loop to generate grain spectra and store it for interpolation        
        
            all_grain_spectras = [[] for i in range(num_grains)]
            for i in range(1, num_grains):
                
                grain = seismogram[(i - 1) * grain_width : i * grain_width]
                grain_spectra = np.fft.fft(grain)            
                all_grain_spectras[i-1] = grain_spectra
                
            # Tier 3 loop to generate interpolated grain spectra
    
            interpolated_grain_spectras = [[] for i in range(len(all_grain_spectras[0]))]
            for j in range(len(all_grain_spectras[0])):
                
                # Tier 4 loop to create spectra timeseries across all grains            
                
                spectral_timeseries = []            
                for i in range(1, num_grains-1):    
                    
                    spectral_timeseries.append(all_grain_spectras[ i - 1 ][j])
                
                spectra_function = scipy.interpolate.interp1d(range(len(spectral_timeseries)), spectral_timeseries)
                interpolated_spectra_timeseries = spectra_function(np.linspace(0, len(spectral_timeseries)-1, (stretch_factor * len(spectral_timeseries))))             
                interpolated_grain_spectras[j] = interpolated_spectra_timeseries
            
            # Tier 3 loop to regenerate seismogram from grain spectra
            
            interpolated_seismogram = []
            reduced_IS = []
            for i in range(len(interpolated_grain_spectras[0])):
                
                # Tier 4 loop to gather frequency data for each grain
                
                interpolated_grain_spectra = []
                for j in range(len(interpolated_grain_spectras)):
                    
                    interpolated_grain_spectra.append(interpolated_grain_spectras[j][i])
                    
                interpolated_seismogram.extend(np.real(np.fft.ifft(interpolated_grain_spectra)))
    
                # Reduced interpolated seismogram to parallel sampling frequency of original seismogram
#                if (i % stretch_factor) == 0:
#                    reduced_IS.extend(np.real(np.fft.ifft(interpolated_grain_spectra)))

            # Seismogram interpolation QC plots
                
    #        plt.plot(interpolated_seismogram)
    #        plt.figure()
    #        plt.plot(reduced_IS)
    #        plt.figure()
    #        plt.plot(seismogram)
    #        plt.show()
                
            # Convert seismogram to 16bit integer format so it plays correctly
            # with .wav file handling audio software, then save as .wav file
            # with speed-up to audiofy the seismogram.            
            
            interpolated_seismogram = 0.99 * ( interpolated_seismogram / max(abs(min(interpolated_seismogram)), max(interpolated_seismogram)) )
            interpolated_seismogram_16bit = np.int16(np.array(interpolated_seismogram) * float(2**15))
            wav.write(station+'_'+doy+'_'+str(t)+'.wav', (stretch_factor * sampling_rate), interpolated_seismogram_16bit)

        break
    
    break