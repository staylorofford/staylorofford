# -*- coding: utf-8 -*-
"""
Take a single waveform and decompose it,
overlapping the composite frequencies with
each composite waveform having its amplitude
from a FFT of the input waveform.

If there is more than one waveform in the
input file, up to 9 waveforms will be decomposed.

If there is more than one component in the 
input file, up to 3 components will be decomposed.
"""

# Import packages

import argparse
import obspy
import numpy
import matplotlib.pyplot as plt

# Get arguments

parser = argparse.ArgumentParser(description = \
        'Decompose a waveform into a superposition of its "composite waveforms" ' + \
        'generated from the output of a Fourier transform.')
parser.add_argument('file_path', nargs = 1, help = 'File path of waveform', type = str)
parser.add_argument('pre_event_time', nargs = 1, help = 'Surplus data before event waveform', type = float)
parser.add_argument('post_event_time', nargs = 1, help = 'Surplus data after event waveform', type = float)

try:
    args = parser.parse_args()
except:
    parser.print_help()
    
file_path = args.file_path[0]
pret = args.pre_event_time[0]
poet = args.post_event_time[0]

# Load waveform and prepare it for FFT application

waveforms = obspy.read(file_path)
sampling_rate = waveforms[0].stats.sampling_rate
waveform_duration = len(waveforms[0].data)
waveforms.filter('highpass', freq = 1 / float(waveform_duration))
waveforms.detrend(type = 'demean')
waveforms.detrend(type = 'simple')

# Calculate # waveforms to process

for component in ['Z', ['N', '1'], ['E', '2']]:

    wmax = -1
    for waveform in waveforms:
    
        if waveform.stats.station[3] == 'M': continue #ignore moraine sites
# HARDCODE !!!
        if waveform.stats.channel[-1] not in component: continue
        wmax +=1
        
for component in ['Z', ['N', '1'], ['E', '2']]:

    plt.figure()

    w = -1
    for waveform in waveforms:
    
        if waveform.stats.station[3] == 'M': continue #ignore moraine sites
# HARDCODE !!!
        if waveform.stats.channel[-1] not in component: continue
        w +=1
        
        # Cut waveform to event
        
        try:
            waveform.data = waveform.data[int(pret * sampling_rate) : waveform_duration - int(poet * sampling_rate)]
        except:
            pass
        
        # Optional: show waveform
        
        #plt.plot(waveform[0])
        #plt.show()
        
        # Calculate FFT of modified waveform
        
        spectrum = numpy.fft.fft(waveform.data).real.tolist()
        frequencies = numpy.fft.fftfreq(len(waveform.data), 1 / float(sampling_rate)).tolist()
        
        # Optional: trim to real component and plot
        
        nyquist_index = frequencies.index(-sampling_rate / 2)
        spectrum = spectrum[1 : nyquist_index]
        frequencies = frequencies[1 : nyquist_index]
        #plt.plot(frequencies, spectrum)
        #plt.show()
        
        # Normalise spectrum values
        
        maxval = max(spectrum)
        for i in range(len(spectrum)):
            
            spectrum[i] = spectrum[i] / maxval
        
        # Decompose waveform
        
        composite_waveforms = [[] for i in range(len(frequencies))]
        t = numpy.linspace(0, len(waveform.data) / sampling_rate, len(waveform.data))
        
        for i in range(int(sampling_rate / 2)):
            
            sinusoid = spectrum[i] * numpy.sin(2 * numpy.pi * frequencies[i] * t)
            composite_waveforms[i] = sinusoid
        
        # Plot all frequencies, showing their relative amplitude differences in
        # the waveforms.
        
        ## Determine subplot layout

        if wmax / 3 > 1:
            rows = 3
            cols = 2
        else:
            rows = wmax
            cols = 1
            
        if wmax / 6 > 1:
            rows = 3
            cols = 3
            
        c = -1
        plt.subplot(str(rows) + str(cols) + str(w))
        for i in range(int(sampling_rate / 2)):
            c += 1    
            
            plt.plot(t[:int(sampling_rate)], composite_waveforms[i][:int(sampling_rate)] + c,
                    color = 'k', alpha = 0.5)
    
#        plt.xlabel('Time (s)', labelpad = 15, fontsize = 10)
        plt.xticks([])
        plt.ylabel('Frequency (Hz)', labelpad = 10, fontsize = 10)
        plt.title(waveform.stats.station + ' - ' + waveform.stats.channel + \
        ' - ' + str(waveform.stats.starttime + pret), fontsize = 10)
        
plt.show()