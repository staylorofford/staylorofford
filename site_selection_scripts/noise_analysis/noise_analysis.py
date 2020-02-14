# -*- coding: utf-8 -*-
"""
Combine local data from test sites [to do: with data from a reference site (queried via FDSN)]
and perform noise analysis and plotting on it.
"""

import os
import obspy
from obspy.core.stream import Stream
from obspy.imaging.cm import pqlx
from obspy.io.xseed import Parser
from obspy.signal import PPSD

# Define root data directory

day_file_directory_root = '/home/samto/PROCESSING/SS_geartest/'

# Define sensor + datalogger metadata (dataless SEED) directory

metadata_directory = '/home/samto/git/staylorofford/site_selection_scripts/metadata/'

# Parse all miniseed files under the root data directory into an obspy Stream object
sites = []
all_streams = Stream()
for root, dirs, files in os.walk(day_file_directory_root):
    for file in files:
        if file[-3:] == '.ms':
            sites.append(file.split('.')[1])
            all_streams += obspy.read(root + '/' + file)
sites = list(set(sites))  # Define site list
all_streams.merge(method=1)  # Merge streams to one per channel

# Prepare list of streams so stream object can be separated by site
streams = []
for n in range(len(sites)):
    streams.append(Stream())

# Split stream object into site streams
for stream in all_streams:
    site = stream.stats.station
    streams[sites.index(site)] += stream

# Build probabilistic power spectral density objects for each trace
all_ppsds = []
all_ppsd_names = []
for stream in streams:
    ppsds = []
    ppsd_names = []
    for tr in stream:
        metadata = Parser(metadata_directory + tr.stats.station + tr.stats.channel[-1:] + '.seed')
        ppsd = PPSD(tr.stats, metadata)
        ppsd.add(tr)
        ppsds.append(ppsd)
        ppsd_names.append(tr.stats.station + '_' + tr.stats.channel + '_PPSD')
    all_ppsds.extend(ppsds)
    all_ppsd_names.extend(ppsd_names)

# Plot PPSD data for each trace in 3 views
for n in range(len(all_ppsds)):
    all_ppsds[n].plot(show_coverage=True,
                      show_noise_models=True,
                      xaxis_frequency=True,
                      cmap=pqlx,
                      filename=all_ppsd_names[n] + '.png',
                      show=False)
    all_ppsds[n].plot_spectrogram(filename=all_ppsd_names[n] + '_spectrogram.png',
                                  show=False)
