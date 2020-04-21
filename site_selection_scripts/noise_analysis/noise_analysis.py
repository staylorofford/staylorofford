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

day_file_directory_root = '/mnt/hgfs/VMSHARE/KUZ testing/'

# Define sensor + datalogger metadata (dataless SEED) directory

metadata_directory = '/home/samto/git/staylorofford/site_selection/site_selection_scripts/metadata/'

# Find all miniseed files under the root data directory into lists split by site
print('Finding all miniSEED files in the root directory...')
sites = []
channels = []
all_files = []
for root, dirs, files in os.walk(day_file_directory_root):
    for file in files:
        if file[-3:] == '.ms':
            sites.append(file.split('.')[0])
            channels.append(file.split('.')[2][-3:])
            all_files.append(root + '/' + file)
sites = list(set(sites))  # Define site list
channels = list(set(channels))  # Define channel list
split_files = [[[] for n in range(len(channels))] for m in range(len(sites))]
for file in all_files:
    print(file)
    site_idx = sites.index(file.split('/')[-1].split('.')[0])
    channel_idx = channels.index(file.split('/')[-1].split('.')[2][-3:])
    split_files[site_idx][channel_idx].append(file)

# Load in data for each site and perform noise analysis
for m, site in enumerate(sites):
    print('Parsing data for site: ' + site)
    for n, channel in enumerate(channels):
        streams = Stream()
        for o in range(len(split_files[m][n])):
            print('Parsing file:')
            print(split_files[m][n][o])
            streams += obspy.read(split_files[m][n][o])
        print('Merging streams...')
        streams.merge()
        print('Current data is:')
        print(streams)

        # Build probabilistic power spectral density objects for each trace
        all_ppsds = []
        all_ppsd_names = []
        for stream in streams:
            print('Calculating PPSDs for stream:')
            print(stream)
            ppsds = []
            ppsd_names = []
            metadata = Parser(metadata_directory + stream.stats.station + stream.stats.channel[-1:] + '.seed')
            ppsd = PPSD(stream.stats, metadata)
            ppsd.add(stream)
            ppsds.append(ppsd)
            ppsd_names.append(stream.stats.station + '_' + stream.stats.channel + '_PPSD')
        all_ppsds.extend(ppsds)
        all_ppsd_names.extend(ppsd_names)

        # Plot PPSD data for each trace in 3 views
        print('Plotting PPSD data...')
        for n in range(len(all_ppsds)):
            all_ppsds[n].plot(show_coverage=True,
                              show_noise_models=True,
                              xaxis_frequency=True,
                              cmap=pqlx,
                              filename=all_ppsd_names[n] + '.png',
                              show=False)
            all_ppsds[n].plot_temporal(filename=all_ppsd_names[n] + '_temporal_fbands.png',
                                       show=False)
            all_ppsds[n].plot_spectrogram(filename=all_ppsd_names[n] + '_spectrogram.png',
                                          show=False)
            all_ppsds[n].save_npz(all_ppsd_names[n])
