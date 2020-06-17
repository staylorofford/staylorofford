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

day_file_directory_root = '/home/samto/PROCESSING/2020/'
parse_channels = ['HNZ', 'HHZ']  # Which channels to parse: you must have response data for each channel.

# Define sensor + datalogger metadata (dataless SEED) directory

metadata_directory = '/home/samto/PROCESSING/2020/NZ/METADATA/'

# Find all miniseed files under the root data directory into lists split by site. Note: one channel per file only!
print('Finding all miniSEED files in the root directory...')
sites = []
location_codes = []
channels = []
all_files = []
for root, dirs, files in os.walk(day_file_directory_root):
    for file in files:
        # See if it's readable
        try:
            stream = obspy.read(root + '/' + file)
        except:
            continue

        if stream[0].stats.channel in parse_channels:
            print('File ' + file + ' found')
            sites.append(stream[0].stats.station)
            location_codes.append(stream[0].stats.location)
            channels.append(stream[0].stats.channel)
            all_files.append(root + '/' + file)

print('Indexing files for processing...')
sites = list(set(sites))  # Define site list
location_codes = list(set(location_codes))  # Define location codes
channels = list(set(channels))  # Define channel list
split_files = [[[[] for o in range(len(channels))] for n in range(len(location_codes))] for m in range(len(sites))]
for file in all_files:
    stream = obspy.read(file)
    site_idx = sites.index(stream[0].stats.station)
    loc_idx = location_codes.index(stream[0].stats.location)
    channel_idx = channels.index(stream[0].stats.channel)
    split_files[site_idx][loc_idx][channel_idx].append(file)

# Load in data for each site and perform noise analysis
for m, site in enumerate(sites):
    print('Parsing data for site: ' + site)
    for n, loc in enumerate(location_codes):
        print('Parsing data for loc: ' + loc)
        for o, channel in enumerate(channels):
            print('Parsing data for channel: ' + channel)
            streams = Stream()
            if len(split_files[m][n][o]) == 0:
                print('No data for this combination')
                continue
            for p in range(len(split_files[m][n][o])):
                print('Parsing file:')
                print(split_files[m][n][o][p])
                streams += obspy.read(split_files[m][n][o][p])
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

                try:
                    metadata = Parser(metadata_directory + stream.stats.station + '_' + stream.stats.location + '.txt')
                except OSError:
                    print('No metadata for this site')
                    continue

                ppsd = PPSD(stream.stats, metadata)
                ppsd.add(stream)
                ppsds.append(ppsd)
                ppsd_names.append(stream.stats.station + '_' + stream.stats.location + '_' + stream.stats.channel + '_PPSD')
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
                # all_ppsds[n].plot_temporal(filename=all_ppsd_names[n] + '_temporal_fbands.png',
                #                            show=False, period=1)
                all_ppsds[n].plot_spectrogram(filename=all_ppsd_names[n] + '_spectrogram.png',
                                              show=False)
                all_ppsds[n].save_npz(all_ppsd_names[n])
