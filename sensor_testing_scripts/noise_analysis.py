# -*- coding: utf-8 -*-
"""
Combine local data from test sites [to do: with data from a reference site (queried via FDSN)]
and perform noise analysis and plotting on it.
"""

import matplotlib.pyplot as plt
import os
import obspy
from obspy.core.stream import Stream
from obspy.imaging.cm import pqlx
from obspy.io.xseed import Parser
from obspy.signal import PPSD

# Define root data directory

day_file_directory_root = '/scratch/SDRFP/archive_150720'

# Define sensor + datalogger metadata (dataless SEED) directory

metadata_directory = '/scratch/SDRFP/METADATA/'

# Find all miniseed files under the root data directory into lists split by site
print('Finding all miniSEED files in the root directory...')
sites = []
channels = []
locations = []
all_files = []
for root, dirs, files in os.walk(day_file_directory_root):
    for file in files:
        try:
            st = obspy.read(root + '/' + file)
            if 'HHZ' in st[0].stats.channel or 'HNZ' in st[0].stats.channel:
                sites.append(st[0].stats.station)
                locations.append(st[0].stats.location)
                channels.append(st[0].stats.channel)
                all_files.append(root + '/' + file)
        except:
            continue
sites = list(set(sites))  # Define site list
channels = list(set(channels))  # Define channel list
locations = list(set(locations))  # " location
split_files = [[[[] for o in range(len(locations))] for n in range(len(channels))] for m in range(len(sites))]
for file in all_files:
    st = obspy.read(file)
    site_idx = sites.index(st[0].stats.station)
    channel_idx = channels.index(st[0].stats.channel)
    location_idx = locations.index(st[0].stats.location)
    split_files[site_idx][channel_idx][location_idx].append(file)

# Load in data for each site and perform noise analysis
for m, site in enumerate(sites):
    print('Parsing data for site: ' + site)
    for n, channel in enumerate(channels):
        for o, location in enumerate(locations):
            streams = Stream()
            for p in range(len(split_files[m][n][o])):
                # Save daily data as dayplot
                print('Parsing file:')
                print(split_files[m][n][o][p])

                day_stream = obspy.read(split_files[m][n][o][p])
                day_stream.plot(type='dayplot', size=(1200, 800),
                                outfile=site + '_' + channel + '_' + day_stream[0].stats.starttime.isoformat()[:10])

                streams += day_stream
            if len(streams) == 0:
                # Skip if there is no data for the given station/channel/location combination
                continue
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
                metadata_file = metadata_directory + \
                                stream.stats.station + '_' + \
                                stream.stats.location + '_' + \
                                stream.stats.channel + '.seed'
                if not os.path.exists(metadata_file):
                    print('Could not find metadata as file ' + metadata_file + ' will not produce PPSDs for stream')
                    continue
                metadata = Parser(metadata_file)
                if 'N' in stream.stats.channel[1:2]:
                    ppsd = PPSD(stream.stats, metadata, special_handling='hydrophone', db_bins=(-200, 20, 1))
                if 'H' in stream.stats.channel[1:2]:
                    ppsd = PPSD(stream.stats, metadata, db_bins=(-200, 20, 1))

                ppsd.add(stream)
                ppsds.append(ppsd)
                ppsd_names.append(stream.stats.station + '_' + stream.stats.channel + '_PPSD')
            all_ppsds.extend(ppsds)
            all_ppsd_names.extend(ppsd_names)

            # Plot PPSD data for each trace
            print('Plotting PPSD data if it exists...')
            for n in range(len(all_ppsds)):
                all_ppsds[n].plot(show_coverage=True,
                                  show_noise_models=True,
                                  xaxis_frequency=True,
                                  cmap=pqlx,
                                  show=False)
                fig = plt.gcf()
                ax = fig.axes[0]
                ax.set_ylim(-200, 20)
                plt.savefig(all_ppsd_names[n] + '.png', dpi=300, format='png')
                all_ppsds[n].plot_spectrogram(filename=all_ppsd_names[n] + '_spectrogram.png',
                                              show=False)
                all_ppsds[n].save_npz(all_ppsd_names[n])
