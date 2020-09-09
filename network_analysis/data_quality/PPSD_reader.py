"""
Open obspy PPSDs from .npy files and produce plots for the operator to inspect
Intention is to integrate this script and seismic_site_data_quality.py once both are mature
"""

from obspy.clients.fdsn import Client
from obspy.signal import PPSD
from obspy.core.utcdatetime import UTCDateTime
from glob import glob
import numpy as np


def get_metadata(filename):

    """
    Get a site's metadata for PPSD processing from the PPSD filename
    """

    client = Client("https://service.geonet.org.nz")
    metadata = client.get_stations(network='NZ',
                                   station=filename.split('_')[0],
                                   location=filename.split('_')[1],
                                   channel=filename.split('_')[2],
                                   starttime=UTCDateTime(filename.split('_')[3]),
                                   endtime=UTCDateTime(filename.split('_')[3]) + 86399,  # Assume day-long PPSD
                                   level='response')

    return metadata

# Set script parameters

file_directory = '/home/samto/git/staylorofford/network_analysis/data_quality/'

# Collect PPSD file paths and prepare output lists
filepaths = glob(file_directory + '*npz')
sites, locations, channels = [], [], []
for filepath in filepaths:
    filename = filepath.split('/')[-1]
    sites.append(filename.split('_')[0])
    locations.append(filename.split('_')[1])
    channels.append(filename.split('_')[2])
sites, locations, channels = zip(*sorted(zip(sites, locations, channels)))  # Order sites

# Order lists for data parsing
unique_sites = list(set(sites))  # Get unique list of sites
unique_site_locations = []
unique_site_channel_locations = []
for n, site in enumerate(unique_sites):
    site_locations = []
    for m, location in enumerate(locations):
        if sites[m] == site:
            site_locations.append(location)  # Collect all location codes for the given site
    site_locations = list(set(site_locations))  # Get unique list of location codes for the given site
    unique_site_locations.append(site_locations)  # Save to unique list

    site_location_channels = [[] for m in range(len(site_locations))]
    for m, location in enumerate(locations):
        for l, channel in enumerate(channels):
            if sites[l] == site and locations[l] == location:
                if channel not in site_location_channels[site_locations.index(location)]:
                    site_location_channels[site_locations.index(location)].append(channel)
    unique_site_channel_locations.append(site_location_channels)

# Parse PPSD data
parsed_files = []
for n, site in enumerate(unique_sites):
    for m, location in enumerate(unique_site_locations[n]):
        for l, channel in enumerate(unique_site_channel_locations[n][m]):
            # Initialise PPSD
            num = 0
            for filepath in filepaths:
                filename = filepath.split('/')[-1]
                if (filename.split('_')[0] == site and
                        filename.split('_')[1] == location and
                        filename.split('_')[2] == channel and
                        filepath not in parsed_files):
                    parsed_files.append(filepath)
                    metadata = get_metadata(filename)
                    if num == 0:
                        all_ppsd = PPSD.load_npz(filepath, metadata)
                        min_td = all_ppsd._times_data[0][0]
                        max_td = all_ppsd._times_data[0][1]
                        num += 1
                    else:
                        ppsd = PPSD.load_npz(filepath, metadata)
                        # Merge PPSDs across all time
                        # To do this, extend one PPSD's data lists by those of all others
                        # (this is done by tweaking attribute values of obspy PPSD objects - a tricky business!)
                        # NOTE: assumption is that processing parameters are the same between PPSDs!
                        all_ppsd._times_processed.extend(ppsd._times_processed)
                        all_ppsd._binned_psds.extend(ppsd._binned_psds)
                        all_ppsd._times_gaps.extend(ppsd._times_gaps)
                        if ppsd._times_data[0][0] < min_td:
                            min_td = ppsd._times_data[0][0]
                        if ppsd._times_data[0][1] > max_td:
                            max_td = ppsd._times_data[0][1]
            # Update time window of data
            all_ppsd._times_data = [[min_td, max_td]]
            # Sort data by time
            all_ppsd._times_processed, all_ppsd._binned_psds = zip(*sorted(zip(all_ppsd._times_processed, all_ppsd._binned_psds)))
            # Plot PPSD data
            all_ppsd.plot()
