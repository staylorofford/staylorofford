"""
Working script for seismic site data quality. Initially, this script will have functionality to:
1. Query data over a given time window for site codes in a csv file, on the vertical channel,
2. Filter data for a given filter(s),
3. Calculate data completeness (% samples) over the time window, and produce average data values using the RMS.
"""

import datetime
from obspy.clients.fdsn import Client
from obspy.clients.fdsn.header import FDSNNoDataException
from obspy.core.stream import Stream
from obspy.core.utcdatetime import UTCDateTime
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def query_fdsn(station, location, channel, starttime, endtime, client="https://service.geonet.org.nz"):
    """
    Query FDSN for timeseries data. Input parameter can be specific or include ? wildcard
    :param: station: site code to get data from (string)
    :param: location: location code to get data from (string)
    :param: channel: channel code to get data from (string)
    :param: starttime: UTC start time of desired data (ISO8601 string)
    :param: endtime: UTC end time of desired data (ISO8601 string)
    :param: client: FDSN webservice to use as the client:
            https://service.geonet.org.nz (archive) or https://service-nrt.geonet.org.nz (last week)
    :return: obspy stream object containing requested data
    """
    client = Client(client)
    stream = client.get_waveforms(network='NZ',
                                  station=station,
                                  channel=channel,
                                  location=location,
                                  starttime=starttime,
                                  endtime=endtime)
    return stream


# Give the path to the site file: a csv containing in each row the station code and location code of each site to use
# in the data quality calculations. The station codes should be in the first column of this file, and the location
# code can be in any other column so long as the file header specifies which by the location of the "Location" string
# in the header.
site_file = '/home/samto/git/staylorofford/site_noise_worktree/network_analysis/data_quality/open_sites.csv'

# Give parameters for which data to query
channel_code = '?HZ'  # Channel code to query, FDSN wildcard functionality is OK, but only one channel should be returned.
starttime = '86400'  # starttime can either be an ISO8601 string, or an integer number of seconds to
                                    # query data for before the endtime.
endtime = 'now'  # endtime can be an ISO8601 string or 'now', if 'now' then endtime
                 # will be set to the current UTC time.

# Set how to filter the data, if at all. Use filter_type=None to negate filtering. Filter types are those in obspy.
filter_type = 'bandpass'
minimum_frequency = 2
maximum_frequency = 15

# Parse parameters
if endtime == 'now':
    endtime = datetime.datetime.utcnow()
    if 'T' not in starttime:
        starttime = endtime - datetime.timedelta(seconds=int(starttime))
elif endtime != 'now' and 'T' not in starttime:
    endtime = datetime.datetime.strptime(endtime, '%Y-%m-%dT%H:%M:%SZ')
    starttime = endtime - datetime.timedelta(seconds=int(starttime))
else:
    endtime = datetime.datetime.strptime(endtime, '%Y-%m-%dT%H:%M:%SZ')
    starttime = datetime.datetime.strptime(starttime, '%Y-%m-%dT%H:%M:%SZ')

# Determine which FDSN client to use
splittime = datetime.datetime.utcnow() - datetime.timedelta(days=7)
if splittime <= starttime and splittime <= endtime:
    # All data is within the nrt-client window
    client = 'https://service-nrt.geonet.org.nz'
elif starttime < splittime and splittime <= endtime:
    # Query starts before nrt client window, but endtime is in nrt client window
    client = 'both'
elif starttime < splittime and endtime < splittime:
    # All data is within the archive-client window
    client = 'https://service.geonet.org.nz'

# Parse site code metadata from file
site_metadata = pd.read_csv(site_file, header=0, index_col=0)

# Query data for each site code
outrows = []
for m in range(len(site_metadata)):
    # Get site information
    site_code = site_metadata.index[m]
    location_code = site_metadata.iloc[m]['Location']

    # Query data from FDSN
    if client == 'both':
        data_arc, data_nrt = None, None
        try:
            data_arc = query_fdsn(station=site_code,
                                  location=location_code,
                                  channel=channel_code,
                                  starttime=UTCDateTime(starttime),
                                  endtime=UTCDateTime(splittime),
                                  client='https://service.geonet.org.nz')
        except FDSNNoDataException:
            print('No data exists for: ')
            print('Site code: ' + site_code)
            print('Location code: ' + location_code)
            print('Channel code: ' + channel_code)
            print('Starttime: ' + starttime.isoformat())
            print('Endtime: ' + splittime.isoformat())
            print('')
        try:
            data_nrt = query_fdsn(station=site_code,
                                  location=location_code,
                                  channel=channel_code,
                                  starttime=UTCDateTime(splittime),
                                  endtime=UTCDateTime(endtime),
                                  client='https://serive-nrt.geonet.org.nz')
        except FDSNNoDataException:
            print('No data exists for: ')
            print('Site code: ' + site_code)
            print('Location code: ' + location_code)
            print('Channel code: ' + channel_code)
            print('Starttime: ' + splittime.isoformat())
            print('Endtime: ' + endtime.isoformat())
            print('')

        # Combine data
        if data_arc is not None and data_nrt is not None:
            data = data_arc + data_nrt
        elif data_arc is not None and data_nrt is None:
            data = data_arc
        elif data_arc is None and data_nrt is not None:
            data = data_nrt
        else:
            data = Stream()
    else:
        try:
            data = query_fdsn(station=site_code,
                              location=location_code,
                              channel=channel_code,
                              starttime=UTCDateTime(starttime),
                              endtime=UTCDateTime(endtime),
                              client=client)
        except FDSNNoDataException:
            print('No data exists for: ')
            print('Site code: ' + site_code)
            print('Location code: ' + location_code)
            print('Channel code: ' + channel_code)
            print('Starttime: ' + starttime.isoformat())
            print('Endtime: ' + endtime.isoformat())
            print('')
            data = Stream()

    if len(data) == 0:
        continue

    # Subtract a running mean of 10 second duration from the data to remove the influence of long periods
    if (data[0].stats.endtime - data[0].stats.starttime) < 10:
        running_mean = np.mean(data[0].data)
        for n in range(len(data[0].data)):
            data[0].data[n] -= running_mean
    else:
        running_mean = []
        for n in range(5 * int(data[0].stats.sampling_rate),
                       len(data[0].data) - 5 * int(data[0].stats.sampling_rate)):
            running_mean.append(np.mean(data[0].data[n - 5 * int(data[0].stats.sampling_rate):
                                                     n + 5 * int(data[0].stats.sampling_rate)]))
        for n in range(5 * int(data[0].stats.sampling_rate),
                       len(data[0].data) - 5 * int(data[0].stats.sampling_rate)):
            data[0].data[n] -= running_mean[n - 5 * int(data[0].stats.sampling_rate)]
        # Remove data that was not adjusted by the running mean subtraction
        data[0].data = data[0].data[5 * int(data[0].stats.sampling_rate):
                                    len(data[0].data) - 5 * int(data[0].stats.sampling_rate)]

    # Apply a filter to the data
    if filter_type:
        data = data.filter(type=filter_type,
                           freqmin=minimum_frequency,
                           freqmax=maximum_frequency)

    # Calculate the RMS value for the data
    RMS = 0
    actual_data = 0
    mean = np.mean(data[0].data)
    for n in range(len(data[0].data)):
        RMS += math.pow((mean - data[0].data[n]), 2)
    RMS = math.sqrt(RMS / len(data[0].data))

    # Calculate data completeness
    expected_data = (data[0].stats.endtime - data[0].stats.starttime) * data[0].stats.sampling_rate
    actual_data = len(data[0].data)

    # Save results

    outrows.append(site_code + ',' +
                   location_code + ',' +
                   channel_code + ',' +
                   starttime.isoformat() + 'Z,' +
                   endtime.isoformat() + 'Z,' +
                   str(100 * expected_data / actual_data)[:6] + ',' +
                   str(RMS)[:6] + '\n')

with open('seismic_site_data_quality_output.csv', 'w') as outfile:
    outfile.write('Station,Location,Channel,Start Time,End Time,Completeness,RMS\n')
with open('seismic_site_data_quality_output.csv', 'a') as outfile:
    for row in outrows:
        outfile.write(row)
