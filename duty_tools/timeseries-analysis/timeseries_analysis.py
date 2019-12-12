# -*- coding: utf-8 -*-
"""
Query timeseries data from FDSN and perform analysis and plotting on it.
"""

import datetime
import obspy
from obspy.clients.fdsn import Client
import matplotlib.pyplot as plt


def query_fdsn(station, location, channel, starttime, endtime):

    """
    Query FDSN for timeseries data. Input parameter can be specific or include ? wildcard
    :param: station: site code to get data from (string)
    :param: location: location code to get data from (string)
    :param: channel: channel code to get data from (string)
    :param: starttime: UTC start time of desired data (ISO8601 string)
    :param: endtime: UTC end time of desired data (ISO8601 string)
    :return: obspy stream object containing requested data
    """

    client = Client("https://service-nrt.geonet.org.nz")
    stream = client.get_waveforms(network='NZ',
                                  station=station,
                                  channel=channel,
                                  location=location,
                                  starttime=starttime,
                                  endtime=endtime)
    return stream


# Set parameters for which data to analyse using lists of strings

stations = ['WIZ',
            'WSRZ']
locations = ['??']
channels = ['H??']
starttime = '2019-12-09T01:09:47.000'
endtime = '2019-12-09T01:16:47.000'

# Set how to filter the data, if at all. Use filter_type=None to negate filtering. Filter types are those in obspy.

filter_type = None
minimum_frequency = 0.1
maximum_frequency = 1

# Set whether to normalise each trace. Use None or False to negate normalisation.

normalise = False

# Reference times to overlay on plots. Use None to negate reference time plotting. Time format is ISO8601.

reference_times = ['2019-12-09T01:11:47Z']

if reference_times:
    for n in range(len(reference_times)):
        reference_times[n] = datetime.datetime.strptime(reference_times[n],
                                                        '%Y-%m-%dT%H:%M:%SZ')

# Query the data from FDSN

stream = obspy.core.stream.Stream()
for station in stations:
    for location in locations:
        for channel in channels:
            query_stream = query_fdsn(station,
                                      location,
                                      channel,
                                      starttime,
                                      endtime)
        stream += query_stream

# Filter the data as desired

if filter_type:
    stream = stream.filter(type=filter_type,
                           freqmin=minimum_frequency,
                           freqmax=maximum_frequency)

# Find min/max values of each trace

min_max_values = [[0] * len(stream),
                  [0] * len(stream),
                  [''] * len(stream)]
for n in range(len(stream)):
    if normalise:
        stream[n].data = stream[n].data / max(abs(min(stream[n].data)),
                                              max(stream[n].data))
    min_max_values[0][n] = min(stream[n].data)
    min_max_values[1][n] = max(stream[n].data)
    min_max_values[2][n] = 'NZ.' + stream[n].stats.station + '.' + stream[n].stats.location + '.' + \
                           stream[n].stats.channel

# Plot the data and overlay reference times if desired

fig = plt.figure()
stream.plot(fig=fig)
for n in range(len(fig.axes)):
    axis = fig.axes[n]
    plotting_stream = axis.texts[0]._text
    min_max_index = min_max_values[2].index(plotting_stream)
    axis.set_ylim(min_max_values[0][min_max_index],
                  min_max_values[1][min_max_index])
    axis.locator_params(axis='y',
                        nbins=3)
    if reference_times:
        axis.vlines(x=reference_times,
                    ymin=min_max_values[0][min_max_index],
                    ymax=min_max_values[1][min_max_index],
                    linestyles='dashed',
                    colors='red')
plt.show()



