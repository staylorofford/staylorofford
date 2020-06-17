# -*- coding: utf-8 -*-
"""
Combine local data from test sites with data from a reference site (queried via FDSN)
and perform analysis and plotting on it. 
Code is derived from timeseries_analysis.py and earthquake_location.py by Sam Taylor-Offord.
"""

import datetime
import io
import obspy
from obspy.clients.fdsn import Client
from obspy.io.quakeml.core import Unpickler
import os
import matplotlib.pyplot as plt
import numpy as np
import pycurl

# Set up objects to use imported modules
quakeml_reader = Unpickler()


def curl(curlstr):
    """
    Perform curl with curlstr
    :param curlstr: string to curl
    :return: curl output
    """

    buffer = io.BytesIO()
    c = pycurl.Curl()
    c.setopt(c.URL, curlstr)
    c.setopt(c.WRITEDATA, buffer)
    c.perform()
    c.close()

    return buffer.getvalue()


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
    client = Client("https://service.geonet.org.nz")
    stream = client.get_waveforms(network='NZ',
                                  station=station,
                                  channel=channel,
                                  location=location,
                                  starttime=starttime,
                                  endtime=endtime)
    return stream


def FDSN_event_query(eventID):
    """
    Use obspy with pycurl to query event details via FDSN..

    :param eventID: FDSN earthquake eventID
    :return: obspy object containing earthquake details
    """

    query = 'https://service.geonet.org.nz/fdsnws/event/1/query?eventid=' + eventID
    queryresult = curl(query)
    event_details = quakeml_reader.loads(queryresult)

    return event_details


# Give eventID to analyse
eventID = '2020p391429'

# Set parameters for which data to collect from FDSN
stations = ['INSS']
locations = ['20']
channels = ['HNZ']  # the same channels will be loaded for test sites

# Set parameters for where to find test site data
test_data_dir = '/home/samto/PROCESSING/2020/'  # root directory under which all test data exists
test_sites = ['TES1', 'TES2', 'TES3', 'TES4']
location_codes = ['20', '20', '20', '20']

# Set how to filter the data, if at all. Use filter_type=None to negate filtering. Filter types are those in obspy.
filter_type = 'bandpass'
minimum_frequency = 2
maximum_frequency = 15

# Set whether to plot spectrograms instead of seismograms. Use None to negative spectrogram plotting.
# Note: if both spectrogram and normalise are true, then time domain waveforms will be plotted over spectrograms.
spectrogram = False

# Set whether to normalise each trace. Use None or False to negate normalisation.
# Note: if both spectrogram and normalise are true, then time domain waveforms will be plotted over spectrograms.
normalise = False

# Reference times to overlay on plots. Use None to negate reference time plotting. Time format is ISO8601.
reference_times = None

# Query event start time from FDSN
start_time = FDSN_event_query(eventID)[0].origins[0].time + 0
end_time = start_time + 50

# Format reference time
if reference_times:
    for n in range(len(reference_times)):
        reference_times[n] = datetime.datetime.strptime(reference_times[n],
                                                        '%Y-%m-%dT%H:%M:%SZ')

# Query data from FDSN for the reference site(s)
stream = obspy.core.stream.Stream()
for station in stations:
    for location in locations:
        for channel in channels:
            query_stream = query_fdsn(station,
                                      location,
                                      channel,
                                      start_time,
                                      end_time)
        stream += query_stream

# Find which days to load data for
start_date = start_time.datetime.timetuple().tm_yday
end_date = end_time.datetime.timetuple().tm_yday
if start_date != end_date:
    data_doys = [start_date, end_date]
else:
    data_doys = [start_date]

# Parse local data over the period of the event
test_data = obspy.core.stream.Stream()
for root, dirs, files in os.walk(test_data_dir):
    for file in files:
        for test_site in test_sites:
            for channel in channels:
                for doy in data_doys:
                    if test_site in file and channel in file and str(doy) in file:
                        loaded_data = obspy.read(root + '/' + file)
                        if loaded_data[0].stats.starttime <= start_time and loaded_data[0].stats.endtime >= end_time:
                            loaded_data.merge(method=1,
                                              fill_value='interpolate')  # Ensure loaded data is in a single stream
                            test_data += loaded_data

# Trim local data to event period
test_data.trim(starttime=start_time,
               endtime=end_time)

# Combine reference and local data
stream += test_data

# Filter the data as desired
if filter_type:
    stream = stream.detrend(type='simple')
    stream = stream.filter(type=filter_type,
                           freqmin=minimum_frequency,
                           freqmax=maximum_frequency,
                           corners=4,
                           zerophase=True)
    stream.trim(starttime=start_time,
                endtime=end_time)

print(stream)

# Find min/max values of each trace
min_max_values = [[0] * len(stream),
                  [0] * len(stream),
                  [0] * len(stream),
                  [0] * len(stream),
                  [''] * len(stream)]
for n in range(len(stream)):
    if normalise and not spectrogram:
        stream[n].data = stream[n].data / max(abs(min(stream[n].data)),
                                              max(stream[n].data))
    if normalise and spectrogram:
        # First shift the data so the mean sits at the middle of the nyquist frequency range
        mean_value = np.mean(stream[n].data)
        shift_value = mean_value - 25
        stream[n].data = stream[n].data - shift_value

        # Then find the max-min range in the normalised data and scale the data to the nyquist frequency range
        max_min_range = max(stream[n].data) - min(stream[n].data)
        scaling_factor = 50 / max_min_range
        stream[n].data = stream[n].data * scaling_factor

        # Then find how much the data needs to be shifted by to sit at the centre of the nyquist frequency range
        shift_value = max(50 - max(stream[n].data), 0 - min(stream[n].data))
        stream[n].data = stream[n].data + shift_value
    min_max_values[0][n] = stream[n].stats.starttime
    min_max_values[1][n] = stream[n].stats.endtime
    min_max_values[2][n] = min(stream[n].data)
    min_max_values[3][n] = max(stream[n].data)
    min_max_values[4][n] = stream[n].stats.network + '.' + stream[n].stats.station + '.' + stream[
        n].stats.location + '.' + \
                           stream[n].stats.channel
min_plot_time = min(min_max_values[0])
max_plot_time = max(min_max_values[1])

# Pad out traces to align all on minimum time
stream.trim(starttime=min_plot_time,
            endtime=max_plot_time,
            pad=True,
            fill_value=None)

# Update reference times to be relative to minimum data start time
if reference_times:
    for m in range(len(reference_times)):
        reference_times[m] = (reference_times[m] - min_plot_time.datetime).total_seconds()

# Plot the data and overlay reference times if desired
fig = plt.figure(figsize=(12, 9))
stream.plot(fig=fig,
            type='relative',
            zorder=2)
for n in range(len(fig.axes)):
    axis = fig.axes[n]
    plotting_stream = axis.texts[0]._text
    if not spectrogram:
        min_max_index = min_max_values[4].index(plotting_stream)
        axis.set_ylim(min_max_values[2][min_max_index] - min_max_values[2][min_max_index] / 10,
                      min_max_values[3][min_max_index] + min_max_values[3][min_max_index] / 10)
        axis.locator_params(axis='y',
                            nbins=3)
        if reference_times:
            axis.vlines(x=reference_times,
                        ymin=min_max_values[2][min_max_index],
                        ymax=min_max_values[3][min_max_index],
                        linestyles='dashed',
                        colors='red')
    else:
        stream[n].spectrogram(axes=axis,
                              zorder=1)
        axis.set_ylim(0, 50)
        axis.locator_params(axis='y',
                            nbins=3)
        if reference_times:
            axis.vlines(x=reference_times,
                        ymin=0,
                        ymax=50,
                        linestyles='dashed',
                        colors='red')
if not filter_type:
    filter_type = 'no'
plt.savefig(
    eventID + '_' + filter_type + '_filter_' + str(minimum_frequency) + '-' + str(maximum_frequency) + '_Hz' + '.png',
    dpi=300)
plt.show()