"""
Script to find a sensor's orientation using teleseismic events. Requires events be queryable from FDSN, and that
stations to orient have metadata and data in FDSN for the events.

Note: functionality with multiple stations to orient has not been tested.
"""

import datetime
import glob
import fnmatch
import io
import math
import matplotlib.pyplot as plt
import multiprocessing
import numpy as np
import numpy.ma as ma
import obspy
from obspy.clients.fdsn import Client
from obspy.taup import TauPyModel
from obspy.io.quakeml.core import Unpickler
from obspy.core.utcdatetime import UTCDateTime
import os
import pycurl
import time
import xml.etree.ElementTree as ET

quakeml_reader = Unpickler()
spherical_velocity_model = TauPyModel(model="iasp91")


def parse_csv(csv_path, header=True):

    """
    Parse the contents of CSV file at csv_path into nested lists
    :param csv_path: full path to CSV file
    :return: list of CSV header columns, nested lists of CSV row data
    """

    data = []
    with open(csv_path, 'r') as openfile:
        if header:
            rc = -1
            for row in openfile:
                if rc == -1:
                    header = row[:-1].split(',')
                    rc = 0
                else:
                    data.append(row[:-1].split(','))
            return header, data
        else:
            for row in openfile:
                data.append(row[:-1].split(','))
            return data


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


def get_event(eventID):

    """
    Get a single event from the GeoNet FDSN service
    :param eventID: GeoNet eventID
    :return: FDSN event
    """

    query = "https://service.geonet.org.nz/fdsnws/event/1/query?eventid=" + eventID
    queryresult = curl(query)
    event = quakeml_reader.loads(queryresult)
    return event


def query_fdsn(station, location, channel, starttime, endtime):

    """
    Query FDSN for timeseries data. Input parameter can be specific or include ? wildcard.
    Note: only works with archive data, not near-real-time.
    Note: instrument response is included with data.
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
                                  endtime=endtime,
                                  attach_response=True),  # For reasons unknown, this comma is very important!
    return stream


def find_lag_time(stream, reference_stream, phase):

    """
    Find the lag time between a given stream and a reference stream using cross-correlation of the
    total energy, which is independent of component alignment.
    :param stream: obspy stream object of the seismogram to calculate lag time for
    :param reference_stream: object stream object of the seismogram to use as reference for lag time calculation
    :param phase: whether to use total energy on vertical channel (P phase) or horizontal channels (S phase)
    :return: lag time between the two sensors in seconds relative to the reference stream
    """

    # Create normalised amplitude envelopes of the data using horizontal total energy
    stream_envelope = calculate_total_energy(stream, phase)
    max_se = ma.max(stream_envelope)
    stream_envelope /= max_se
    ref_envelope = calculate_total_energy(reference_stream, phase)
    max_re = ma.max(ref_envelope)
    ref_envelope /= max_re

    # Ensure all data are masked arrays
    if not ma.is_masked(stream_envelope):
        stream_envelope = ma.masked_array(stream_envelope)
    if not ma.is_masked(ref_envelope):
        ref_envelope = ma.masked_array(ref_envelope)

    # Find the lag time from the maximum cross-correlation value between the two waveforms
    xcorr_values = []
    ref_envelope = ref_envelope.filled(0).tolist() + len(stream_envelope) * [0]
    ref_envelope = np.asarray(
        smooth_data(ref_envelope, int(round(1/float(values[parameters.index('upper_frequency')])))))
    for m in range(2 * len(stream_envelope)):

        # Shift the stream
        if m <= len(stream_envelope):
            shifted_stream_envelope = (len(stream_envelope) - m) * [0] + stream_envelope.filled(0).tolist() + m * [0]
        else:
            shifted_stream_envelope = max(0, len(stream_envelope) - m) * [0] + \
                                      stream_envelope[:len(stream_envelope) - m].filled(0).tolist() + m * [0]
        shifted_stream_envelope = np.asarray(
            smooth_data(shifted_stream_envelope, int(round(1/float(values[parameters.index('upper_frequency')])))))

        # Perform cross-correlation
        try:
            xcorr_value = np.corrcoef(shifted_stream_envelope, ref_envelope)[0][1]
        except ValueError:
            print('Cross-correlation failed! Perhaps one stream is a data point different to the other? '
                  'This can occur for certain corner frequency and data length combinations... It is a bug.')
        xcorr_values.append(xcorr_value)

    # Find lag time from highest cross-correlation value
    max_xcorr_value = max(xcorr_values)
    print('Correlation value at best alignment of total energy waveforms is: ' + str(max_xcorr_value))
    lag_time = (1 / stream[0].stats.sampling_rate) * (len(stream_envelope) - xcorr_values.index(max_xcorr_value))

    return lag_time


def calculate_total_energy(stream, phase):

    """
    Produce the total energy trace for a given 3-component seismogram.
    :param stream: obspy stream object containing 3-component seismogram to calculate total energy for
    :return: list containing total energy waveform
    """

    # Remove undesired data for the given phase
    for m in range(len(stream)):
        if stream[m].stats.channel[-1] == 'Z':
            if phase.lower() == 'p':
                stream = [stream[m]]
            elif phase.lower() == 's':
                stream.pop(m)

    # Calculate total energy waveform
    te_waveform = []
    for m in range(len(stream[0].data)):
        try:
            te_waveform.append(stream[0].data[m] ** 2 + stream[1].data[m] ** 2)
        except:
            # Catch IndexError and TypeError
            try:
                te_waveform.append(stream[0].data[m] ** 2)
            except TypeError:
                # Finally, catch TypeError if data is a NaN
                te_waveform.append(float('nan'))

    return te_waveform


def smooth_data(data, N):

    """
    Apply a running mean to data using the N samples around a given point to produce mean values.
    :param data: list of data to apply running mean to
    :param N: number of samples either side of each data point to use in the mean calculation
    :return: smoothed data in list
    """

    # Replace nan values in data with 0s
    for m in range(len(data)):
        if np.isnan(data[m]):
            data[m] = 0

    # Smooth data
    smoothed_data = [0] * len(data)
    for m in range(N, len(data) - N):
        smoothed_data[m] = np.nanmean(data[m - N: m + N])

    return smoothed_data


def find_rotation_angle(shifted_seismogram, reference_seismogram, phase):

    """
    Calculate the angle between the horizontal components of each sensor and those of the reference sensor.
    Uses the method of rotating the seismogram horizontal components by a given angle and cross-correlating
    the result against the reference seismogram. The rotation angle that produces the best correlation is
    taken as the rotation angle and returned to the operator.
    :param shifted_seismogram: the seismogram from the station with unknown orientation, filtered, time-shifted, and
           without instrument response.
    :param reference_seismogram: the seismogram from the reference station, filtered and without instrument response.
    :param phase that is being used: if P then only correlate the 1* and -1* shifted stream values, if S do the full
           routine.
    :return: counterclockwise angle between the horizontal components of the sensors at the two stations,
            and list containing cross-correlations values at each angle between 0 and 180.
    """

    normalised_xcorr_values = []
    if phase.lower() == 'p':
        for factor in [1, -1]:
            x_corr_mean = 0
            for n in range(len(reference_seismogram)):
                shifted_seismogram[n].data = factor * shifted_seismogram[n].data
                x_mean = np.nanmean(shifted_seismogram[n].data)
                y_mean = np.nanmean(reference_seismogram[n].data)
                x_var = np.nanvar(shifted_seismogram[n].data)
                y_var = np.nanvar(reference_seismogram[n].data)
                sum = 0
                for m in range(min(len(reference_seismogram[n].data),
                                   len(shifted_seismogram[n].data))):
                    if np.isnan(reference_seismogram[n][m]) or np.isnan(shifted_seismogram[n][m]):
                        # Skip any nan values. These are unlikely, but can get through the window trimming due to
                        # edge cases where downsampled data contain no nan in the window, but "full" data do.
                        continue
                    sum += ((shifted_seismogram[n][m] - x_mean) *
                            (reference_seismogram[n][m] - y_mean))
                normalised_xcorr_value = (1 / min(len(reference_seismogram[n].data),
                                                  len(shifted_seismogram[n])) *
                                          sum / math.sqrt(x_var * y_var))
                x_corr_mean += normalised_xcorr_value
            normalised_xcorr_values.append(x_corr_mean)
        # Give the max corr value
        max_xcorr_value = max(normalised_xcorr_values)
        print('Cross-correlation values for p-phase correlation (1 * data, -1 * data) are: ')
        print(normalised_xcorr_values)
        max_xcorr_value_idx = [1, -1][normalised_xcorr_values.index(max_xcorr_value)]  # The factor to scale by

    elif phase.lower() == 's':
        for angle in range(0, 360):
            rad = math.radians(angle)
            north_component = []  # "north" component after rotation
            east_component = []  # "east" component after rotation

            # Spin the seismogram components clockwise by angle
            for m in range(len(shifted_seismogram[0])):
                north_component.append(math.cos(rad) * shifted_seismogram[0][m] +
                                       math.sin(rad) * shifted_seismogram[1][m])
                east_component.append(-math.sin(rad) * shifted_seismogram[0][m] +
                                      math.cos(rad) * shifted_seismogram[1][m])

            # Make sure channel order in spun seismogram matches that in reference seismogram
            if reference_seismogram[0].stats.channel[-1] == 'N' and reference_seismogram[1].stats.channel[-1] == 'E':
                spun_seismogram = [north_component, east_component]
            elif reference_seismogram[0].stats.channel[-1] == 'E' and reference_seismogram[1].stats.channel[-1] == 'N':
                spun_seismogram = [east_component, north_component]
            else:
                print('Mate, what is that reference seismogram? It has not got N and E channels? How can it be a '
                      'reference seismogram if it is not oriented? I am not even going to try, I am going to exit.')
                exit()

            # Cross-correlate the seismogram with the reference seismogram
            x_corr_mean = 0
            for n in range(len(reference_seismogram)):
                x_mean = np.nanmean(spun_seismogram[n])
                y_mean = np.nanmean(reference_seismogram[n].data)
                x_var = np.nanvar(spun_seismogram[n])
                y_var = np.nanvar(reference_seismogram[n].data)
                sum = 0
                for m in range(min(len(reference_seismogram[n].data),
                                   len(spun_seismogram[n]))):
                    if np.isnan(reference_seismogram[n][m]) or np.isnan(spun_seismogram[n][m]):
                        # Skip any nan values. These are unlikely, but can get through the window trimming due to
                        # edge cases where downsampled data contain no nan in the window, but "full" data do.
                        continue
                    sum += ((spun_seismogram[n][m] - x_mean) *
                            (reference_seismogram[n][m] - y_mean))
                normalised_xcorr_value = (1 / min(len(reference_seismogram[n].data),
                                                  len(spun_seismogram[n])) *
                                          sum / math.sqrt(x_var * y_var))
                x_corr_mean += normalised_xcorr_value
            normalised_xcorr_values.append(x_corr_mean / 2)

        # Give the max corr value
        max_xcorr_value = max(normalised_xcorr_values)
        max_xcorr_value_idx = normalised_xcorr_values.index(max_xcorr_value)

    return max_xcorr_value_idx, normalised_xcorr_values


def FDSN_station_query(station):

    """
    Use pycurl to query station details via FDSN.

    :param station: station site code
    :return: station latitude (dec. deg.), longitude (dec. deg.), and depth (+ve direction is down, in m)
    """

    # Query station information
    query = 'https://service.geonet.org.nz/fdsnws/station/1/query?station=' + station
    queryresult = curl(query)

    # Parse station information
    root = ET.fromstring(queryresult.decode('utf-8'))
    latitude = float(root[4][3][2].text)
    longitude = float(root[4][3][3].text)
    depth = -1 * float(root[4][3][4].text)

    return latitude, longitude, depth


def get_data(event, arrival_times, query_stations, phase, parameters, values):

    """
    Get data for the given stations and input parameters for the given phase in arrival times.
    """

    # Assume non-scattered phase arrivals do not occur after the phase window length after the
    # last first arrival at any station and vice versa for the early arrivals

    start_time = min(arrival_times) - datetime.timedelta(seconds=
                                                         int(values[parameters.index(phase + '_window')]))
    end_time = max(arrival_times) + datetime.timedelta(seconds=
                                                       int(values[parameters.index(phase + '_window')]))

    # Query waveform data from GeoNet FDSN for each station between the start and end times defined for the event
    stations = []
    streams = []
    downsampled_streams = []
    min_times, max_times = [], []
    for station in query_stations:
        stations.append(station)
        # If the reference data is local, get this data from local storage
        if station == values[parameters.index('reference_station')] and \
                values[parameters.index('local_reference')] == 'True':
            reference_station_files = glob.glob(values[parameters.index('reference_data_dir')] + '*' +
                                                values[parameters.index('reference_station')] + '*')
            # Open all reference station files in turn and note those which contain the event data
            files_with_event = []
            for file in reference_station_files:
                stream = obspy.read(file)
                if stream[0].stats.starttime <= start_time or stream[0].stats.endtime >= end_time:
                    files_with_event.append(file)
            # Load files containing event
            station_stream = obspy.core.stream.Stream()
            for file in files_with_event:
                stream = obspy.read(file,
                                    starttime=UTCDateTime(start_time),
                                    endtime=UTCDateTime(end_time))
                for tr in stream:
                    if fnmatch.fnmatch(tr.stats.channel, values[parameters.index('reference_channels')]) is True:
                        station_stream += tr
            station_stream.merge()
        elif station == values[parameters.index('reference_station')]:
            station_stream = query_fdsn(station,
                                        '??',
                                        values[parameters.index('reference_channels')],
                                        UTCDateTime(start_time),
                                        UTCDateTime(end_time))[0]
        else:
            station_stream = query_fdsn(station,
                                        '??',
                                        values[parameters.index('station_channels')],
                                        UTCDateTime(start_time),
                                        UTCDateTime(end_time))[0]

        # Plot the data before any manipulation
        fig = plt.figure()
        station_stream.plot(fig=fig)
        plt.savefig(station + '_' + event + '_' + phase + '_1-raw-data.png')
        plt.clf()

        # Filter the waveforms of all events to increase waveform similarly at all sensors:
        # Filter corner frequency satisfies the condition that it is much smaller than the lowest seismic velocity
        # in the propagation medium of all events scaled by the linear distance between each pair of sensors.
        station_stream.filter(type='bandpass',
                              freqmin=float(values[parameters.index('lower_frequency')]),
                              freqmax=float(values[parameters.index('upper_frequency')]))
        # Cut first and last 10 seconds of data to remove edge effects introduced by filtering
        station_stream.trim(station_stream[0].stats.starttime + 10,
                            station_stream[0].stats.endtime - 10)

        # Plot the data after filtering
        fig = plt.figure()
        station_stream.plot(fig=fig)
        plt.savefig(station + '_' + event + '_' + phase + '_2-filtered-data.png')
        plt.clf()

        # Resample a copy of the stream to twice the corner frequency
        downsampled_station_stream = station_stream.copy()
        for trace in downsampled_station_stream:
            trace.resample(int(2 * round(1/float(values[parameters.index('upper_frequency')]))))
            min_times.append(trace.stats.starttime)
            max_times.append(trace.stats.endtime)

        streams.append(station_stream)
        downsampled_streams.append(downsampled_station_stream)

    # Pad out traces to ensure an equal number of data in each
    for m in range(len(streams)):

        streams[m].trim(starttime=min(min_times),
                        endtime=max(max_times),
                        pad=True,
                        fill_value=None)
        downsampled_streams[m].trim(starttime=min(min_times),
                                    endtime=max(max_times),
                                    pad=True,
                                    fill_value=None)

    reference_station_idx = stations.index(values[parameters.index('reference_station')])
    stations.pop(reference_station_idx)
    reference_station_stream = streams.pop(reference_station_idx)
    downsampled_rss = downsampled_streams.pop(reference_station_idx)

    return reference_station_stream, downsampled_rss, streams, downsampled_streams


def align_seismograms(stations, arrival_times, streams, downsampled_streams, reference_station_stream, downsampled_rss,
                      phase, parameters, values):

    """
    Align all seismograms for each event to facilitate cross-correlation:
    Use the lag time that produces the maximum cross-correlation value between each sensor and the reference
    sensor's energy traces. In this process, convert all numpy masked arrays to numpy arrays with nan
    values as mask fill values.
    """

    nandices = [0, None]
    shifted_streams = streams
    shifted_downsampled_streams = downsampled_streams
    for m in range(len(streams)):
        # Cut the lag time streams to 5 seconds before and after the arrival at each site
        station_stream = downsampled_streams[m].copy()
        try:
            station_stream.trim(starttime=UTCDateTime(arrival_times[stations.index(stations[m])] -
                                                      datetime.timedelta(seconds=5)),
                                endtime=UTCDateTime(arrival_times[stations.index(stations[m])] +
                                                    datetime.timedelta(seconds=5)))
        except ValueError:
            # When there is no data for the station, use the minimum and maximum from the other sites
            station_stream.trim(starttime=UTCDateTime(min(arrival_times) - datetime.timedelta(seconds=5)),
                                endtime=UTCDateTime(max(arrival_times) + datetime.timedelta(seconds=5)))
        reference_stream = downsampled_rss.copy()
        try:
            reference_stream.trim(starttime=UTCDateTime(arrival_times[stations.index(values[parameters.index(
                'reference_station')])] - datetime.timedelta(seconds=5)),
                                    endtime=UTCDateTime(arrival_times[stations.index(values[parameters.index(
                                        'reference_station')])] + datetime.timedelta(seconds=5)))
        except ValueError:
            # When there is no data for the station, use the minimum and maximum from the other sites
            reference_stream.trim(starttime=UTCDateTime(min(arrival_times) - datetime.timedelta(seconds=5)),
                                  endtime=UTCDateTime(max(arrival_times) + datetime.timedelta(seconds=5)))
        lag_time = find_lag_time(station_stream, reference_stream, phase)
        shift_idx = int(abs(lag_time * streams[m][0].stats.sampling_rate))
        downsampled_shift_idx = int(abs(lag_time * downsampled_streams[m][0].stats.sampling_rate))
        for n in range(len(streams[m])):

            # Ensure all data are masked arrays
            if not ma.is_masked(streams[m][n].data):
                streams[m][n].data = ma.masked_array(streams[m][n].data)
            if not ma.is_masked(downsampled_streams[m][n].data):
                downsampled_streams[m][n].data = ma.masked_array(downsampled_streams[m][n].data)

            # Apply shift
            if lag_time > 0:
                nandices[0] = downsampled_shift_idx
                shifted_streams[m][n].data = np.asarray([float('nan')] * shift_idx +
                                                        streams[m][n].data[:-shift_idx].filled(
                                                            float('nan')).tolist())
                shifted_downsampled_streams[m][n].data = np.asarray([float('nan')] * downsampled_shift_idx +
                                                                    downsampled_streams[m][n].data[
                                                                    :-downsampled_shift_idx].filled(
                                                                        float('nan')).tolist())
            else:
                nandices[1] = len(downsampled_streams[m][n].data) - downsampled_shift_idx
                shifted_streams[m][n].data = np.asarray(streams[m][n].data[shift_idx:].filled(
                    float('nan')).tolist() + [float('nan')] * shift_idx)
                shifted_downsampled_streams[m][n].data = np.asarray(downsampled_streams[m][n].data[
                                                                    downsampled_shift_idx:].filled(
                    float('nan')).tolist() + [float('nan')] * downsampled_shift_idx)

        print(shifted_streams[m][0].stats.station + ' seismograms have been aligned to the reference station by '
                                                    'appling a shift of ' + str(lag_time) + ' seconds')

    for m in range(len(reference_station_stream)):
        if ma.is_masked(reference_station_stream[m].data):
            reference_station_stream[m].data = reference_station_stream[m].data.filled(float('nan'))
        if ma.is_masked(downsampled_rss[m].data):
            downsampled_rss[m].data = downsampled_rss[m].data.filled(float('nan'))

    return nandices, reference_station_stream, downsampled_rss, shifted_streams, shifted_downsampled_streams


def find_xcorr_window(shifted_downsampled_streams, downsampled_rss, nandices, phase, parameters, values):

    """
    Find the time window of each event for which the correlation of the vertical waveforms
    at the sensors is highest: This window is that for which the normalised cross-correlation of the
    total energy traces of each sensor pair is highest.
    """

    xcorr_window = []

    # Calculate horizontal total energy for the reference stream
    reference_total_energy_waveform = calculate_total_energy(downsampled_rss, phase)
    reference_total_energy_waveform = np.asarray(smooth_data(
        reference_total_energy_waveform, int(round(1/float(values[parameters.index('lower_frequency')])))))
    for m in range(len(shifted_downsampled_streams)):
        # Calculate horizontal total energy for the station stream
        downsampled_sss = shifted_downsampled_streams[m].copy()
        shifted_stream_total_energy_waveform = calculate_total_energy(downsampled_sss, phase)
        shifted_stream_total_energy_waveform = np.asarray(smooth_data(
            shifted_stream_total_energy_waveform, int(round(1/float(values[parameters.index('lower_frequency')])))))
        # Initiate one loop to work through each possible start time in the waveform
        normalised_xcorr_values = [0] * len(shifted_stream_total_energy_waveform)
        window_length = ((int(values[parameters.index(phase + '_window')]) - 10) *
                         int(shifted_downsampled_streams[m][0].stats.sampling_rate))
        for n in range(nandices[0], len(shifted_stream_total_energy_waveform) - window_length):
            if nandices[1] and n > nandices[1]:  # Don't do cross-correlation past the data
                break
            #  Investigate each search window  at least 1 of the phase search windows long after the current n
            #  (minus 10 seconds on either due to trimming post-filtering).
            o = n + window_length
            if nandices[1] and o > nandices[1]:  # Don't do cross-correlation past the data
                break
            # Calculate mean, variance for data in the given window
            x_mean = np.nanmean(shifted_stream_total_energy_waveform[n:o])
            y_mean = np.nanmean(reference_total_energy_waveform[n:o])
            x_var = np.nanvar(shifted_stream_total_energy_waveform[n:o])
            y_var = np.nanvar(reference_total_energy_waveform[n:o])
            if x_var == 0 or y_var == 0:
                continue
            if np.isnan(x_mean) or np.isnan(y_mean) or np.isnan(x_var) or np.isnan(y_var):
                continue
            # Iterate through all values in the given window
            sum = 0
            for p in range(n, o):
                # Skip this window if there are any nan values
                if (np.isnan(shifted_stream_total_energy_waveform[p]) or
                        ma.is_masked(shifted_stream_total_energy_waveform[p]) or
                        np.isnan(reference_total_energy_waveform[p]) or
                        ma.is_masked(reference_total_energy_waveform[p])):
                    break
                sum += ((shifted_stream_total_energy_waveform[p] - x_mean) *
                        (reference_total_energy_waveform[p] - y_mean))
            normalised_xcorr_value = (1 / len(shifted_stream_total_energy_waveform[n:o]) *
                                      sum / math.sqrt(x_var * y_var))
            # Store normalised cross-correlation values in a nested list where the first index is the window
            # start index and the second index is the window end index.
            normalised_xcorr_values[n] = normalised_xcorr_value
        if np.nanmax(normalised_xcorr_values) == 0:
            print('Seismograms failed to find any suitable correlation window!')
            return np.nan
        max_normalised_xcorr_value_idx = normalised_xcorr_values.index(max(normalised_xcorr_values))
        print('For station ' + shifted_downsampled_streams[m][0].stats.station +
              ' maximum normalised cross-correlation value occurs between times ' +
              str(downsampled_sss[0].times(type='utcdatetime')[max_normalised_xcorr_value_idx]) +
              ' (index ' + str(max_normalised_xcorr_value_idx) + ' in downsampled data) - ' +
              str(downsampled_sss[0].times(type='utcdatetime')[max_normalised_xcorr_value_idx + window_length]) +
              ' in the aligned data (index ' + str(max_normalised_xcorr_value_idx + window_length) +
              ' in the downsampled data)')
        print('Maximum cross-correlation value is: ' + str(max(normalised_xcorr_values)))
        if nandices[0] > max_normalised_xcorr_value_idx:
            print('There are ' + str(nandices[0] - max_normalised_xcorr_value_idx) +
                  ' NaN values at the front of the cross-correlation window in the downsampled aligned data')
        if nandices[1] and nandices[1] < max_normalised_xcorr_value_idx + window_length:
            print('There are ' + str(max_normalised_xcorr_value_idx + window_length - nandices[1]) +
                  ' NaN values at the end of the cross-correlation window in the downsampled aligned data')
        xcorr_window.append([downsampled_sss[0].times(type='utcdatetime')[
                                 max_normalised_xcorr_value_idx],
                             downsampled_sss[0].times(type='utcdatetime')[
                                 max_normalised_xcorr_value_idx + window_length]])

    return xcorr_window


def generate_polarity_xcorr_values(args):

    """
    Generate the range of polarity values through cross-correlation of event data at a given station
    with that of the reference station.
    """

    event, station, parameters, values = args

    print('\nevent is ' + str(event))
    # Query event details from GeoNet FDSN
    event_details = get_event(event)[0]
    p_times, s_times = [], []
    p_stations, s_stations = [], []
    for pick in event_details.picks:
        if '.' + station + '.' in str(pick.resource_id):
            if pick.phase_hint == 'P':
                p_times.append(pick.time.datetime)
                p_stations.append(station)
            elif pick.phase_hint == 'S':
                s_times.append(pick.time.datetime)
                s_stations.append(station)

    # If event does not have picked phase arrivals, create theoretical ones
    hypocentre = [event_details.origins[0].latitude,
                  event_details.origins[0].longitude,
                  event_details.origins[0].depth]
    origin_time = event_details.origins[0].time
    try:
        station_location = FDSN_station_query(station)
    except:
        # If this fails, the station is not in FDSN. Ignore it.
        print('Failed to query station ' + station + ' from FDSN; will not include in analysis')
        return
    # Calculate epicentral distance in degrees
    delta = math.degrees(2 * (math.asin(((math.sin(1 / 2 * math.radians(abs(hypocentre[0] -
                                                                            station_location[0])))) ** 2 +
                                         math.cos(math.radians(hypocentre[0])) *
                                         math.cos(math.radians(station_location[0])) *
                                         (math.sin(1 / 2 * math.radians(abs(hypocentre[1] -
                                                                            station_location[1])))) ** 2) **
                                        (1 / 2))))
    # Calculate theoretical travel times
    arrivals = spherical_velocity_model.get_travel_times(source_depth_in_km=hypocentre[2] / 1000.0,
                                                         receiver_depth_in_km=max(
                                                             0.0, (station_location[2] / 1000.0)),
                                                         distance_in_degree=delta,
                                                         phase_list=['p', 'P', 's', 'S'])
    # Extract travel times
    p_tt, s_tt = None, None
    for arrival in arrivals:
        if arrival.name == 'p' or arrival.name == 'P' and p_tt is None:
            p_tt = arrival.time
        if arrival.name == 's' or arrival.name == 'S' and s_tt is None:
            s_tt = arrival.time
    # Save travel times, if none exist for the station
    if station not in p_stations:
        p_times.append(origin_time + datetime.timedelta(seconds=p_tt))
        p_stations.append(station)
    if station not in s_stations:
        s_times.append(origin_time + datetime.timedelta(seconds=s_tt))
        s_stations.append(station)

    # Get data around P phase
    phase = 'p'
    reference_station_stream, downsampled_rss, streams, downsampled_streams = \
        get_data(event, p_times, [values[parameters.index('reference_station')], station], phase, parameters, values)

    print('All P-phase data has been loaded into the script. '
          'Finding shift times to apply to each station to align them with the reference station...')

    # Align P phase data
    nandices, reference_station_stream, downsampled_rss, shifted_streams, shifted_downsampled_streams = \
        align_seismograms(p_stations, p_times, streams, downsampled_streams, reference_station_stream,
                          downsampled_rss, phase, parameters, values)

    print('All P-phase seismograms have now been aligned. '
          'Searching now for the optimal cross-correlation windows...')

    xcorr_window = find_xcorr_window(shifted_downsampled_streams, downsampled_rss, nandices, phase, parameters,
                                     values)
    if xcorr_window is np.nan:
        print('This event will not continue through the analysis.')
        return

    print('All P-phase cross-correlation windows are now defined. '
          'The orientation angles will now be calculated...')

    # Perform the orientation angle calculation in the cross-correlation window
    for m in range(len(shifted_streams)):
        print('Station is ' + str(shifted_streams[m][0].stats.station))

        # Make copies of data streams for shifting
        trimmed_and_shifted_streams = shifted_streams[m].copy()
        trimmed_reference_stream = reference_station_stream.copy()

        # Trim data to window
        window_start = xcorr_window[m][0]
        window_end = xcorr_window[m][1]
        trimmed_and_shifted_streams.trim(starttime=window_start,
                                         endtime=window_end)
        trimmed_reference_stream.trim(starttime=window_start,
                                      endtime=window_end)

        # Show xcorr window over trimmed and shifted data

        plt.subplot(411)
        plt.plot(reference_station_stream[0].times(type='utcdatetime'), reference_station_stream[0].data)
        plt.plot(trimmed_reference_stream[0].times(type='utcdatetime'), trimmed_reference_stream[0].data, color='red')
        plt.subplot(412)
        plt.plot(reference_station_stream[1].times(type='utcdatetime'), reference_station_stream[1].data)
        plt.plot(trimmed_reference_stream[1].times(type='utcdatetime'), trimmed_reference_stream[1].data, color='red')
        plt.subplot(413)
        plt.plot(shifted_streams[m][0].times(type='utcdatetime'), shifted_streams[m][0].data)
        plt.plot(trimmed_and_shifted_streams[0].times(type='utcdatetime'), trimmed_and_shifted_streams[0].data, color='red')
        plt.subplot(414)
        plt.plot(shifted_streams[m][1].times(type='utcdatetime'), shifted_streams[m][1].data)
        plt.plot(trimmed_and_shifted_streams[1].times(type='utcdatetime'), trimmed_and_shifted_streams[1].data, color='red')
        plt.savefig(station + '_' + event + '_' + phase + '_3-xcorr-data.png')
        plt.clf()
        time.sleep(30)  # Allow plotting in parallel threads to finish

        # Retain only z component data
        for n, tr in enumerate(trimmed_and_shifted_streams):
            if tr.stats.channel[-1] == 'Z':
                trimmed_and_shifted_streams = [tr]
        for n, tr in enumerate(trimmed_reference_stream):
            if tr.stats.channel[-1] == 'Z':
                trimmed_reference_stream = [tr]

        # Calculate angle
        polarity, xcorr_values = find_rotation_angle(trimmed_and_shifted_streams, trimmed_reference_stream, phase)
        print('Polarity is ' + str(polarity))
        np.save(shifted_streams[m][0].stats.station + '_' + event + '_' + phase + '_xcorr_values.npy', xcorr_values)


def generate_rotation_xcorr_angles(args):

    """
    Generate the range of rotation correlation values through cross-correlation of event data at a given station
    with that of the reference station.
    """

    event, station, polarity, parameters, values = args

    print('\nevent is ' + str(event))
    # Query event details from GeoNet FDSN
    event_details = get_event(event)[0]
    p_times, s_times = [], []
    p_stations, s_stations = [], []
    for pick in event_details.picks:
        if '.' + station + '.' in str(pick.resource_id):
            if pick.phase_hint == 'P':
                p_times.append(pick.time.datetime)
                p_stations.append(station)
            elif pick.phase_hint == 'S':
                s_times.append(pick.time.datetime)
                s_stations.append(station)

    # If event does not have picked phase arrivals, create theoretical ones
    hypocentre = [event_details.origins[0].latitude,
                  event_details.origins[0].longitude,
                  event_details.origins[0].depth]
    origin_time = event_details.origins[0].time
    try:
        station_location = FDSN_station_query(station)
    except:
        # If this fails, the station is not in FDSN. Ignore it.
        print('Failed to query station ' + station + ' from FDSN; will not include in analysis')
        return
    # Calculate epicentral distance in degrees
    delta = math.degrees(2 * (math.asin(((math.sin(1 / 2 * math.radians(abs(hypocentre[0] -
                                                                            station_location[0])))) ** 2 +
                                         math.cos(math.radians(hypocentre[0])) *
                                         math.cos(math.radians(station_location[0])) *
                                         (math.sin(1 / 2 * math.radians(abs(hypocentre[1] -
                                                                            station_location[1])))) ** 2) **
                                        (1 / 2))))
    # Calculate theoretical travel times
    arrivals = spherical_velocity_model.get_travel_times(source_depth_in_km=hypocentre[2] / 1000.0,
                                                         receiver_depth_in_km=max(
                                                             0.0, (station_location[2] / 1000.0)),
                                                         distance_in_degree=delta,
                                                         phase_list=['p', 'P', 's', 'S'])
    # Extract travel times
    p_tt, s_tt = None, None
    for arrival in arrivals:
        if arrival.name == 'p' or arrival.name == 'P' and p_tt is None:
            p_tt = arrival.time
        if arrival.name == 's' or arrival.name == 'S' and s_tt is None:
            s_tt = arrival.time
    # Save travel times, if none exist for the station
    if station not in p_stations:
        p_times.append(origin_time + datetime.timedelta(seconds=p_tt))
        p_stations.append(station)
    if station not in s_stations:
        s_times.append(origin_time + datetime.timedelta(seconds=s_tt))
        s_stations.append(station)

    # Get data around S phase
    phase = 's'
    reference_station_stream, downsampled_rss, streams, downsampled_streams = \
        get_data(event, s_times, [values[parameters.index('reference_station')], station], phase, parameters, values)

    # Scale data by polarity
    for tr in reference_station_stream:
        tr.data *= polarity
    for tr in downsampled_rss:
        tr.data *= polarity
    for stream in streams:
        for tr in stream:
            tr.data *= polarity
    for stream in downsampled_streams:
        for tr in stream:
            tr.data *= polarity

    print('All S-phase data has been loaded into the script. '
          'Finding shift times to apply to each station to align them with the reference station...')

    # Align S phase data
    nandices, reference_station_stream, downsampled_rss, shifted_streams, shifted_downsampled_streams = \
        align_seismograms(s_stations, s_times, streams, downsampled_streams, reference_station_stream,
                          downsampled_rss, phase, parameters, values)

    print('All S-phase seismograms have now been aligned. '
          'Searching now for the optimal cross-correlation windows...')

    xcorr_window = find_xcorr_window(shifted_downsampled_streams, downsampled_rss, nandices, phase, parameters,
                                     values)
    if xcorr_window is np.nan:
        print('This event will not continue through the analysis.')
        return

    print('All S-phase cross-correlation windows are now defined. '
          'The orientation angles will now be calculated...')

    # Perform the orientation angle calculation in the cross-correlation window
    for m in range(len(shifted_streams)):
        print('Station is ' + str(shifted_streams[m][0].stats.station))

        # Make copies of data streams for shifting
        trimmed_and_shifted_streams = shifted_streams[m].copy()
        trimmed_reference_stream = reference_station_stream.copy()

        # Trim data to window
        window_start = xcorr_window[m][0]
        window_end = xcorr_window[m][1]
        trimmed_and_shifted_streams.trim(starttime=window_start,
                                         endtime=window_end)
        trimmed_reference_stream.trim(starttime=window_start,
                                      endtime=window_end)

        # Show xcorr window over trimmed and shifted data

        plt.subplot(411)
        plt.plot(reference_station_stream[0].times(type='utcdatetime'), reference_station_stream[0].data)
        plt.plot(trimmed_reference_stream[0].times(type='utcdatetime'), trimmed_reference_stream[0].data, color='red')
        plt.subplot(412)
        plt.plot(reference_station_stream[1].times(type='utcdatetime'), reference_station_stream[1].data)
        plt.plot(trimmed_reference_stream[1].times(type='utcdatetime'), trimmed_reference_stream[1].data, color='red')
        plt.subplot(413)
        plt.plot(shifted_streams[m][0].times(type='utcdatetime'), shifted_streams[m][0].data)
        plt.plot(trimmed_and_shifted_streams[0].times(type='utcdatetime'), trimmed_and_shifted_streams[0].data, color='red')
        plt.subplot(414)
        plt.plot(shifted_streams[m][1].times(type='utcdatetime'), shifted_streams[m][1].data)
        plt.plot(trimmed_and_shifted_streams[1].times(type='utcdatetime'), trimmed_and_shifted_streams[1].data, color='red')
        plt.savefig(station + '_' + event + '_' + phase + '_3-xcorr-data.png')
        plt.clf()
        time.sleep(30)  # Allow plotting in parallel threads to finish

        # Remove all z component data
        for n, tr in enumerate(trimmed_and_shifted_streams):
            if tr.stats.channel[-1] == 'Z':
                trimmed_and_shifted_streams.pop(n)
        for n, tr in enumerate(trimmed_reference_stream):
            if tr.stats.channel[-1] == 'Z':
                trimmed_reference_stream.pop(n)

        # Calculate angle
        orientation_angle, xcorr_values = find_rotation_angle(trimmed_and_shifted_streams,
                                                              trimmed_reference_stream, phase)
        print('Rotation angle is ' + str(orientation_angle))
        np.save(shifted_streams[m][0].stats.station + '_' + event + '_' + phase + '_xcorr_values.npy', xcorr_values)


if __name__ == "__main__":

    # Parse data from config file
    parameters = ['eventid_file', 'station_file', 'station_channels', 'lower_frequency', 'upper_frequency', 'p_window',
                  's_window', 'reference_station', 'reference_channels', 'local_reference', 'reference_data_dir']
    values = [0] * len(parameters)
    config_file = os.path.dirname(os.path.realpath(__file__)) + '/config-find_sensor_orientation.txt'
    with open(config_file, 'r') as openfile:
        for row in openfile:
            if row.split(' ')[0] in parameters:
                idx = parameters.index(row.split(' ')[0])
                values[idx] = row.split(' ')[2]
    print('Input parameter values are: ')
    for n in range(len(parameters)):
        print(parameters[n] + ' = ' + values[n])

    # Parse files from paths given in config file
    stations_to_orient = parse_csv(values[parameters.index('station_file')], header=False)
    stations_to_orient = [item for sublist in stations_to_orient for item in sublist]
    events = parse_csv(values[parameters.index('eventid_file')], header=False)
    events = [item for sublist in events for item in sublist]
    events.sort()

    # Define list of stations to query data from FDSN for
    query_stations = []
    if values[parameters.index('reference_station')] not in stations_to_orient:
        query_stations.append(values[parameters.index('reference_station')])
        query_stations.extend(stations_to_orient)
    else:
        query_stations = stations_to_orient

    # For each station, derive a polarity from all events

    args = []
    for station in stations_to_orient:
        for event in events:
            # Send jobs to the pool of workers
            args.append([event, station, parameters, values])

    pool = multiprocessing.Pool(8)
    pool.map(generate_polarity_xcorr_values, args)
    pool.close()
    pool.join()

    print('\nCalculating stacked polarities for stations...')
    site_polarities = []
    for m, station in enumerate(stations_to_orient):
        site_polarity_xcorr_values = [0] * 2
        nevents = 0
        for n, event in enumerate(events):
            try:
                xcorr_values = np.load(station + '_' + event + '_p_xcorr_values.npy')
                nevents += 1
            except FileNotFoundError:
                print('Failed to find file ' + station + '_' + event + '_p_xcorr_values.npy')
                continue
            for o in range(len(xcorr_values)):
                site_polarity_xcorr_values[o] += xcorr_values[o]

        print('Number of events used for polarity determination is ' + str(nevents))
        plt.scatter([-1, 1], site_polarity_xcorr_values)
        plt.savefig(stations_to_orient[m] + '_p_xcorr-values.png')

        # Give value from all data
        max_xcorr_value = max(site_polarity_xcorr_values)
        max_xcorr_value_idx = site_polarity_xcorr_values.index(max_xcorr_value)
        print('Stacked peak polarity is: ' + str([1, -1][max_xcorr_value_idx]))
        print('\n')
        polarity = [1, -1][max_xcorr_value_idx]
        site_polarities.append(polarity)

    # For each station, derive an orientation angle from all events

    args = []
    for m, station in enumerate(stations_to_orient):
        for event in events:
            # Send jobs to the pool of workers
            args.append([event, station, site_polarities[m], parameters, values])

    pool = multiprocessing.Pool(8)
    pool.map(generate_rotation_xcorr_angles, args)
    pool.close()
    pool.join()

    # For each station, derive a rotation angle from all events and save to file

    print('\nCalculating stacked orientation angles for stations...')
    site_rotation_angles = []
    for m, station in enumerate(stations_to_orient):
        site_rotation_xcorr_values = [0] * 360
        nevents = 0
        for n, event in enumerate(events):
            try:
                xcorr_values = np.load(station + '_' + event + '_s_xcorr_values.npy')
                nevents += 1
            except FileNotFoundError:
                print('Failed to find file ' + station + '_' + event + '_s_xcorr_values.npy')
                continue
            for o in range(len(xcorr_values)):
                site_rotation_xcorr_values[o] += xcorr_values[o]
            plt.scatter(list(range(len(xcorr_values))), xcorr_values)
            plt.savefig(station + '_' + event + '_4-xcorr-values.png')
            plt.clf()

        print('Number of events used for rotation determination is ' + str(nevents))
        plt.scatter(list(range(0, 360)), site_rotation_xcorr_values)
        plt.savefig(stations_to_orient[m] + '_s_xcorr-values.png')

        # Give value from all data
        max_xcorr_value = max(site_rotation_xcorr_values)
        max_xcorr_value_idx = site_rotation_xcorr_values.index(max_xcorr_value)
        print('Stacked peak rotation value is: ' + str(max_xcorr_value_idx))
        site_rotation_angles.append(max_xcorr_value_idx)

    # Write results to file
    with open('find_sensor_orientation_output.csv', 'w') as outfile:
        outfile.write('Station,Polarity,Orientation Angle\n')
    with open('find_sensor_orientation_output.csv', 'a') as outfile:
        for m in range(len(stations_to_orient)):
            outfile.write(str(stations_to_orient[m]) + ',' + str(site_polarities[m]) + ',' +
                          str(site_rotation_angles[m]) + '\n')
