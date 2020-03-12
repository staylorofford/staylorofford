"""
Script to find a sensor's orientation using teleseismic events.
"""

import datetime
import io
import math
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import obspy
from obspy.clients.fdsn import Client
from obspy.io.quakeml.core import Unpickler
import os
import pycurl
import scipy
from scipy.fftpack import hilbert

quakeml_reader = Unpickler()


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
                                  endtime=endtime),
                                  # attach_response=True)
    return stream


def calculate_envelope(waveform):

    """
    Find the amplitude envelope for a waveform. Waveform can be in a numpy masked array.
    :param waveform: waveform data in a list, numpy array, or numpy masked array
    :return: amplitude envelope of data, including masks if input
    """

    # Find period of data with no mask values and calculate the envelope for this data
    maskdices = [None, None]
    for n in range(len(waveform) - 1):
        if ma.is_masked(waveform[n]) and not ma.is_masked(waveform[n + 1]):
            maskdices[0] = n + 1
    for n in range(1, len(waveform)):
        if ma.is_masked(waveform[n]) and not ma.is_masked(waveform[n - 1]):
            maskdices[1] = n
    if maskdices[0] is not None and maskdices[1] is not None:
        hilb = hilbert(waveform[maskdices[0]:maskdices[1]])
        waveform_envelope = (waveform[maskdices[0]:maskdices[1]] ** 2 + hilb ** 2) ** 0.5
        waveform[maskdices[0]:maskdices[1]] = waveform_envelope
    elif maskdices[0] is not None:
        hilb = hilbert(waveform[maskdices[0]:])
        waveform_envelope = (waveform[maskdices[0]:] ** 2 + hilb ** 2) ** 0.5
        waveform[maskdices[0]:] = waveform_envelope
    else:
        hilb = hilbert(waveform)
        waveform_envelope = (waveform ** 2 + hilb ** 2) ** 0.5
        waveform = waveform_envelope
    return waveform


def find_lag_time(stream, reference_stream):

    """
    Find the lag time between a given stream and a reference stream using cross-correlation of the vertical traces.
    :param stream: obspy stream object of the seismogram to calculate lag time for
    :param reference_stream: object stream object of the seismogram to use as reference for lag time calculation
    :return: lag time between the two sensors in seconds relative to the reference stream
    """

    print(stream, reference_stream)

    # Get vertical component data
    for tr in stream:
        if tr.stats.channel[-1] == 'Z':
            stream_z = tr.data
    for tr in reference_stream:
        if tr.stats.channel[-1] == 'Z':
            ref_z = tr.data

    # Create normalised amplitude envelopes of the data
    stream_envelope = calculate_envelope(stream_z)
    max_se = ma.max(stream_envelope)
    stream_envelope /= max_se
    ref_envelope = calculate_envelope(ref_z)
    max_re = ma.max(ref_envelope)
    ref_envelope /= max_re

    plt.plot(ref_envelope, color='b')
    plt.plot(stream_envelope, color='r', alpha=0.4)
    plt.show()

    # Find the lag time from the maximum cross-correlation value between the two waveforms
    xcorr_values = []
    ref_envelope = ref_envelope.tolist() + len(stream_envelope) * [float('nan')]
    for m in range(2 * len(stream_envelope)):

        # Shift the stream
        if m <= len(stream_envelope):
            shifted_stream_envelope = (len(stream_envelope) - m) * [float('nan')] + \
                                      stream_envelope.filled(float('nan')).tolist() + \
                                      m * [float('nan')]
        else:
            shifted_stream_envelope = max(0, len(stream_envelope) - m) * [float('nan')] + \
                                      stream_envelope[:len(stream_envelope) - m].filled(float('nan')).tolist() + \
                                      m * [float('nan')]

        # Perform cross-correlation
        xcorr_value = ma.corrcoef(shifted_stream_envelope, ref_envelope)[0][1]
        xcorr_values.append(xcorr_value)

    # Find lag time from highest cross-correlation value
    max_xcorr_value = max(xcorr_values)
    print(xcorr_values)
    print(max_xcorr_value)
    lag_time = (1 / stream[0].stats.sampling_rate) * (len(stream_envelope) - xcorr_values.index(max_xcorr_value))

    print(lag_time)
    plt.plot(ref_envelope, color='b')
    plt.plot(stream_envelope, color='r', alpha=1)
    plt.plot(xcorr_values, color='purple', alpha=1)
    plt.show()

    return lag_time


def calculate_horizontal_total_energy(stream):

    """
    Produce the horizontal total energy trace for a given 3-component seismogram.
    :param stream: obspy stream object containing 3-component seismogram to calculate horizontal total energy for
    :return: list containing horizontal total energy waveform
    """

    # Remove vertical component data from stream, if it exists
    for m in range(len(stream)):
        if stream[m].stats.channel[-1] == 'Z':
            stream.pop(m)
            break

    # Calculate horizontal total energy waveform
    hte_waveform = []
    for m in range(len(stream[0].data)):
        try:
            hte_waveform.append(stream[0].data[m] ** 2 + stream[1].data[m] ** 2)
        except TypeError:
            hte_waveform.append(np.nan)
    return hte_waveform


def calculate_angle(shifted_seismogram, reference_seismogram):

    """
    Calculate the angle between the horizontal components of each sensor and those of the reference sensor.
    Uses the method of Zheng and McMehan (2006).
    :param shifted_seismogram: the seismogram from the station with unknown orientation, filtered, time-shifted, and 
           without instrument response.
    :param reference_seismogram: the seismogram from the reference station, filtered and without instrument response. 
    :return: angle between the horizontal components of the sensors at the two stations.
    """

    # First, calculate the ratio of numerator to denominator terms in equation A-5 of Zheng and McMechan (2006).
    numerator = 0
    denominator = 0
    for m in range(len(reference_seismogram)):  # Reference and shifted seismograms should have the same length
        try:
            numerator += ((reference_seismogram[m][0] * shifted_seismogram[m][0]) +
                          (reference_seismogram[m][1] * shifted_seismogram[m][1]))
            denominator += ((reference_seismogram[m][1] * shifted_seismogram[m][0]) -
                            (reference_seismogram[m][0] * shifted_seismogram[m][1]))
        except TypeError:  # Catch when nan values exist and skip operations for such data
            continue
    ratio = numerator / denominator

    # Now, calculate the relative angle from the inverse-tangent of the ratio
    relative_angle = math.degrees(math.atan(ratio))
    return relative_angle


if __name__ == "__main__":

    # Parse data from config file
    parameters = ['eventid_file', 'station_file', 'channels', 'corner_frequency', 'reference_station', 'sampling_rate']
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
    stations_to_orient = parse_csv(values[parameters.index('station_file')], header=False)[0]
    events = parse_csv(values[parameters.index('eventid_file')], header=False)[0]

    # Define list of stations to query data from FDSN for
    query_stations = []
    if values[parameters.index('reference_station')] not in stations_to_orient:
        query_stations.append(values[parameters.index('reference_station')])
        query_stations.extend(stations_to_orient)
    else:
        query_stations = stations_to_orient

    # For each event, produce the orientation angle for each station
    all_orientation_angles = []
    for event in events:
        # Query event details from GeoNet FDSN
        event_details = get_event(event)[0]
        start_times = []
        for pick in event_details.picks:
            for station in stations_to_orient:
                if '.' + station + '.' in str(pick.resource_id):
                    start_times.append(pick.time.datetime - datetime.timedelta(seconds=0))
        # Assume event duration does not exceed 2 minutes at any station
        start_time = min(start_times)
        end_time = start_time + datetime.timedelta(seconds=180)

        # Query waveform data from GeoNet FDSN for each station between the start and end times defined for the event
        stations = []
        streams = []
        min_times, max_times = [], []
        for station in query_stations:
            stations.append(station)
            station_stream = query_fdsn(station,
                                        '??',
                                        values[parameters.index('channels')],
                                        start_time.isoformat(),
                                        end_time.isoformat())[0]

            # Remove instrument response so seismograms are sensor-independent
            # station_stream.remove_response(output='VEL')

            # Filter the waveforms of all events to increase waveform similarly at all sensors:
            # Filter corner frequency satisfies the condition that it is much smaller than the lowest seismic velocity
            # in the propagation medium of all events scaled by the linear distance between each pair of sensors.
            station_stream.filter(type='lowpass',
                                  freq=float(values[parameters.index('corner_frequency')]))

            # Resample stream if sampling rate is not equal to that desired
            # if station_stream[0].stats.sampling_rate != values[parameters.index('sampling_rate')]:

            # Resample stream to twice the corner frequency
            for trace in station_stream:
                trace.resample(sampling_rate=2 * int(values[parameters.index('corner_frequency')]))
                min_times.append(trace.stats.starttime)
                max_times.append(trace.stats.endtime)

            streams.append(station_stream)

        # Pad out traces to ensure an equal number of data in each
        for stream in streams:
            stream.trim(starttime=min(min_times),
                        endtime=max(max_times),
                        pad=True,
                        fill_value=None)

        reference_station_idx = stations.index(values[parameters.index('reference_station')])
        stations.pop(reference_station_idx)
        reference_station_stream = streams.pop(reference_station_idx)

        # Align all seismograms for each event to facilitate cross-correlation:
        # Use the lag time that produces the maximum cross-correlation value between each sensor and the reference sensor's
        # total energy traces (sum of squares of horizontal seismograms).
        nandices = [0, None]
        shifted_streams = streams.copy()
        for m in range(len(streams)):
            lag_time = find_lag_time(streams[m], reference_station_stream)
            shift_idx = int(abs(lag_time * streams[m][0].stats.sampling_rate))
            for n in range(len(streams[m])):
                if lag_time >= 0:
                    nandices[0] = shift_idx
                    shifted_streams[m][n].data = np.asarray([float('nan')] * shift_idx +
                                                            streams[m][n].data[:shift_idx].tolist())
                else:
                    nandices[1] = len(streams[m][n].data) - shift_idx
                    shifted_streams[m][n].data = np.asarray(streams[m][n].data[shift_idx:].tolist() +
                                                            [float('nan')] * shift_idx)

        print(lag_time)
        plt.plot(reference_station_stream[0].data, color='b')
        plt.plot(shifted_streams[0][0].data, alpha=0.4, color='r')
        plt.show()

        # Find the time window of each event for which the correlation of the horizontal waveforms at the sensors is highest:
        # This window is that for which the normalised cross-correlation of the total energy traces of each sensor pair is
        # highest. Total energy traces are independent of sensor orientation, so the value of this normalised cross-correlation
        # is indicative of waveform similarity alone.
        xcorr_window_idx = []
        reference_horizontal_total_energy_waveform = calculate_horizontal_total_energy(reference_station_stream)
        for m in range(len(shifted_streams)):
            shifted_stream_horizontal_total_energy_waveform = calculate_horizontal_total_energy(shifted_streams[m])
            # Initiate one loop to work through each possible start time in the waveform
            normalised_xcorr_values = [([0] * len(reference_horizontal_total_energy_waveform))
                                       for y in range(len(shifted_stream_horizontal_total_energy_waveform))]
            for n in range(nandices[0], len(shifted_stream_horizontal_total_energy_waveform)):
                if nandices[1] and n > nandices[1]:  # Don't do cross-correlation past the data
                    break
                # Initiate a second loop to work through each possible end time in the waveform,
                # so that all possible windows are tested, BUT require that windows are at least 1 second in length.
                for o in range(n + int(shifted_streams[m][0].stats.sampling_rate) + 1,
                               len(reference_horizontal_total_energy_waveform)):
                    if nandices[1] and o > nandices[1]:  # Don't do cross-correlation past the data
                        break
                    # Calculate mean, variance for data in the given window
                    x_mean = np.nanmean(shifted_stream_horizontal_total_energy_waveform[n:o])
                    y_mean = np.nanmean(reference_horizontal_total_energy_waveform[n:o])
                    x_var = np.nanvar(shifted_stream_horizontal_total_energy_waveform[n:o])
                    y_var = np.nanvar(reference_horizontal_total_energy_waveform[n:o])
                    if np.isnan(x_mean) or np.isnan(y_mean) or np.isnan(x_var) or np.isnan(y_var):
                        continue
                    # Iterate through all values in the given window
                    sum = 0
                    for p in range(o - n):
                        if (np.isnan(shifted_stream_horizontal_total_energy_waveform[n + p]) or
                                ma.is_masked(shifted_stream_horizontal_total_energy_waveform[n + p]) or
                                np.isnan(reference_horizontal_total_energy_waveform[n + p]) or
                                ma.is_masked(reference_horizontal_total_energy_waveform[n + p])):
                            continue
                        sum += ((shifted_stream_horizontal_total_energy_waveform[n + p] - x_mean) *
                                (reference_horizontal_total_energy_waveform[n + p] - y_mean))
                    normalised_xcorr_value = (1 / len(shifted_stream_horizontal_total_energy_waveform[n:o]) *
                                              sum / math.sqrt(x_var * y_var))
                    # Store normalised cross-correlation values in a nested list where the first index is the window
                    # start index and the second index is the window end index.
                    normalised_xcorr_values[n][o] = normalised_xcorr_value
            if np.nanmax(normalised_xcorr_values) == 0:
                print('Seismograms failed to find any suitable correlation window!')
                exit()
            normalised_xcorr_values = np.asarray(normalised_xcorr_values)
            max_normalised_xcorr_value_idx = np.unravel_index(np.nanargmax(normalised_xcorr_values),
                                                              normalised_xcorr_values.shape)
            print('For station ' + stations[m] + ' maximum normalised cross-correlation value occurs between samples ' +
                  str(max_normalised_xcorr_value_idx[0]) + '-' + str(max_normalised_xcorr_value_idx[1]) +
                  ' in the aligned data')
            if nandices[0] > max_normalised_xcorr_value_idx[0]:
                print('There are ' + str(nandices[0] - max_normalised_xcorr_value_idx[0]) +
                      ' NaN values at the front of the cross-correlation window in the aligned data')
            if nandices[1] and nandices[1] < max_normalised_xcorr_value_idx[1]:
                print('There are ' + str(max_normalised_xcorr_value_idx[1] - nandices[1]) +
                      ' NaN values at the end of the cross-correlation window in the aligned data')
            xcorr_window_idx.append(max_normalised_xcorr_value_idx)

        # Perform the orientation angle calculation in the cross-correlation window
        orientation_angles = []
        for m in range(len(shifted_streams)):
            # Make copies of data streams for shifting
            trimmed_and_shifted_streams = shifted_streams[m].copy()
            trimmed_reference_stream = reference_station_stream.copy()
            # Convert window indices to start and end times for stream trimming prior to angle calculation
            window_start = trimmed_and_shifted_streams[0].stats.starttime + \
                           (1 / trimmed_and_shifted_streams[0].stats.sampling_rate * xcorr_window_idx[m][0])
            window_end = trimmed_and_shifted_streams[0].stats.starttime + \
                         (1 / trimmed_and_shifted_streams[0].stats.sampling_rate * xcorr_window_idx[m][1])
            # Trim data to window
            trimmed_and_shifted_streams.trim(starttime=window_start,
                                             endtime=window_end)
            trimmed_reference_stream.trim(starttime=window_start,
                                          endtime=window_end)
            # Calculate angle
            orientation_angle = calculate_angle(trimmed_and_shifted_streams, trimmed_reference_stream)
            orientation_angles.append(orientation_angle)

        # Save orientation angles for the event
        all_orientation_angles.append(orientation_angles)

    # Calculate average orientation angle and error from angle distributions
    site_orientation_angles = [[] for n in range(len(stations))]
    for m in range(len(all_orientation_angles)):
        for n in range(len(all_orientation_angles[m])):
            site = stations[n]
            site_orientation_angles[n].append(all_orientation_angles[m][n])
    print('\nNumber of events used for orientation is ' + str(len(events)))
    for n in range(len(site_orientation_angles)):
        print('Site is: ' + str(stations[n]))
        print('Mean orientation angle is: ' + str(np.mean(site_orientation_angles[n])))
        print('Median orientation angle is: ' + str(np.median(site_orientation_angles[n])))
        print('Orientation angle mode is: ' + str(scipy.stats.mode(site_orientation_angles[n])))
        print('Orientation angle stdev is: ' + str(np.std(site_orientation_angles[n])))
        print('\n')
