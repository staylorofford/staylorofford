"""
Script to find a sensor's orientation using teleseismic events.
"""

import datetime
import io
import math
import numpy as np
import obspy
from obspy.clients.fdsn import Client
from obspy.io.quakeml.core import Unpickler
import os
import pycurl

quakeml_reader = Unpickler()


def parse_csv(csv_path):

    """
    Parse the contents of CSV file at csv_path into nested lists
    :param csv_path: full path to CSV file
    :return: list of CSV header columns, nested lists of CSV row data
    """

    data = []
    with open(csv_path, 'r') as openfile:
        rc = -1
        for row in openfile:
            if rc == -1:
                header = row[:-1].split(',')
                rc = 0
            else:
                data.append(row[:-1].split(','))
    return header, data


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

    query = "https://service.geonet.org.nz/fdsnws/dataselect/1/query?network=NZ&" + eventID
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
                                  attach_response=True)
    return stream


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
        hte_waveform.append(stream[0].data[m] ** 2 + stream[1].data[m] ** 2)
    return hte_waveform


def find_lag_time(stream, reference_stream):

    """
    Find the lag time between a given stream and a reference stream using cross-correlation of the horizontal total
    energy trace in each stream.
    :param stream: obspy stream object of the seismogram to calculate lag time for
    :param reference_stream: object stream object of the seismogram to use as reference for lag time calculation
    :return: lag time between the two sensors in seconds relative to the reference stream
    """

    # Calculate horizontal total energy waveforms
    hte_waveform = calculate_horizontal_total_energy(stream)
    hte_reference_waveform = calculate_horizontal_total_energy(reference_stream)

    # Find the lag time from the maximum cross-correlation value between the two waveforms
    xcorr_values = []
    for m in range(len(hte_waveform)):
        # Shift the non-reference waveform temporally
        shifted_hte_waveform = [float('nan')] * m + hte_waveform[:m]
        # Cross-correlate the waveforms
        xcorr_value = sum(np.correlate(shifted_hte_waveform, hte_reference_waveform))
        xcorr_values.append(xcorr_value)
    max_xcorr_value = max(xcorr_values)
    lag_time = (1 / stream[0].stats.sampling_rate) * xcorr_values.index(max_xcorr_value)

    # Perform the same operation backwards to be sure there is not a better match
    xcorr_values = []
    for m in range(len(hte_waveform)):
        # Shift the non-reference waveform temporally
        shifted_hte_waveform = hte_waveform[m:] + [float('nan')] * m
        # Cross-correlate the waveforms
        xcorr_value = sum(np.correlate(shifted_hte_waveform, hte_reference_waveform))
        xcorr_values.append(xcorr_value)
    if max(xcorr_values) > max_xcorr_value:
        lag_time = (-1 / stream[0].stats.sampling_rate) * xcorr_values.index(max(xcorr_values))

    return lag_time


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
        numerator += ((reference_seismogram[m][0] * shifted_seismogram[m][0]) +
                      (reference_seismogram[m][1] * shifted_seismogram[m][1]))
        denominator += ((reference_seismogram[m][1] * shifted_seismogram[m][0]) -
                        (reference_seismogram[m][0] * shifted_seismogram[m][1]))
    ratio = numerator / denominator

    # Now, calculate the relative angle from the inverse-tangent of the ratio
    relative_angle = math.atan(ratio)
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
    stations_to_orient = parse_csv(values[parameters.index('station_file')])
    events = parse_csv(values[parameters.index('eventid_file')])

    # Define list of stations to query data from FDSN for
    query_stations = []
    if values[parameters.index('reference_stations')] not in stations_to_orient:
        query_stations.append(values[parameters.index('reference_stations')])
        query_stations.extend(stations_to_orient)
    else:
        query_stations = stations_to_orient

    # For each event, produce the orientation angle for each station
    all_orientation_angles = []
    for event in events:
        # Query event details from GeoNet FDSN
        event_details = get_event(event)[0]
        start_time = event_details.origins[0].time.datetime
        # Assume event duration does not exceed 5 minutes at any station
        end_time = start_time + datetime.datetime.timedelta(seconds=300)

        # Query waveform data from GeoNet FDSN for each station between the start and end times defined for the event
        stations = []
        streams = []
        for station in query_stations:
            stations.append(station)
            station_stream = query_fdsn(station,
                                        '??',
                                        values[parameters.index('channels')],
                                        start_time,
                                        end_time)
            # Remove instrument response so seismograms are sensor-independent
            station_stream.remove_response(output='VEL')

            # Resample stream if sampling rate is not equal to that desired
            if station_stream[0].stats.sampling_rate != values[parameters.index('sampling_rate')]:
                for trace in station_stream:
                    trace.resample(sampling_rate=values[parameters.index('sampling_rate')])

            # Low pass the waveforms of all events to increase waveform similarly at all sensors:
            # Filter corner frequency satisfies the condition that it is much smaller than the lowest seismic velocity in the
            # propagation medium of all events scaled by the linear distance between each pair of sensors
            station_stream.filter(type='lowpass',
                                  freq=values[parameters.index('corner_frequency')])
            streams.append(station_stream)
        reference_station_idx = stations.index(values[parameters.index('reference_station')])
        reference_station_stream = streams.pop(reference_station_idx)

        # Align all seismograms for each event to facilitate cross-correlation:
        # Use the lag time that produces the maximum cross-correlation value between each sensor and the reference sensor's
        # total energy traces (sum of squares of horizontal seismograms).
        shifted_streams = streams.copy()
        for m in range(len(streams)):
            lag_time = find_lag_time(streams[m], reference_station_stream)
            shift_idx = abs(lag_time * streams[m][0].stats.sampling_rate)
            for n in range(len(streams[m])):
                if lag_time >= 0:
                    shifted_streams[m][n].data = [float('nan')] * shift_idx + streams[m][n].data[:shift_idx]
                else:
                    shifted_streams[m][n].data = streams[m][n].data[shift_idx:] + [float('nan') * shift_idx]

        # Find the time window of each event for which the correlation of the horizontal waveforms at the sensors is highest:
        # This window is that for which the normalised cross-correlation of the total energy traces of each sensor pair is
        # highest. Total energy traces are independent of sensor orientation, so the value of this normalised cross-correlation
        # is indicative of waveform similarity alone.
        xcorr_window_idx = []
        reference_horizontal_total_energy_waveform = calculate_horizontal_total_energy(reference_station_stream)
        for m in range(len(shifted_streams)):
            shifted_stream_horizontal_total_energy_waveform = calculate_horizontal_total_energy(shifted_streams[m])
            # Initiate one loop to work through each possible start time in the waveform
            normalised_xcorr_values = [[[] for x in range(len(reference_horizontal_total_energy_waveform))]
                                       for y in range(len(shifted_stream_horizontal_total_energy_waveform))]
            for n in range(len(shifted_stream_horizontal_total_energy_waveform)):
                # Initiate a second loop to work through each possible end time in the waveform,
                # so that all possible windows are tested.
                for o in range(len(shifted_stream_horizontal_total_energy_waveform)):
                    # Calculate mean, variance for data in the given window
                    x_mean = np.nanmean(shifted_stream_horizontal_total_energy_waveform[n:o])
                    y_mean = np.nanmean(reference_horizontal_total_energy_waveform[n:o])
                    x_var = np.nanvar(shifted_stream_horizontal_total_energy_waveform[n:o])
                    y_var = np.nanvar(reference_horizontal_total_energy_waveform[n:o])
                    # Iterate through all values in the given window
                    sum = 0
                    for p in range(o - n):
                        sum += ((shifted_stream_horizontal_total_energy_waveform[p] - x_mean) *
                                (reference_horizontal_total_energy_waveform[p] - y_mean))
                    normalised_xcorr_value = (1 / len(shifted_stream_horizontal_total_energy_waveform[n:o]) *
                                              sum / math.sqrt(x_var * y_var))
                    # Store normalised cross-correlation values in a nested list where the first index is the window
                    # start index and the second index is the window end index.
                    normalised_xcorr_values[n][o] = normalised_xcorr_value
            normalised_xcorr_values = np.array(normalised_xcorr_values)  # Convert to numpy array to use numpy functions
            max_normalised_xcorr_value = np.nanmax(normalised_xcorr_values)
            max_normalised_xcorr_value_idx = np.unravel_index(max_normalised_xcorr_value, normalised_xcorr_values.shape)
            xcorr_window_idx.append(max_normalised_xcorr_value_idx)

        # Perform the orientation angle calculation in the cross-correlation window
        orientation_angles = []
        for m in range(len(shifted_streams)):
            # Make copies of data streams for shifting
            trimmed_and_shifted_streams = shifted_streams[m].copy()
            trimmed_reference_stream = reference_station_stream.copy()
            # Convert window indices to start and end times for stream trimming prior to angle calculation
            window_start = 1 / trimmed_and_shifted_streams[0].stats.sampling_rate * xcorr_window_idx[m][0]
            window_end = 1 / trimmed_and_shifted_streams[0].stats.sampling_rate * xcorr_window_idx[m][1]
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
    site_orientation_angles = [[] for n in range(len(query_stations))]
    for m in range(len(all_orientation_angles)):
        for n in range(len(all_orientation_angles[m])):
            site = query_stations[n]
            site_orientation_angles[n].append(all_orientation_angles[m][n])
    print('\nNumber of events used for orientation is ' + str(len(events)))
    for n in range(len(site_orientation_angles)):
        print('Site is:' + str(query_stations[n]))
        print('Mean orientation angle is: ' + str(np.mean(site_orientation_angles)))
        print('Median orientation angle is: ' + str(np.median(site_orientation_angles)))
        print('Orientation angle mode is: ' + str(np.mode(site_orientation_angles)))
        print('Orientation angle stdev is: ' + str(np.std(site_orientation_angles)))
        print('\n')
