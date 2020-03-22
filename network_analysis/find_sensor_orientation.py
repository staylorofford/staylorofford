"""
Script to find a sensor's orientation using teleseismic events. Requires events be queryable from FDSN, and that the
stations to orient contain picks and waveform data for the events available in FDSN.
"""

import datetime
import glob
import fnmatch
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
    The input can contain masked and nan elements, but nan elements take priority.
    :param waveform: waveform data in a numpy array, or numpy masked array
    :return: amplitude envelope of data, including masks if input
    """

    # Find period of data with no mask values and/or nan values
    maskdices = [None, None]
    nandices = [None, None]
    for n in range(len(waveform) - 1):
        if ma.is_masked(waveform[n]) and not ma.is_masked(waveform[n + 1]):
            maskdices[0] = n + 1
        if np.isnan(waveform[n]) and not np.isnan(waveform[n + 1]):
            nandices[0] = n + 1
    for n in range(1, len(waveform)):
        if ma.is_masked(waveform[n]) and not ma.is_masked(waveform[n - 1]):
            maskdices[1] = n
        if np.isnan(waveform[n]) and not np.isnan(waveform[n - 1]):
            nandices[1] = n

    # Calculate the envelope for the data that is not nan or masked
    if nandices[0] is None and nandices[1] is None:
        if maskdices[0] is not None and maskdices[1] is not None:
            hilb = hilbert(waveform[maskdices[0]:maskdices[1]])
            waveform_envelope = (waveform[maskdices[0]:maskdices[1]] ** 2 + hilb ** 2) ** 0.5
            waveform[maskdices[0]:maskdices[1]] = waveform_envelope
        elif maskdices[0] is not None:
            hilb = hilbert(waveform[maskdices[0]:])
            waveform_envelope = (waveform[maskdices[0]:] ** 2 + hilb ** 2) ** 0.5
            waveform[maskdices[0]:] = waveform_envelope
        elif maskdices[1] is not None:
            hilb = hilbert(waveform[:maskdices[1]])
            waveform_envelope = (waveform[:maskdices[1]] ** 2 + hilb ** 2) ** 0.5
            waveform[:maskdices[1]] = waveform_envelope
        else:
            hilb = hilbert(waveform)
            waveform_envelope = (waveform ** 2 + hilb ** 2) ** 0.5
            waveform = waveform_envelope
    else:
        if nandices[0] is not None and nandices[1] is not None:
            hilb = hilbert(waveform[nandices[0]:nandices[1]])
            waveform_envelope = (waveform[nandices[0]:nandices[1]] ** 2 + hilb ** 2) ** 0.5
            waveform[nandices[0]:nandices[1]] = waveform_envelope
        elif nandices[0] is not None:
            hilb = hilbert(waveform[nandices[0]:])
            waveform_envelope = (waveform[nandices[0]:] ** 2 + hilb ** 2) ** 0.5
            waveform[nandices[0]:] = waveform_envelope
        elif nandices[1] is not None:
            hilb = hilbert(waveform[:nandices[1]])
            waveform_envelope = (waveform[:nandices[1]] ** 2 + hilb ** 2) ** 0.5
            waveform[:nandices[1]] = waveform_envelope
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

    # Ensure all data are masked arrays
    if not ma.is_masked(stream_envelope):
        stream_envelope = ma.masked_array(stream_envelope)
    if not ma.is_masked(ref_envelope):
        ref_envelope = ma.masked_array(ref_envelope)

    # Find the lag time from the maximum cross-correlation value between the two waveforms
    xcorr_values = []
    ref_envelope = ref_envelope.filled(0).tolist() + len(stream_envelope) * [0]
    ref_envelope = np.asarray(smooth_data(ref_envelope, int(values[parameters.index('corner_frequency')])))
    for m in range(2 * len(stream_envelope)):

        # Shift the stream
        if m <= len(stream_envelope):
            shifted_stream_envelope = (len(stream_envelope) - m) * [0] + stream_envelope.filled(0).tolist() + m * [0]
        else:
            shifted_stream_envelope = max(0, len(stream_envelope) - m) * [0] + stream_envelope[:len(stream_envelope) - m].filled(0).tolist() + m * [0]
        shifted_stream_envelope = np.asarray(smooth_data(shifted_stream_envelope,
                                                         int(values[parameters.index('corner_frequency')])))

        # Perform cross-correlation (using only the first 30 seconds of the data to focus on first arriving data)
        try:
            xcorr_value = np.corrcoef(shifted_stream_envelope, ref_envelope)[0][1]
        except ValueError:
            print('Cross-correlation failed! Perhaps one stream is a data point different to the other? '
                  'This can occur for certain corner frequency and data length combinations... It is a bug.')
        xcorr_values.append(xcorr_value)

    # Find lag time from highest cross-correlation value
    max_xcorr_value = max(xcorr_values)
    lag_time = (1 / stream[0].stats.sampling_rate) * (len(stream_envelope) - xcorr_values.index(max_xcorr_value))

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
            hte_waveform.append(float('nan'))
    return hte_waveform


def calculate_horizontal_total_energy_from_list(data):

    """
    Produce the horizontal total energy trace for a list containing lists of the 2 horizontal components
    :param data: list containing 2 horizontal components of a seismogram to calculate horizontal total energy for
    :return: list containing horizontal total energy waveform
    """

    # Calculate horizontal total energy waveform
    hte_waveform = []
    for m in range(len(data[0])):
        try:
            hte_waveform.append(data[0][m] ** 2 + data[1][m] ** 2)
        except TypeError:
            hte_waveform.append(float('nan'))
    return hte_waveform


def smooth_data(data, N):

    """
    Apply a running mean to data using the N samples around a given point to produce mean values.
    :param data: list of data to apply running mean to
    :param N: number of samples either side of each data point to use in the mean calculation
    :return: smoothed data in list
    """

    smoothed_data = [0] * len(data)
    for m in range(N, len(data) - N):
        smoothed_data[m] = np.nanmean(data[m - N: m + N])
    return smoothed_data


def calculate_angle(shifted_seismogram, reference_seismogram):

    """
    Calculate the angle between the horizontal components of each sensor and those of the reference sensor.
    Uses the method of Zheng and McMehan (2006).
    :param shifted_seismogram: the seismogram from the station with unknown orientation, filtered, time-shifted, and 
           without instrument response.
    :param reference_seismogram: the seismogram from the reference station, filtered and without instrument response. 
    :return: angle between the horizontal components of the sensors at the two stations.
    """

    # First, calculate the c1 and c2 terms in equation A-5 of Zheng and McMechan (2006).
    c1 = 0
    c2 = 0
    for m in range(len(reference_seismogram)):  # Reference and shifted seismograms should have the same length
        try:
            c1 += ((reference_seismogram[m][0] * shifted_seismogram[m][0]) +
                   (reference_seismogram[m][1] * shifted_seismogram[m][1]))
            c2 += ((reference_seismogram[m][1] * shifted_seismogram[m][0]) -
                   (reference_seismogram[m][0] * shifted_seismogram[m][1]))
        except TypeError:  # Catch when nan values exist and skip operations for such data
            continue

    # Now, calculate the relative angle
    relative_angle = math.degrees(math.atan2(c2, c1))
    return relative_angle


def find_rotation_angle(shifted_seismogram, reference_seismogram):

    """
    Calculate the angle between the horizontal components of each sensor and those of the reference sensor.
    Uses the method of rotating the seismogram horizontal components by a given angle and cross-correlating
    the result against the reference seismogram. The rotation angle that produces the best correlation is
    taken as the rotation angle and returned to the operator.
    :param shifted_seismogram: the seismogram from the station with unknown orientation, filtered, time-shifted, and
           without instrument response.
    :param reference_seismogram: the seismogram from the reference station, filtered and without instrument response.
    :return: angle between the horizontal components of the sensors at the two stations.
    """

    normalised_xcorr_values = []
    for angle in range(0, 360):
        rad = math.radians(angle)
        north_component = []  # "north" component after rotation
        east_component = []  # "east" component after rotation

        # Spin the seismogram components clockwise by angle
        for m in range(len(shifted_seismogram[0])):
            north_component.append(math.cos(rad) * shifted_seismogram[0][m] +
                                   math.sin(rad) * shifted_seismogram[1][m])
            east_component.append(math.sin(rad) * shifted_seismogram[0][m] +
                                  math.cos(rad) * shifted_seismogram[1][m])

        # Cross-correlate the seismogram with the reference seismogram using horizontal total energy traces
        reference_horizontal_total_energy_waveform = \
            calculate_horizontal_total_energy_from_list(reference_seismogram)
        shifted_stream_horizontal_total_energy_waveform = \
            calculate_horizontal_total_energy_from_list([north_component, east_component])
        x_mean = np.nanmean(shifted_stream_horizontal_total_energy_waveform)
        y_mean = np.nanmean(reference_horizontal_total_energy_waveform)
        x_var = np.nanvar(shifted_stream_horizontal_total_energy_waveform)
        y_var = np.nanvar(reference_horizontal_total_energy_waveform)
        sum = 0
        for m in range(min(len(reference_horizontal_total_energy_waveform),
                           len(shifted_stream_horizontal_total_energy_waveform))):
            sum += ((shifted_stream_horizontal_total_energy_waveform[m] - x_mean) *
                    (reference_horizontal_total_energy_waveform[m] - y_mean))
        normalised_xcorr_value = (1 / min(len(reference_horizontal_total_energy_waveform),
                                          len(shifted_stream_horizontal_total_energy_waveform)) *
                                  sum / math.sqrt(x_var * y_var))
        normalised_xcorr_values.append(normalised_xcorr_value)
    max_xcorr_value = max(normalised_xcorr_values)
    max_xorr_value_idx = normalised_xcorr_values.index(max_xcorr_value)
    print('Rotation angle is ' + str(max_xorr_value_idx))
    return max_xorr_value_idx


if __name__ == "__main__":

    # Parse data from config file
    parameters = ['eventid_file', 'station_file', 'station_channels', 'corner_frequency', 'reference_station',
                  'reference_channels', 'local_reference', 'reference_data_dir']
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
        p_times = []
        s_times = []
        p_stations = []
        s_stations = []
        for pick in event_details.picks:
            for station in stations_to_orient:
                if '.' + station + '.' in str(pick.resource_id):
                    if pick.phase_hint == 'P':
                        p_times.append(pick.time.datetime)
                        p_stations.append(station)
                    elif pick.phase_hint == 'S':
                        s_times.append(pick.time.datetime)
                        s_stations.append(station)

        # Assume event duration does not exceed 3 minutes at any station
        start_time = min(p_times) - datetime.timedelta(seconds=15)
        end_time = start_time + datetime.timedelta(seconds=195)

        # Query waveform data from GeoNet FDSN for each station between the start and end times defined for the event
        stations = []
        streams = []
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
                                        starttime=obspy.core.utcdatetime.UTCDateTime(start_time),
                                        endtime=obspy.core.utcdatetime.UTCDateTime(end_time))
                    for tr in stream:
                        if fnmatch.fnmatch(tr.stats.channel, values[parameters.index('reference_channels')]) is True:
                            station_stream += tr
                station_stream.merge()
            elif station == values[parameters.index('reference_station')]:
                station_stream = query_fdsn(station,
                                            '??',
                                            values[parameters.index('reference_channels')],
                                            start_time.isoformat(),
                                            end_time.isoformat())[0]
            else:
                station_stream = query_fdsn(station,
                                            '??',
                                            values[parameters.index('station_channels')],
                                            start_time.isoformat(),
                                            end_time.isoformat())[0]

            # Filter the waveforms of all events to increase waveform similarly at all sensors:
            # Filter corner frequency satisfies the condition that it is much smaller than the lowest seismic velocity
            # in the propagation medium of all events scaled by the linear distance between each pair of sensors.
            station_stream.filter(type='bandpass',
                                  freqmin=1,  # In lieu of removing response, filter data to above 1 Hz so it is more
                                  # comparable between sensors of different type
                                  freqmax=float(values[parameters.index('corner_frequency')]))
            # Cut first and last 10 seconds of data to remove edge effects introduced by filtering
            station_stream.trim(station_stream[0].stats.starttime + 10,
                                station_stream[0].stats.endtime - 10)

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

        print('All data has been loaded into the script. '
              'Finding shift times to apply to each station to align them with the reference station...')

        # Align all seismograms for each event to facilitate cross-correlation:
        # Use the lag time that produces the maximum cross-correlation value between each sensor and the reference
        # sensor's vertical energy traces. In this process, convert all numpy masked arrays to numpy arrays with nan
        # values as mask fill values.
        nandices = [0, None]
        shifted_streams = streams
        for m in range(len(streams)):
            # Create a copy of the station stream and cut it to be between the P arrival time and halfway into the
            # P wave coda. Use this subset of the stream to find the time lag.
            try:
                p_time = p_times[p_stations.index(stations[m])]
            except ValueError:  # If the station doesn't have a P pick, assume it occurs at the minimum of all P picks
                p_time = min(p_times)
            try:
                s_time = s_times[s_stations.index(stations[m])]
            except ValueError:  # If the station doesn't have an S pick, assume it is 60 seconds after P pick
                s_time = p_time + datetime.timedelta(seconds=60)
            end_time = p_time + datetime.timedelta(seconds=(s_time - p_time).total_seconds() / 2)
            stream_p = streams[m].copy()
            stream_p.trim(starttime=obspy.core.utcdatetime.UTCDateTime(p_time) - 5,
                          endtime=obspy.core.utcdatetime.UTCDateTime(end_time))
            stream_p.trim(starttime=streams[m][0].stats.starttime,
                          endtime=streams[m][0].stats.endtime,
                          pad=True,
                          fill_value=None)
            # Perform the same operation to the reference stream, assuming wave arrivals are within 5 seconds at each
            ref_p = reference_station_stream.copy()
            ref_p.trim(starttime=obspy.core.utcdatetime.UTCDateTime(p_time) - 5,
                       endtime=obspy.core.utcdatetime.UTCDateTime(end_time))
            ref_p.trim(starttime=streams[m][0].stats.starttime,
                       endtime=streams[m][0].stats.endtime,
                       pad=True,
                       fill_value=None)
            print(stream_p)
            print(ref_p)
            plt.plot(ref_p[0], color='b')
            plt.plot(stream_p[0], color='r', alpha=0.8)
            plt.savefig('p.png')
            plt.clf()
            lag_time = find_lag_time(stream_p, ref_p)
            print(lag_time)
            shift_idx = int(abs(lag_time * streams[m][0].stats.sampling_rate))
            for n in range(len(streams[m])):

                # Ensure all data are masked arrays
                if not ma.is_masked(streams[m][n].data):
                    streams[m][n].data = ma.masked_array(streams[m][n].data)

                # Apply shift
                if lag_time >= 0:
                    nandices[0] = shift_idx
                    shifted_streams[m][n].data = np.asarray([float('nan')] * shift_idx +
                                                            streams[m][n].data[:-shift_idx].filled(
                                                                float('nan')).tolist())
                else:
                    nandices[1] = len(streams[m][n].data) - shift_idx
                    shifted_streams[m][n].data = np.asarray(streams[m][n].data[shift_idx:].filled(
                        float('nan')).tolist() + [float('nan')] * shift_idx)

            print('Seismograms have been aligned to the reference station by appling a shift of ' + str(lag_time) +
                  ' seconds to the ' + stations[m] + ' stream')

        for m in range(len(reference_station_stream)):
            reference_station_stream[m].data = reference_station_stream[m].data.filled(float('nan'))

        print('All seismograms have now been aligned. Searching now for the optimal cross-correlation windows...')

        plt.plot(reference_station_stream[0], color='b')
        plt.plot(shifted_streams[0][0], color='r', alpha=0.8)
        plt.savefig('shifted_plot.png')
        plt.clf()

        # Find the time window of each event for which the correlation of the vertical waveforms
        # at the sensors is highest: This window is that for which the normalised cross-correlation of the
        # total energy traces of each sensor pair is highest.
        xcorr_window_idx = []

        # Get vertical component data and calculate its envelope for reference stream
        # for tr in reference_station_stream:
            # if tr.stats.channel[-1] == 'Z':
            #     reference_horizontal_total_energy_waveform = calculate_envelope(tr.data)

        reference_horizontal_total_energy_waveform = calculate_horizontal_total_energy(reference_station_stream)
        reference_horizontal_total_energy_waveform = np.asarray(smooth_data(
            reference_horizontal_total_energy_waveform, int(values[parameters.index('corner_frequency')])))

        for m in range(len(shifted_streams)):
            # Get vertical component data and calculate its envelope for shifted stream
            # for tr in shifted_streams[m]:
            #     if tr.stats.channel[-1] == 'Z':
            #         shifted_stream_horizontal_total_energy_waveform = calculate_envelope(tr.data)

            shifted_stream_horizontal_total_energy_waveform = calculate_horizontal_total_energy(shifted_streams[m])
            shifted_stream_horizontal_total_energy_waveform = np.asarray(smooth_data(
                shifted_stream_horizontal_total_energy_waveform, int(values[parameters.index('corner_frequency')])))

            plot1 = reference_horizontal_total_energy_waveform / max(reference_horizontal_total_energy_waveform)
            plot2 = shifted_stream_horizontal_total_energy_waveform / max(shifted_stream_horizontal_total_energy_waveform)
            plt.plot(plot1, color='b')
            plt.plot(plot2, color='r', alpha=0.8)
            plt.savefig('smoothed_shifted_waveforms.png')
            plt.clf()

            # Initiate one loop to work through each possible start time in the waveform
            normalised_xcorr_values = [([0] * len(reference_horizontal_total_energy_waveform))
                                       for y in range(len(shifted_stream_horizontal_total_energy_waveform))]
            for n in range(nandices[0], len(shifted_stream_horizontal_total_energy_waveform)):
                if nandices[1] and n > nandices[1]:  # Don't do cross-correlation past the data
                    break
                # Initiate a second loop to work through each possible end time in the waveform,
                # so that all possible windows are tested, BUT require that windows are at least 1 seconds in length.
                for o in range(n + 1 * int(shifted_streams[m][0].stats.sampling_rate) + 1,
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
                    for p in range(n, o):
                        # Skip this window if there are any nan values
                        if (np.isnan(shifted_stream_horizontal_total_energy_waveform[p]) or
                                ma.is_masked(shifted_stream_horizontal_total_energy_waveform[p]) or
                                np.isnan(reference_horizontal_total_energy_waveform[p]) or
                                ma.is_masked(reference_horizontal_total_energy_waveform[p])):
                            break
                        sum += ((shifted_stream_horizontal_total_energy_waveform[p] - x_mean) *
                                (reference_horizontal_total_energy_waveform[p] - y_mean))
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

        print('All cross-correlation windows are now defined. The orientation angles will now be calculated...')

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

            plt.plot(trimmed_reference_stream[0], color='k')
            plt.plot(trimmed_and_shifted_streams[0], color='k')
            plt.savefig('trimmed_and_shifted_plot.png')
            plt.clf()

            # Remove all z component data
            for n, tr in enumerate(trimmed_and_shifted_streams):
                if tr.stats.channel[-1] == 'Z':
                    trimmed_and_shifted_streams.pop(n)
            for n, tr in enumerate(trimmed_reference_stream):
                if tr.stats.channel[-1] == 'Z':
                    trimmed_reference_stream.pop(n)

            # Calculate angle
            # orientation_angle = calculate_angle(trimmed_and_shifted_streams, trimmed_reference_stream)
            orientation_angle = find_rotation_angle(trimmed_and_shifted_streams, trimmed_reference_stream)
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
