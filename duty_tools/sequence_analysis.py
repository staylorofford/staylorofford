# -*- coding: utf-8 -*-
"""
Parse data from GHA Kermadec Spreadsheet or csv of eventIDs and produce plots of data/analysis
"""

import argparse
import datetime
from io import BytesIO
import math
import matplotlib
import matplotlib.pyplot as plt
from obspy.io.quakeml.core import Unpickler
import pycurl
import subprocess

quakeml_reader = Unpickler()

def get_event(eventID):

    """
    Get a single event from the GeoNet FDSN service
    :param eventID: GeoNet eventID
    :return: FDSN event
    """

    query = "https://service.geonet.org.nz/fdsnws/event/1/query?eventID=" + eventID
    queryresult = curl(query)
    event = quakeml_reader.loads(queryresult)
    return event


def curl(curlstr):

    """
    Perform curl with curlstr
    :param curlstr: string to curl
    :return: curl output
    """

    buffer = BytesIO()
    c = pycurl.Curl()
    c.setopt(c.URL, curlstr)
    c.setopt(c.WRITEDATA, buffer)
    c.perform()
    c.close()

    return buffer.getvalue()


def parse_spreadsheet_csv(file_path, data_columns):

    """
    Parse contents of sequence csv created by GHAs into data lists
    :param file_path: path to sequence csv file
    :return: lists containing relevant data from sequence csv
    """

    data_lists = [[] for i in range(len(data_columns))]
    with open(file_path, 'r') as openfile:
        for row in openfile:
            cols = row.split(',')

            try:  # Determine if the row is a header or not
                data_lists[0].append(datetime.datetime.strptime(cols[0] + cols[1], '%Y/%m/%d%H:%M:%S'))
            except:
                continue

            if cols[2] == 'Y':
                data_lists[1].append(True)
            elif cols[2] == 'N':
                data_lists[1].append(False)

            if cols[3] == 'Y':
                data_lists[2].append(True)
            else:
                data_lists[2].append(False)

            data_lists[3].append(cols[4])

            try:
                data_lists[4].append(float(cols[5]))
            except:
                data_lists[4].append(float('nan'))

            try:
                data_lists[5].append(float(cols[6]))
            except:
                data_lists[5].append(float('nan'))

            if cols[7] == 'Y':
                data_lists[6].append(True)
            elif cols[7] == 'N':
                data_lists[6].append(False)

            try:
                data_lists[7].append(float(cols[8]))
            except:
                data_lists[7].append(float('nan'))

    return data_lists


def get_eventID_data(eventIDs, data_columns, event_fields):

    """
    Gather data for each event in the eventID list
    """

    data_lists = [[] for i in range(len(data_columns))]
    event_details = [[] for i in range(len(event_fields))]
    for eventID in eventIDs:
        event = get_event(eventID)[0]

        # Parse data into data_lists format

        data_lists[0].append(event.origins[0].time.datetime)
        data_lists[1].append(True)  # can't tell from FDSN, so just say Yes
        data_lists[2].append(True)  # these are all events
        data_lists[3].append(eventID)
        MLv = None
        mB = None
        for magnitude in event.magnitudes:
            if magnitude.magnitude_type == 'MLv':
                MLv = True
                data_lists[4].append(float(magnitude.mag))
            elif magnitude.magnitude_type == 'mB':
                mB = True
                data_lists[5].append(float(magnitude.mag))
        if not MLv:
            data_lists[4].append(float('nan'))
        if not mB:
            data_lists[5].append(float('nan'))
        data_lists[6].append(False)  # Do not look for USGS matches
        data_lists[7].append(float('nan'))  # As above

        # Parse data into event_details format

        event_details[0].append(event.resource_id)
        event_details[1].append(event.origins[0].time)
        event_details[4].append(event.origins[0].longitude)
        event_details[5].append(event.origins[0].latitude)
        event_details[6].append(event.origins[0].depth / 1000.0)
        MLv = False
        mB = False
        for magnitude in event.magnitudes:
            if magnitude.magnitude_type == 'MLv':
                MLv = True
                event_details[2].append(magnitude.mag)
            elif magnitude.magnitude_type == 'mB':
                mB = True
                event_details[3].append(magnitude.mag)
        if not MLv:
            event_details[2].append(float('nan'))
        if not mB:
            event_details[3].append(float('nan'))

    return data_lists, event_details


def get_eventIDs(parameter_file):

    """
    Parse parameter values from parameter file, then query FDSN for eventIDs satisfying those parameters
    :param parameter_file: text file containing parameter=value pairs for the earthquake search
    :return: a list of eventIDs
    """

    # Get parameter values from file

    parameters = ['starttime', 'endtime', 'minlatitude', 'maxlatitude', 'minlongitude', 'maxlongitude', 'mindepth',
                  'maxdepth', 'minmagnitude', 'maxmagnitude']
    values = ['nan'] * len(parameters)
    with open(parameter_file, 'r') as openfile:
        for row in openfile:
            rc = parameters.index(row.split('=')[0])
            values[rc] = row.split('=')[1][:-1]

    # Build GeoNet FDSN Query

    query = 'https://service.geonet.org.nz/fdsnws/event/1/query?'
    print('\nQuery parameters are:')
    for n in range(len(parameters)):
        print(parameters[n] + '=' + values[n])
        if values[n] != 'nan':
            query += parameters[n] + '=' + values[n] + '&'
    query = query[:-1]

    # Query FDSN

    print('\nQuerying...')
    queryresult = curl(query)
    try:
        events = quakeml_reader.loads(queryresult)
    except:
        print('There is an issue with your parameter values, please check them over.')
        exit()

    # Strip eventIDs from query
    # NOTE: this is inefficient programming! All data that is later queried is already contained in this query result!

    eventIDs = []
    for event in events:
        eventIDs.append(str(event.resource_id).split('/')[-1])
    print(str(len(eventIDs)) + ' events found within query parameters\n')

    return eventIDs


def calculate_distribution(data_list, round_to=0.1):

    """
    Calculate the distribution of a list of floats, where floats are rounded to round_to decimal places.
    :param data_list: list of data values
    :param round_to: value to round data to for distribution calculation
    :return: data values in distribution, number of data at each data value
    """

    nn_data_list = []
    for i in range(len(data_list)):
        try:
            if round_to <= 1:
                nn_data_list.append(int(round((1 / round_to) * data_list[i])))
            else:
                nn_data_list.append(int(round(data_list[i] / round_to)))
        except:  # Fails on nan values
            pass

    D = [0] * int(max(nn_data_list) + 1)
    data_values = range(max(nn_data_list) + 1)

    for i in data_values:
        if nn_data_list.count(i) > 0:
            D[i] = nn_data_list.count(i)

    data_values = list(data_values)

    for i in range(len(data_values)):
        if round_to <= 1:
            data_values[i] /= (1 / round_to)

    return data_values, D


def distance(p1, p2):

    """

    Calculate the distance between two points on a sphere, assuming no elevation difference exists between the points.

    :param p1: point 1 horizontal position [lon, lat]
    :param p2: point 2 horizontal position [lon, lat]
    :return: distance between p1 and p2
    """

    # Convert longitude and latitude to radians

    for n in range(len(p1)):
        p1[n] = math.radians(p1[n])
        p2[n] = math.radians(p2[n])

    # Calculate distance

    r = 6371000
    d = 2 * r * math.asin(math.sqrt((math.sin((p2[1] - p1[1]) / 2.0)) ** 2 + math.cos(p1[0]) * math.cos(p2[0]) *
                          (math.sin((p2[0] - p1[0]) / 2.0) ** 2)))
    return d / 1000


# Set data column details (for data_lists)

data_columns_in = ['date', 'GLKZ arrival time (UTC)', 'autolocation?', 'GHA event created?', 'eventID',
                   'MLv', 'mB', 'USGS location?', 'USGS Mww', 'Comments']
data_columns = data_columns_in[1:9]

# Set event details (for event_details)

event_fields = ['eventID', 'origin_time', 'MLv', 'mB', 'longitude', 'latitude', 'depth',
                'relative_origin_time', 'relative_origin_distance']

# Set plot details for use with magnitudes

colors = ['k', 'r', 'b']
labels = ['all detections', 'manual detections', 'automatic detections']
magnitude_indices = [4, 5]

# Parse arguments

parser = argparse.ArgumentParser()
parser.add_argument('--params', type=str, help='Path to text file containing earthquake search parameter=value pairs')
parser.add_argument('--csv', type=str, help='Path to csv file containing eventIDs in sequence')
parser.add_argument('--xlsx', type=str, help='Path to excel spreadsheet containing all data of sequence')
parser.add_argument('--t0', type=str, help='Optional: time to plot data relative to in ISO8601 format, default is'
                                           'time of largest earthquake in data')
args = parser.parse_args()

# Parse data

if args.xlsx:

    # If the input is a spreadsheet, convert to a csv and then parse data in from it

    p = subprocess.Popen('ssconvert ' + args.xlsx + ' sequence_analysis.csv', shell=True)
    p.wait()
    file_path = 'sequence_analysis.csv'
    data_lists = parse_spreadsheet_csv(file_path, data_columns)

    # Get event details from GeoNet FDSN for events in the spreadsheet

    event_details = [[] for i in range(len(event_fields))]

    for eventID in data_lists[3]:
        try:
            event = get_event(eventID)[0]
        except:  # Fails on N/A value
            continue
        event_details[0].append(event.resource_id)
        event_details[1].append(event.origins[0].time)
        event_details[4].append(event.origins[0].longitude)
        event_details[5].append(event.origins[0].latitude)
        event_details[6].append(event.origins[0].depth / 1000.0)
        for magnitude in event.magnitudes:
            if magnitude.magnitude_type == 'MLv':
                event_details[2].append(magnitude.mag)
            elif magnitude.magnitude_type == 'mB':
                event_details[3].append(magnitude.mag)

elif args.csv:

    # If the input is a csv file containing a list of eventIDs, get the relevant data for these events

    eventIDs = []
    with open(args.csv, 'r') as openfile:
        for row in openfile:
            eventIDs.append(row[:row.index('\n')])

    data_lists, event_details = get_eventID_data(eventIDs, data_columns, event_fields)

elif args.params:

    # If the input is a parameter file, find the relevant eventIDs and gather data for these events

    eventIDs = get_eventIDs(args.params)
    data_lists, event_details = get_eventID_data(eventIDs, data_columns, event_fields)

else:
    print('You need to give an input argument, use -h flag to see which options are available.')
    exit()

# Ensure data are sorted by time

[data_lists[0], data_lists[1], data_lists[2], data_lists[3],
data_lists[4], data_lists[5], data_lists[6], data_lists[7]] = zip(*sorted(zip(data_lists[0], data_lists[1],
                                                                              data_lists[2], data_lists[3],
                                                                              data_lists[4], data_lists[5],
                                                                              data_lists[6], data_lists[7])))

# Get details of largest event (MLv) in sequence

max_mag = 0
for i in range(len(data_lists[magnitude_indices[0]])):
    if data_lists[magnitude_indices[0]][i] > max_mag:
        max_mag, max_time, max_n = data_lists[magnitude_indices[0]][i], data_lists[0][i], i

# Set relative time to plot data against

tn = datetime.datetime.utcnow()  # get current time
if args.t0:
    t0 = datetime.datetime.strptime(args.t0, '%Y-%m-%dT%H:%M:%SZ')
else:
    t0 = max_time
tn_rel = (tn - t0).total_seconds() / 86400.0  # get relative current time

# Calculate time, distance of events from largest event in sequence

times, distances = [], []
for i in range(len(event_details[0])):
    event_details[7].append((event_details[1][i].datetime - t0).total_seconds() / 86400)
    event_details[8].append(distance([event_details[4][i], event_details[5][i]],
                              [event_details[4][max_n], event_details[5][max_n]]))

# Now that the event_details lists are full, ensure the data lists are sorted by time

[event_details[1], event_details[0], event_details[2], event_details[3],
event_details[4], event_details[5], event_details[6]] = zip(*sorted(zip(event_details[1], event_details[0],
                                                                        event_details[2], event_details[3],
                                                                        event_details[4], event_details[5],
                                                                        event_details[6])))

# Build automatic and manual detection data lists

automatic_data_lists = [[] for i in range(len(data_columns))]
manual_data_lists = [[] for i in range(len(data_columns))]
dll = [data_lists, manual_data_lists, automatic_data_lists]

indices = [False, True]
for i in range(len(data_lists[0])):
    try:
        idx = indices.index(data_lists[1][i])
    except:  # Fails for "non-locatable" events (those that have detection data missing)
        idx = 1
    try:
        for j in range(len(data_lists)):
            dll[idx + 1][j].append(data_lists[j][i])
    except:  # This fails when trying to get detection data that doesn't exist, but only after getting time data
        continue

# Plot event locations in 2D: needs work to show a third variable: mag, depth, time? Perhaps 3 subplots...

fig = plt.figure(figsize=(9, 6))
cmap = plt.cm.get_cmap('jet')
EQ_loc_data = [event_details[7], event_details[2], event_details[3]]
EQ_loc_labels = ['relative time', 'MLv', 'mB']
for j in range(len(EQ_loc_data)):
    if j == 0:
        ax1 = plt.subplot(str(len(EQ_loc_data)) + '1' + str(j + 1))
    else:
        ax2 = plt.subplot(str(len(EQ_loc_data)) + '1' + str(j + 1), sharex=ax1)

    plot = plt.scatter(event_details[4], event_details[5], s=9, c=EQ_loc_data[j], cmap=cmap)

    plt.xlabel('longitude', labelpad=10)
    plt.ylabel('latitude', labelpad=10)

    cbar = plt.colorbar(plot)
    cbar.ax.get_yaxis().labelpad=15
    cbar.ax.set_ylabel(EQ_loc_labels[j], rotation=270)

# Only have x-axis tick labels on last subplot

for ax in plt.gcf().axes:
    try:
        ax.label_outer()
    except:
        pass

plt.subplots_adjust(left=0.13, right=0.98, bottom=0.12, top=0.9)
fig.suptitle('Earthquake Locations, t0=' + str(t0), y=0.95)
plt.show()

# Plot event location details over time

fig = plt.figure(figsize=(9, 6))
EQ_loc_data = [event_details[4], event_details[5], event_details[6]]
EQ_loc_labels = ['longitude', 'latitude', 'depth']
for j in range(len(EQ_loc_data)):
    if j == 0:
        ax1 = plt.subplot(str(len(EQ_loc_data)) + '1' + str(j + 1))
    else:
        ax2 = plt.subplot(str(len(EQ_loc_data)) + '1' + str(j + 1), sharex=ax1)

    plt.scatter(event_details[7], EQ_loc_data[j], color='k', s=9)
    plt.axvline(0, color='k', linestyle='--', alpha=0.5)
    plt.axvline(tn_rel, color='r', alpha=0.5)
    plt.xlabel('days from t0', labelpad=10)
    plt.ylabel(EQ_loc_labels[j], labelpad=10)

# Only have x-axis tick labels on last subplot

for ax in plt.gcf().axes:
    try:
        ax.label_outer()
    except:
        pass

plt.subplots_adjust(left=0.13, right=0.98, bottom=0.12, top=0.9)
fig.suptitle('Earthquake Locations over Time, t0=' + str(t0), y=0.95)
plt.show()

# Plot relative distance of events over time

fig = plt.figure(figsize=(9, 6))
plt.scatter(event_details[7], event_details[8], color='k', s=9)
plt.axvline(0, color='k', linestyle='--', alpha=0.5)
plt.axvline(tn_rel, color='r', alpha=0.5)
plt.xlabel('days from t0', labelpad=10)
plt.ylabel('distance from largest earthquake (km)', labelpad=10)
plt.subplots_adjust(left=0.13, right=0.98, bottom=0.12, top=0.9)
fig.suptitle('Relative Time and Distance of Earthquakes in Sequence, t0=' + str(t0), y=0.95)
plt.show()

# Plot cumulative event number and cumulative moment release over time

fig = plt.figure(figsize=(9, 6))

ax1 = plt.subplot('211')
for i in range(len(dll)):
    if len(dll[i][0]) != 0:
        plt.plot(dll[i][0], list(range(len(dll[i][0]))), color=colors[i], label=labels[i])
plt.legend()
plt.axvline(t0, color='k', linestyle='--', alpha=0.5)
plt.axvline(tn, color='r', linestyle='-', alpha=0.5)
plt.ylabel('cumulative\nevent count', labelpad=10)

ax2 = plt.subplot('212', sharex=ax1)
for i in range(len(dll)):
    if len(dll[i][0]) != 0:
        cumulative_moments, mB_times = [], []
        for j in range(len(dll[i][magnitude_indices[1]])):  # Calculate moment for each event
            if str(dll[i][magnitude_indices[1]][j]) != 'nan':
                try:
                    cumulative_moments.append(cumulative_moments[-1] + 10 ** (1.5 * (dll[i][magnitude_indices[1]][j] + 6.06)))
                except:
                    cumulative_moments.append(10 ** (1.5 * (dll[i][magnitude_indices[1]][j] + 6.06)))
                mB_times.append(dll[i][0][j])
        plt.plot(mB_times, cumulative_moments, color=colors[i], label=labels[i])
        if i == 0:  # Give the cumulative moment magnitude
            try:
                print('Cumulative Mw is ' + str((2 / 3) * math.log10(cumulative_moments[-1]) - 6.06))
            except:  # Fails when there are no mB measurements
                pass

plt.axvline(t0, color='k', linestyle='--', alpha=0.5)
plt.axvline(tn, color='r', alpha=0.5)
plt.ylabel('cumulative\nmoment release', labelpad=10)
plt.xlabel('time (UTC)', labelpad=10)
plt.subplots_adjust(left=0.13, right=0.98, bottom=0.12, top=0.9)
fig.suptitle('Earthquake Cumulative Time Series', y=0.95)
fig.autofmt_xdate()
plt.show()

# Plot magnitude time series for all detected events

fig = plt.figure(figsize=(9, 6))
for j in range(len(magnitude_indices)):
    if j == 0:
        ax1 = plt.subplot(str(len(magnitude_indices)) + '1' + str(j + 1))
    else:
        ax2 = plt.subplot(str(len(magnitude_indices)) + '1' + str(j + 1), sharex=ax1)

    for k in range(len(dll[0][magnitude_indices[j]])):
        if str(dll[0][magnitude_indices[j]][k]) != 'nan':
            plt.scatter(dll[0][0][k], dll[0][magnitude_indices[j]][k], color='k', s=3)
    plt.axvline(t0, color='k', linestyle='--', alpha=0.5)
    plt.axvline(tn, color='r', alpha=0.5)
    plt.ylabel(data_columns[magnitude_indices[j]], labelpad=10)

plt.xlabel('time (UTC)', labelpad=10)
plt.xlim(min(dll[0][0]), max(dll[0][0]))
plt.subplots_adjust(left=0.13, right=0.98, bottom=0.12, top=0.9)
fig.suptitle('Earthquake Magnitude Over Time', y=0.95)
fig.autofmt_xdate()
plt.show()

# Calculate magnitude distribution and plot it

fig = plt.figure(figsize=(9, 6))
for j in range(len(magnitude_indices)):

    if j == 0:
        ax1 = plt.subplot(str(len(magnitude_indices)) + '1' + str(j + 1))
    else:
        ax2 = plt.subplot(str(len(magnitude_indices)) + '1' + str(j + 1), sharex=ax1)

    for i in [0, 1]:

        if len(dll[i][0]) != 0:

            # On both axes, plot all detections, and overlay all manual detections

            try:
                mag_values, value_count = calculate_distribution(dll[i][magnitude_indices[j]])
                plt.bar(mag_values, value_count, width=0.05, color=colors[0], label=labels[0])
            except:
                pass

    plt.xlabel(data_columns[magnitude_indices[j]], labelpad=10)
    plt.ylabel('count', labelpad=10)
    try:
        plt.legend()
    except:
        pass

plt.subplots_adjust(left=0.13, right=0.98, bottom=0.12, top=0.9, hspace=0.35)
fig.suptitle('Earthquake Magnitude Distribution', y=0.95)
plt.show()

# Plot interevent time distribution

interevent_times = []
for i in range(1, len(dll[0][0])):
    interevent_times.append((dll[0][0][i] - dll[0][0][i - 1]).total_seconds())
interevent_time_values, it_counts = calculate_distribution(interevent_times, 3600)

fig = plt.figure(figsize=(9, 6))
plt.bar(interevent_time_values, it_counts, color='k')
plt.ylabel('count', labelpad=10)
plt.xlabel('interevent time', labelpad=10)
plt.subplots_adjust(left=0.13, right=0.98, bottom=0.12, top=0.9, hspace=0.35)
fig.suptitle('Earthquake Inter-Event Time Distribution', y=0.95)
plt.show()

# Plot magnitudes diff. USGS Mww

fig = plt.figure(figsize=(9, 6))
for j in range(len(magnitude_indices)):
    if j == 0:
        ax1 = plt.subplot(str(len(magnitude_indices)) + '1' + str(j + 1))
    else:
        ax2 = plt.subplot(str(len(magnitude_indices)) + '1' + str(j + 1), sharex=ax1)

    for k in range(len(dll[0][magnitude_indices[j]])):
        if str(dll[0][magnitude_indices[j]][k]) != 'nan' and str(dll[0][-1][k]) != 'nan':
            diff = dll[0][magnitude_indices[j]][k] - dll[0][-1][k]
            plt.scatter(dll[0][magnitude_indices[j]][k], diff, color='k', s=9)
    plt.ylabel(data_columns[magnitude_indices[j]] + ' - USGS Mww', labelpad=10)
    plt.xlabel(data_columns[magnitude_indices[j]], labelpad=10)
    plt.ylim(-1, 1)

plt.subplots_adjust(left=0.13, right=0.98, bottom=0.12, top=0.9, hspace=0.35)
fig.suptitle('Earthquake Magnitude - Mww', y=0.95)
plt.show()