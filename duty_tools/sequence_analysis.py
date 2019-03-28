# -*- coding: utf-8 -*-
"""
Parse data from csv-converted GHA Kermadec Spreadsheet and produce plots of analysis
"""

import argparse
import datetime
from io import BytesIO
import math
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


def parse_csv(file_path, data_columns_out):

    """
    Parse contents of sequence csv created by GHAs into data lists
    :param file_path: path to sequence csv file
    :return: lists containing relevant data from sequence csv
    """

    data_lists = [[] for i in range(len(data_columns_out))]
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


# Data column details

data_columns_in = ['date', 'GLKZ arrival time (UTC)', 'autolocation?', 'GHA event created?', 'eventID',
                   'MLv', 'mB', 'USGS location?', 'USGS Mww', 'Comments']
data_columns_out = data_columns_in[1:9]

# Parse arguments

parser = argparse.ArgumentParser()
parser.add_argument('-csv', type=str, help='Path to csv-version of Kermadec sequence spreadsheet')
parser.add_argument('-xlsx', type=str, help='Path to excel spreadsheet of Kermadec sequence')
args = parser.parse_args()

file_path = args.csv

# If the input is a spreadsheet, convert to a csv
if not file_path:
    if args.xlsx:
        p = subprocess.Popen('ssconvert ' + args.xlsx + ' sequence_analysis.csv', shell=True)
        p.wait()
        file_path = 'sequence_analysis.csv'
    else:
        print('You need to give an argument, use -h flag to see which options are available.')
        exit()

# Parse data from csv file

data_lists = parse_csv(file_path, data_columns_out)

# Build automatic and manual detection data lists

automatic_data_lists = [[] for i in range(len(data_columns_out))]
manual_data_lists = [[] for i in range(len(data_columns_out))]
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

# Calculate magnitude distribution and plot it

colors = ['g', 'r', 'b']
labels = ['all detections', 'manual detections', 'automatic detections']
magnitude_indices = [4, 5]

fig = plt.figure(figsize=(9, 6))
mag_values, value_count = [[] for i in range(len(magnitude_indices))], [[] for i in range(len(magnitude_indices))]
for j in range(len(magnitude_indices)):

    if j == 0:
        ax1 = plt.subplot(str(len(magnitude_indices)) + '1' + str(j + 1))
    else:
        ax2 = plt.subplot(str(len(magnitude_indices)) + '1' + str(j + 1), sharex=ax1)

    for i in range(1, len(dll)):
        mag_values[i - 1], value_count[i - 1] = calculate_distribution(dll[i][magnitude_indices[j]])

    for n in range(len(mag_values[0])):
        for m in range(len(mag_values[1])):
            if mag_values[0][n] == mag_values[1][m]:
                value_count[1][m] += value_count[0][n]

    for i in range(1, len(dll)):
        plt.bar(mag_values[len(dll) - (i + 1)], value_count[len(dll) - (i + 1)], width=0.05,
                color=colors[len(dll) - i], label=labels[len(dll) - i])

    plt.xlabel(data_columns_out[magnitude_indices[j]], labelpad=10)
    plt.ylabel('count', labelpad=10)
    if j == 0:
        plt.legend()

plt.subplots_adjust(left=0.13, right=0.98, bottom=0.12, top=0.9, hspace=0.35)
fig.suptitle('Earthquake Magnitude Distribution', y=0.95)
plt.show()

# Get details of largest event (mB) in sequence

max_mag = 0
for i in range(len(data_lists[magnitude_indices[1]])):
    if data_lists[magnitude_indices[1]][i] > max_mag:
        max_mag, max_time, max_n = data_lists[magnitude_indices[1]][i], data_lists[0][i], i

# Plot cumulative event number and cumulative moment release over time

fig = plt.figure(figsize=(9, 6))

ax1 = plt.subplot('211')
for i in range(len(dll)):
    plt.plot(dll[i][0], list(range(len(dll[i][0]))), color=colors[i], label=labels[i])
plt.legend()
plt.axvline(max_time, color='k', linestyle='--', alpha=0.5)
plt.ylabel('cumulative\nevent count', labelpad=10)

ax2 = plt.subplot('212', sharex=ax1)
for i in range(len(dll)):
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
        print('Cumulative Mw is ' + str((2 / 3) * math.log10(cumulative_moments[-1]) - 6.06))

plt.axvline(max_time, color='k', linestyle='--', alpha=0.5)
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
    plt.axvline(max_time, color='k', linestyle='--', alpha=0.5)
    plt.ylabel(data_columns_out[magnitude_indices[j]], labelpad=10)

plt.xlabel('time (UTC)', labelpad=10)
plt.xlim(min(dll[0][0]), max(dll[0][0]))
plt.subplots_adjust(left=0.13, right=0.98, bottom=0.12, top=0.9)
fig.suptitle('Earthquake Magnitude Over Time', y=0.95)
fig.autofmt_xdate()
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
    plt.ylabel(data_columns_out[magnitude_indices[j]] + ' - USGS Mww', labelpad=10)
    plt.xlabel(data_columns_out[magnitude_indices[j]], labelpad=10)
    plt.ylim(-1, 1)

plt.subplots_adjust(left=0.13, right=0.98, bottom=0.12, top=0.9, hspace=0.35)
fig.suptitle('Earthquake Magnitude - Mww', y=0.95)
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

# Get event details from GeoNet FDSN

event_fields = ['eventID', 'origin_time', 'MLv', 'mB', 'longitude', 'latitude', 'depth',
                'relative_origin_time', 'relative_origin_distance']
event_details = [[] for i in range(len(event_fields))]

for eventID in dll[0][3]:
    try:
        event = get_event(eventID)[0]
    except:  # Fails on N/A value
        continue
    event_details[0].append(event.resource_id)
    event_details[1].append(event.origins[0].time)
    event_details[4].append(event.origins[0].longitude)
    event_details[5].append(event.origins[0].latitude)
    event_details[6].append(event.origins[0].depth)
    for magnitude in event.magnitudes:
        if magnitude.magnitude_type == 'MLv':
            event_details[2].append(magnitude.mag)
        elif magnitude.magnitude_type == 'mB':
            event_details[3].append(magnitude.mag)

# Calculate time, distance of events from largest event in sequence

times, distances = [], []
for i in range(len(event_details[0])):
    event_details[7].append((event_details[1][i].datetime - event_details[1][max_n].datetime).total_seconds() / 86400)
    event_details[8].append(distance([event_details[4][i], event_details[5][i]],
                              [event_details[4][max_n], event_details[5][max_n]]))

# Plot event details

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
    plt.xlabel('days since largest earthquake', labelpad=10)
    plt.ylabel(EQ_loc_labels[j], labelpad=10)

# Only have x-axis tick labels on last subplot

for ax in plt.gcf().axes:
    try:
        ax.label_outer()
    except:
        pass

plt.subplots_adjust(left=0.13, right=0.98, bottom=0.12, top=0.9)
fig.suptitle('Earthquake Locations over Time', y=0.95)
plt.show()

# Plot relative distance of events over time

fig = plt.figure(figsize=(9, 6))
plt.scatter(event_details[7], event_details[8], color='k', s=9)
plt.axvline(0, color='k', linestyle='--', alpha=0.5)
plt.xlabel('days since largest earthquake', labelpad=10)
plt.ylabel('distance from largest earthquake (km)', labelpad=10)
plt.subplots_adjust(left=0.13, right=0.98, bottom=0.12, top=0.9)
fig.suptitle('Relative Time and Distance of Earthquakes in Sequence', y=0.95)
plt.show()