from obspy.io.quakeml.core import Unpickler
import pycurl
from io import BytesIO
import math
import matplotlib.pyplot as plt
import datetime
import glob

quakeml_reader = Unpickler()

def FDSN_event_query(service, minmagnitude, minlongitude, maxlongitude,
                     minlatitude, maxlatitude, starttime = '0000-01-01T00:00:00',
                     maxmagnitude = 10):
    """
    Use obspy with pycurl to query event catalogs from FDSN members.

    :param service:
    :param minmagnitude:
    :param maxmagnitude:
    :param minlongitude:
    :param maxlongitude:
    :param minlatitude:
    :param maxlatitude:
    :param starttime:
    :return:
    """

    # Adjust format if required

    if "usgs" in service:
        maxlongitude += 360

    # Curl FDSN response and parse quakeml bytes string through obspy
    # using an interative query approach when a single query fails


    factor = 1
    success = False
    while success == False:

        successes = 0
        events = []

        # Build magnitude ranges for query

        magnitude_limits = [minmagnitude]
        for i in range(1, factor + 1):
                magnitude_limits.append(magnitude_limits[-1] + (maxmagnitude - minmagnitude) / factor)

        # Run queries

        for i in range(1, len(magnitude_limits)):

            # Build query

            query = ""
            query = query.join((service, "query?", "minmagnitude=", str(magnitude_limits[i - 1]),
                                "&maxmagnitude=", str(magnitude_limits[i]),
                                "&minlatitude=", str(minlatitude), "&maxlatitude=", str(maxlatitude),
                                "&minlongitude=", str(minlongitude), "&maxlongitude=", str(maxlongitude),
                                "&starttime=", starttime))

            try:

                print('\nAttempting FDSN catalog query for events between M ' + str(magnitude_limits[i - 1]) + ' and ' + str(magnitude_limits[i]))
                queryresult = curl(query)

                print('Processing query result...')
                process_start = datetime.datetime.now()
                catalog = quakeml_reader.loads(queryresult)
                process_end = datetime.datetime.now()
                print('That took ' + str((process_end - process_start).total_seconds()) + ' seconds')

                events.extend(catalog.events)
                print('Catalog now has ' + str(len(events)) + ' events')
                successes += 1

            except:

                print('Failed!')
                process_end = datetime.datetime.now()
                print('That took ' + str((process_end - process_start).total_seconds()) + ' seconds')
                if successes > 0:
                    print('Assuming query failed because no events exist at high magnitude range')
                    successes += 1
                else:
                    factor += 999 # Only fails for huge datasets, so try minimise the size of the first new query
                    break

        if successes == len(magnitude_limits) - 1:
            success = True

    return events


def curl(curlstr):

    """
    Perform curl with curlstr
    :param curlstr: string to curl
    :return: curl output
    """

    print('Querying...')
    process_start = datetime.datetime.now()
    buffer = BytesIO()
    c = pycurl.Curl()
    c.setopt(c.URL, curlstr)
    c.setopt(c.WRITEDATA, buffer)
    c.perform()
    c.close()
    process_end = datetime.datetime.now()
    print('That took ' + str((process_end - process_start).total_seconds()) + ' seconds')

    return buffer.getvalue()


def to_cartesian(latitude, longitude, depth):
    """
    Convert a point on the Earth in spherical coordinates to cartesian coordinates
    :param latitude: point latitude in decimal degrees, south is negative
    :param longitude: point longitude in decimal degrees, west is negative
    :param depth: point depth in metres, down is positive
    :return: cartesian coordinates of the point
    """

    # Ensure variables are floats

    latitude = float(latitude)
    longitude = float(longitude)
    depth = float(depth)

    Re = 6371000  # average Earth radius in m (assuming a spherical Earth)

    # Convert to spherical coordinate conventions

    r = Re - depth
    if latitude < 0:
        latitude += 180

    # Convert to cartesian coordinates

    x = r * math.sin(latitude) * math.cos(longitude)
    y = r * math.sin(latitude) * math.sin(longitude)
    z = r * math.cos(latitude)

    return x, y, z


def match_magnitudes(magnitude_timeseries, timeseries_types, comparison_magnitudes, max_dt, max_dist):

    """
    Match events between two catalogs and save all matched events with magnitude
    types in comparison_magnitudes to file
    :param event_catalog: list of obspy event objects
    :param reference_catalog: list of obspy event objects to match against
    :param comparison_magnitudes: list of two lists containing all magnitudes to match
                                between the catalog and the reference catalog
    :param max_dt: maximum seconds (absolute) between events for them to be matched
    :param max_dist: maximum distance (absolute) between events for them to be matched
    :return: saves matched events with matched magnitudes to a csv file
    """

    datalist = [[[] for i in range(len(timeseries_types))] for j in range(7)]
    for n in range(len(timeseries_types)):
        if timeseries_types[n].split('_')[0] == catalog_names[0].split('_')[0] and \
                timeseries_types[n].split('_')[2] in comparison_magnitudes[0]:
            for m in range(len(timeseries_types)):
                if timeseries_types[m].split('_')[0] == catalog_names[1].split('_')[0] and \
                        timeseries_types[m].split('_')[2] in comparison_magnitudes[1]:
                    for i in range(len(magnitude_timeseries[0][n])):
                        lengths = [[], []]
                        ETi, ELa, ELo, EDe = [datetime.datetime.strptime(magnitude_timeseries[1][n][i], '%Y-%m-%dT%H:%M:%S.%fZ'),
                                              float(magnitude_timeseries[4][n][i]), float(magnitude_timeseries[5][n][i]),
                                              float(magnitude_timeseries[6][n][i])]
                        Ex, Ey, Ez = to_cartesian(ELa, ELo, EDe)
                        for j in range(len(magnitude_timeseries[0][m])):
                            RETi, RELa, RELo, REDe = [datetime.datetime.strptime(magnitude_timeseries[1][m][j], '%Y-%m-%dT%H:%M:%S.%fZ'),
                                              float(magnitude_timeseries[4][m][j]), float(magnitude_timeseries[5][m][j]),
                                              float(magnitude_timeseries[6][m][j])]
                            REx, REy, REz = to_cartesian(RELa, RELo, REDe)

                            temporal_length = abs((ETi - RETi).total_seconds())

                            if temporal_length > max_dt:
                                continue

                            spatial_length = math.sqrt((Ex - REx) ** 2 + (Ey - REy) ** 2 + (Ez - REz) ** 2) / 1000.0
                            if spatial_length > max_dist:
                                continue

                            lengths[0].append(math.sqrt(temporal_length ** 2 + spatial_length ** 2))
                            lengths[1].append(n)

                        print(str(len(lengths[0])) + ' reference events found for the ' + str(
                            i) + 'th event in the catalog')
                        if len(lengths[0]) > 0:
                            match_idx = lengths[1][lengths[0].index(min(lengths[0]))]
                            datalist[0][n].append(magnitude_timeseries[0][n][i])
                            datalist[1][n].append(magnitude_timeseries[3][n][i])
                            datalist[2][n].append(magnitude_timeseries[3][m][match_idx])
                            datalist[3][n].append(lengths[0][lengths[0].index(min(lengths[0]))])
                            datalist[4][n].append(magnitude_timeseries[4][n][i])
                            datalist[5][n].append(magnitude_timeseries[5][n][i])
                            datalist[6][n].append(magnitude_timeseries[6][n][i])

                    print('    ' + str(
                        len(datalist[0][n])) + ' matched events with desired magnitude types were found')

                    # Save data

                    print('\nSaving data to file...')
                    if len(datalist[0][n]) == 0:
                        continue
                    else:
                        with open(timeseries_types[n].split('_')[2] + '_' + timeseries_types[m].split('_')[
                            2] + '_' + 'magnitude_matches.csv', 'w') as outfile:
                            outfile.write(
                                'eventID,' + timeseries_types[n].split('_')[2] + ',' + timeseries_types[m].split('_')[
                                    2] + ',' +
                                'length,longitude,latitude,depth' + '\n')
                        with open(timeseries_types[n].split('_')[2] + '_' + timeseries_types[m].split('_')[
                            2] + '_' + 'magnitude_matches.csv', 'a') as outfile:
                            for i in range(len(datalist[0][n])):
                                outfile.write(
                                    str(datalist[0][n][i]) + ',' + str(datalist[1][n][i]) + ',' + str(datalist[2][n][i])
                                    + ',' + str(datalist[3][n][i]) + ',' + str(datalist[4][n][i]) + ',' + str(
                                        datalist[5][n][i])
                                    + ',' + str(datalist[6][n][i]))


def build_magnitude_timeseries(catalogs, catalog_names, comparison_magnitudes):

    """
    Saves magnitude timeseries to disk for each magnitude type in the comparison_magnitudes
    list corresponding to a given catalog
    :param catalogs: list of lists containing obspy event objects
    :param catalog_names: names of catalogs (for file naming)
    :param comparison_magnitudes: list of two lists containing all magnitudes to match
            between the catalog and the reference catalog
    :return: saves matched events with matched magnitudes to a csv file
    """

    print('\nBuilding magnitude timeseries')
    for i in range(len(catalogs)):
        datalist = [[[] for j in range(len(comparison_magnitudes[i]))] for k in range(7)]
        for j in range(len(comparison_magnitudes[i])):
            for event in catalogs[i]:
                for magnitude in event.magnitudes:
                    if magnitude.magnitude_type == comparison_magnitudes[i][j]:
                        datalist[0][j].append(event.resource_id)
                        datalist[1][j].append(event.origins[0].time)
                        datalist[2][j].append(magnitude.magnitude_type)
                        datalist[3][j].append(magnitude.mag)
                        datalist[4][j].append(event.origins[0].latitude)
                        datalist[5][j].append(event.origins[0].longitude)
                        datalist[6][j].append(event.origins[0].depth)

            if len(datalist[0][j]) == 0:
                continue
            else:
                print('Saving data to file for catalog ' + str(catalog_names[i]) + ' for magnitude type ' +
                      str(comparison_magnitudes[i][j]) + '...')

                with open(catalog_names[i] + '_' + comparison_magnitudes[i][j] + '_timeseries.csv', 'w') as outfile:
                    outfile.write('eventID,origin_time,magnitude_type,magnitude,latitude,longitude,depth\n')
                with open(catalog_names[i] + '_' + comparison_magnitudes[i][j] + '_timeseries.csv', 'a') as outfile:
                    for n in range(len(datalist[0][j])):
                        outfile.write(
                            str(datalist[0][j][n]) + ',' + str(datalist[1][j][n]) + ',' + datalist[2][j][n]
                            + ',' + str(datalist[3][j][n]) + ',' + str(datalist[4][j][n]) + ',' + str(datalist[5][j][n])
                            + ',' + str(datalist[6][j][n]) + '\n')


def GeoNet_Mw(minmagnitude, starttime):

    """
    Query GeoNet Mw catalog from GitHub
    :param minmagnitude:
    :param starttime:
    :return: saves to file datalist in the same format as generate_timeseries function
    """

    print('\nBuilding GeoNet Mw timeseries')

    URL = "https://raw.githubusercontent.com/GeoNet/data/master/moment-tensor/GeoNet_CMT_solutions.csv"
    result = curl(URL).decode('ascii')

    datalist = [[] for i in range(7)]
    rc = 0
    for row in result.split('\n'):
        if rc == 0:
            rc += 1
            continue
        else:
            rowsplit = row.split(',')
            if len(rowsplit) < 2:
                continue
            else:
                time = rowsplit[1]
                if ((datetime.datetime.strptime(time, '%Y%m%d%H%M%S') >= datetime.datetime.strptime(starttime, '%Y-%m-%dT%H:%M:%S'))
                        and (float(rowsplit[11]) >= minmagnitude)):


                    datalist[0].append('smi:nz.org.geonet/' + rowsplit[0])
                    datalist[1].append(time[0:4] + '-' + time[4:6] + '-' + time[6:8] + 'T' + time[8:10] + ':' +
                                       time[10:12] + ':' + time[12:14] + '.000000Z')
                    datalist[2].append('Mw')
                    datalist[3].append(rowsplit[11])
                    datalist[4].append(rowsplit[2])
                    datalist[5].append(rowsplit[3])
                    datalist[6].append(rowsplit[13] + '000')

    if len(datalist[0]) == 0:
        return
    else:
        with open('GeoNet_catalog_Mw_timeseries.csv', 'w') as outfile:
            outfile.write('eventID,origin_time,magnitude_type,magnitude,latitude,longitude,depth\n')
        with open('GeoNet_catalog_Mw_timeseries.csv', 'a') as outfile:
            for n in range(len(datalist[0])):
                outfile.write(
                    str(datalist[0][n]) + ',' + str(datalist[1][n]) + ',' + datalist[2][n]
                    + ',' + str(datalist[3][n]) + ',' + str(datalist[4][n]) + ',' + str(datalist[5][n])
                    + ',' + str(datalist[6][n]) + '\n')


def parse_data(filelist, split_str):

    """
    Parse magnitude timeseries or matches from csv file
    :param filelist: list of files
    :param split_str: string to split filename by to get file type
    :return: list of data in file, list of type of data in files
    """

    data_types = []
    datalist = [[[] for j in range(len(filelist))] for i in range(7)]
    for j in range(len(filelist)):
        file = filelist[j]
        with open(file, 'r') as infile:
            rc = 0
            for row in infile:
                if rc == 0:
                    rc += 1
                    continue
                else:
                    rowsplit = row.split(',')
                    for i in range(len(rowsplit)):
                        datalist[i][j].append(str(rowsplit[i]))
            if len(datalist[0][j]) > 0:
                data_types.append(file.split('/')[-1].split(split_str)[0])

    return datalist, data_types

def cumulative_sum(times):

    """
    Generate a cumulative sum timeseries for event times in times
    :param times: list of ISO8601 format times (strings)
    :return: list of event times (datetime objects), cumulative sum at each event time (list of int)
    """

    cumulative_event_sum = 0
    cumulative_event_sums = []
    event_times = []

    # Generate cumulative sum list
    # and list of event times

    times.sort()
    for time in times:

        cumulative_event_sum += 1
        cumulative_event_sums.append(cumulative_event_sum)

        event_time = datetime.datetime.strptime(time, '%Y-%m-%dT%H:%M:%S.%fZ')
        event_times.append(event_time)

    return event_times, cumulative_event_sums


def plot_timeseries(magnitude_timeseries, timeseries_types):

    """
    Plot cumulative sum timeseries for magnitude types in magnitude_timeseries
    :param magnitude_timeseries:
    :param timeseries_types:
    :return: a single plot containing cumulative sum timeseries of all magnitude types
    """

    for n in range(len(magnitude_timeseries[0])):
        event_times, cumulative_event_sums = cumulative_sum(magnitude_timeseries[1][n])
        plt.scatter(event_times, cumulative_event_sums, label=timeseries_types[n])

    plt.ylabel('cumulative number of events')
    plt.xlabel('time')
    plt.legend()
    plt.show()


# Set script parameters

minmagnitude = 5     # minimum event magnitude to get from catalog
minlatitude, maxlatitude = -60, -13  # minimum and maximum latitude for event search window
minlongitude, maxlongitude = 145, -145  # minimum and maximum longitude for event search window
starttime = '2012-01-01T00:00:00' #event query starttime

max_dt = 100  # maximum time (s) between events in separate catalogs for them to be considered records of the same earthquake
max_dist = 1000  # maximum distance (km) "

# Define catalogs and their associated FDSN webservice event URL

catalog_names = ['GeoNet_catalog', 'USGS_catalog']
services = ["https://service.geonet.org.nz/fdsnws/event/1/", "https://earthquake.usgs.gov/fdsnws/event/1/"]
catalogs = [[] for i in range(len(catalog_names))]

# Define comparison magnitudes: first nested list is from GeoNet catalog, second if from USGS
# Code will do all combinations across the two catalogs

comparison_magnitudes = [['MLv', 'mB', 'Mw(mB)', 'M', 'Mw'], ['mww']]

# # Build event catalogs from FDSN
#
# print('\nSearching earthquake catalogs for events above magnitude ' + str(minmagnitude) +
#       ' between ' + str(minlatitude) + ' and ' + str(maxlatitude) + ' degrees latitude and ' +
#       str(minlongitude) + ' and ' + str(maxlongitude) + ' degrees longitude after ' + starttime)
#
# for n in range(len(catalogs)):
#
#     catalogs[n] = FDSN_event_query(services[n], minmagnitude, minlongitude, maxlongitude,
#                                   minlatitude, maxlatitude, starttime)
#     print('\n' + str(len(catalogs[n])) + ' events were found in catalog ' + str(n + 1))
#
# # Create a timeseries of use of each magnitude type in comparison_magnitudes for each catalog
#
# build_magnitude_timeseries(catalogs, catalog_names, comparison_magnitudes)
#
# # Create GeoNet Mw timeseries from GitHub, if desired
#
# try:
#     if 'Mw' in comparison_magnitudes[catalog_names.index('GeoNet_catalog')]:
#         GeoNet_Mw(minmagnitude, starttime)
# except:
#     pass

# Load timeseries data from files and do event matching

magnitude_timeseries_files = glob.glob('/home/samto/*timeseries.csv')
magnitude_timeseries, timeseries_types = parse_data(magnitude_timeseries_files, '_timeseries')

# match_magnitudes(magnitude_timeseries, timeseries_types, comparison_magnitudes, max_dt, max_dist)

magnitude_match_files = glob.glob('/home/samto/*matches.csv')
magnitude_matches, match_types = parse_data(magnitude_match_files, '_magnitude')

# # Sort data by alphabetic magnitude type
#
# for n in range(len(magnitude_timeseries)):
#     _, magnitude_timeseries[n] = zip(*sorted(zip(timeseries_types, magnitude_timeseries[n])))
# timeseries_types.sort()
#
# # Generate cumulative sum timeseries for each magnitude and plot it
#
# plot_timeseries(magnitude_timeseries, timeseries_types)

# Create matches between catalogs: get all possible magnitudes

datalist = [[] for i in range(5 + len(timeseries_types))]
columns = ['eventID', 'length', 'latitude', 'longitude', 'depth']
for magnitude_type in timeseries_types:
    columns.append(magnitude_type.split('_')[-1])

# Generate list of events with magnitude matches

event_list = []
for n in range(len(magnitude_matches[0])):
    for m in range(len(magnitude_matches[0][n])):
        event_list.append(magnitude_matches[0][n][m])
event_list = list(set(event_list))

# Populate datalist

datalist = [[[] for m in range(len(event_list))] for n in range(len(columns))]
for n in range(len(magnitude_matches)):
    for m in range(len(magnitude_matches[n])):
        column_idices = [0, columns.index(match_types[m].split('_')[0]), columns.index(match_types[m].split('_')[1]), 1, 2, 3, 4]
        for k in range(len(magnitude_matches[n][m])):
            event_index = event_list.index(magnitude_matches[0][m][k])
            if n == len(magnitude_matches) - 1:
                datalist[column_idices[n]][event_index] = str(float(magnitude_matches[n][m][k]))
            else:
                datalist[column_idices[n]][event_index] = magnitude_matches[n][m][k]

# Write datalist to file

with open('magnitude_matches_all.csv', 'w') as outfile:
    header = ""
    for column in columns:
        header += column + ','
    header = header[:-1]
    outfile.write(header + '\n')

with open('magnitude_matches_all.csv', 'a') as outfile:
    for m in range(len(datalist[0])):
        outstr = ""
        for n in range(len(datalist)):
            try:
                outstr += datalist[n][m] + ','
            except:
                outstr += "nan,"
        outfile.write(outstr[:-1] + '\n')

