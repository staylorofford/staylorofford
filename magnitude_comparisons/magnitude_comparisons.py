from obspy.io.quakeml.core import Unpickler
import pycurl
from io import BytesIO
import math
import matplotlib.pyplot as plt
import datetime
import glob
from scipy.odr import Model, Data, ODR

quakeml_reader = Unpickler()

def FDSN_event_query(service, minmagnitude, minlongitude, maxlongitude,
                     minlatitude, maxlatitude, starttime = '0000-01-01T00:00:00',
                     endtime = '9999-01-01T00:00:00', maxmagnitude = 10):
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
    :param endtime:
    :return:
    """

    # Adjust format if required

    if "usgs" in service:
        maxlongitude += 360

    # Curl FDSN response and parse quakeml bytes string through obspy
    # using an interative query approach when a single query fails

    factor = 1
    success = False
    while not success:

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
                                "&starttime=", starttime, "&endtime=", endtime))

            try:

                print('\nAttempting FDSN catalog query for events between M ' + str(magnitude_limits[i - 1]) + ' and '
                      + str(magnitude_limits[i]))
                queryresult = curl(query)

                catalog = quakeml_reader.loads(queryresult)

                events.extend(catalog.events)
                print('Catalog now has ' + str(len(events)) + ' events')
                successes += 1

            except:

                print('Failed!')
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

    buffer = BytesIO()
    c = pycurl.Curl()
    c.setopt(c.URL, curlstr)
    c.setopt(c.WRITEDATA, buffer)
    c.perform()
    c.close()

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


def GeoNet_Mw(minmagnitude, starttime, endtime):

    """
    Query GeoNet Mw catalog from GitHub
    :param minmagnitude:
    :param starttime:
    :param endtime:
    :return: saves to file datalist in the same format as generate_timeseries function
    """

    print('\nBuilding GeoNet Mw timeseries')

    URL = "https://raw.githubusercontent.com/GeoNet/data/master/moment-tensor/GeoNet_CMT_solutions.csv"
    result = curl(URL).decode('ascii')
    print("")

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
                if ((datetime.datetime.strptime(starttime, '%Y-%m-%dT%H:%M:%S')) <=
                        datetime.datetime.strptime(time, '%Y%m%d%H%M%S') <=
                        datetime.datetime.strptime(endtime, '%Y-%m-%dT%H:%M:%S')) and \
                        float(rowsplit[11]) >= minmagnitude:

                    try:
                        URL = "https://service.geonet.org.nz/fdsnws/event/1/query?eventid=" + rowsplit[0]
                        event = quakeml_reader.loads(curl(URL))[0]
                    except:
                        print('Event with eventID ' + rowsplit[0] + ' not in GeoNet catalog?')

                    datalist[0].append('smi:nz.org.geonet/' + rowsplit[0])
                    datalist[2].append('Mw')
                    datalist[3].append(rowsplit[11])

                    # Get timing and location details from equivalent GeoNet catalog event

                    datalist[1].append(event.origins[0].time)
                    datalist[4].append(event.origins[0].latitude)
                    datalist[5].append(event.origins[0].longitude)
                    datalist[6].append(event.origins[0].depth)

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


def parse_data(filelist, split_str, starttime, endtime):

    """
    Parse magnitude timeseries or matches from csv file, using start and end times for time filtering
    :param filelist: list of files
    :param split_str: string to split filename by to get file type
    :param starttime: do not load data from before this time
    :param endtime: do not load data from after this time
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
                else:
                    rowsplit = row.split(',')
                    if (starttime <=
                        datetime.datetime.strptime(rowsplit[1], '%Y-%m-%dT%H:%M:%S.%fZ') <=
                        endtime):
                        for i in range(len(rowsplit)):
                            datalist[i][j].append(str(rowsplit[i]))
                if len(datalist[0][j]) > 0:
                    data_types.append(file.split('/')[-1].split(split_str)[0])

    return datalist, data_types


def match_magnitudes(magnitude_timeseries, timeseries_types, catalog_names, comparison_magnitudes, max_dt, max_dist,
                     show_matching=False):

    """
    Match events between two catalogs and save all events and their matched magnitudes to file
    :param magnitude_timeseries:
    :param timeseries_types:
    :param catalog_names:
    :param comparison_magnitudes: list of two lists containing all magnitudes to match
                                between the catalog and the reference catalog
    :param max_dt: maximum seconds (absolute) between events for them to be matched
    :param max_dist: maximum distance (absolute) between events for them to be matched
    :param show_matching: whether to plot relative distance and time of all events that have matched
    :return: saves matched events with matched magnitudes to a csv file
    """

    # Build column types for output csv (columns)

    columns = ['eventID', 'length', 'latitude', 'longitude', 'depth']
    magnitudes_columns = []
    for magnitude_type in timeseries_types:
        magnitudes_columns.append(magnitude_type.split('_')[-1])
    magnitudes_columns = list(set(magnitudes_columns))
    magnitudes_columns.sort()
    columns.extend(magnitudes_columns)

    # Generate list of events for output csv (rows)

    event_list = []
    for n in range(len(magnitude_timeseries[0])):
        for m in range(len(magnitude_timeseries[0][n])):
            if timeseries_types[n].split('_')[0] in catalog_names[0]: # Only populate the event list with the non-reference catalog
                event_list.append(magnitude_timeseries[0][n][m])
    event_list = list(set(event_list))

    # Pre-populated eventID, location, and length in datalist prior to matching

    datalist = [[[] for m in range(len(event_list))] for n in range(len(columns))]

    for n in range(len(magnitude_timeseries[0])):
        for k in range(len(magnitude_timeseries[0][n])):
            try:
                event_index = event_list.index(magnitude_timeseries[0][n][k])
                datalist[0][event_index] = magnitude_timeseries[0][n][k]
                datalist[1][event_index] = '0'  # Length 0 for internal matches: external matches will overwrite
                datalist[2][event_index] = magnitude_timeseries[4][n][k]
                datalist[3][event_index] = magnitude_timeseries[5][n][k]
                datalist[4][event_index] = str(float(magnitude_timeseries[6][n][k]))  # Remove trailing newline
            except: # Fails when the event is not from the non-reference catalog
                pass

    # Match events between timeseries and fill in magnitude information in the datalist
    complete_pairs = []
    matched_temporal_lengths = []
    matched_spatial_lengths = []
    for n in range(len(timeseries_types)):
        if timeseries_types[n].split('_')[0] == catalog_names[0].split('_')[0] and \
                timeseries_types[n].split('_')[2] in comparison_magnitudes[0]:
                # We have one of our first sets of comparison magnitudes
            for m in range(len(timeseries_types)):

                if str(m) + ',' + str(n) in complete_pairs:
                    continue

                print('Looking for matching events with magnitude types ' + timeseries_types[n] + ' and ' + timeseries_types[m] + '...')

                if timeseries_types[m].split('_')[0] == catalog_names[0].split('_')[0] and \
                        timeseries_types[m].split('_')[2] in comparison_magnitudes[0]:
                        # We have another of our first sets of comparison magnitudes
                        # Find matches and load data into datalist
                        for k in range(len(magnitude_timeseries[0][n])): # Go through all the entires for the nth magnitude type
                            event_index = event_list.index(magnitude_timeseries[0][n][k])
                            for l in range(len(magnitude_timeseries[0][m])):  # Go through all the entires for the mth magnitude type
                                if magnitude_timeseries[0][n][k] == magnitude_timeseries[0][m][l]: # Match based on eventID
                                    datalist[columns.index(timeseries_types[n].split('_')[2])][event_index] = magnitude_timeseries[3][n][k]
                                    datalist[columns.index(timeseries_types[m].split('_')[2])][event_index] = magnitude_timeseries[3][m][l]


                elif timeseries_types[m].split('_')[0] == catalog_names[1].split('_')[0] and \
                        timeseries_types[m].split('_')[2] in comparison_magnitudes[1]:
                        # We have one of our second sets of comparison magnitudes
                    for k in range(len(magnitude_timeseries[0][n])):
                        event_index = event_list.index(magnitude_timeseries[0][n][k])

                        # Calculate 2D length between event and reference events for matching criteria

                        temporal_lengths = []
                        spatial_lengths = []
                        lengths = []
                        indices = []
                        ETi, ELa, ELo, EDe = [datetime.datetime.strptime(magnitude_timeseries[1][n][k], '%Y-%m-%dT%H:%M:%S.%fZ'),
                                              float(magnitude_timeseries[4][n][k]), float(magnitude_timeseries[5][n][k]),
                                              float(magnitude_timeseries[6][n][k])]
                        Ex, Ey, Ez = to_cartesian(ELa, ELo, EDe)

                        for l in range(len(magnitude_timeseries[0][m])):

                            RETi, RELa, RELo, REDe = [datetime.datetime.strptime(magnitude_timeseries[1][m][l], '%Y-%m-%dT%H:%M:%S.%fZ'),
                                              float(magnitude_timeseries[4][m][l]), float(magnitude_timeseries[5][m][l]),
                                              float(magnitude_timeseries[6][m][l])]
                            REx, REy, REz = to_cartesian(RELa, RELo, REDe)

                            temporal_length = abs((ETi - RETi).total_seconds())
                            if temporal_length > max_dt:
                                continue
                            else:
                                temporal_lengths.append(temporal_length)

                            spatial_length = math.sqrt((Ex - REx) ** 2 + (Ey - REy) ** 2 + (Ez - REz) ** 2) / 1000.0
                            if spatial_length > max_dist:
                                continue
                            else:
                                spatial_lengths.append(spatial_length)

                            lengths.append(math.sqrt(temporal_length ** 2 + spatial_length ** 2))
                            indices.append(l)

                        if len(lengths) > 0:
                            match_idx = indices[lengths.index(min(lengths))] # Event match is that with smallest length

                            matched_spatial_lengths.append(spatial_lengths[indices.index(match_idx)])
                            matched_temporal_lengths.append(temporal_lengths[indices.index(match_idx)])

                            datalist[1][event_index] = str(lengths[lengths.index(min(lengths))])
                            datalist[columns.index(timeseries_types[n].split('_')[2])][event_index] = \
                            magnitude_timeseries[3][n][k]
                            datalist[columns.index(timeseries_types[m].split('_')[2])][event_index] = \
                            magnitude_timeseries[3][m][match_idx]

                complete_pairs.append(str(n) + ',' + str(m))

    if show_matching:

        print('\nNOTE: To investigate the spread of matched data in an unconstrained format, ensure maximum limits are'
              '>=1E9\n')

        plt.scatter(matched_temporal_lengths, matched_spatial_lengths, s=2)
        plt.xlabel('relative time (s)', labelpad=15)
        plt.ylabel('relative distance (km)', labelpad=15)
        plt.title('relative distance vs. time for all matched events')
        plt.tight_layout()
        plt.show()

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


def f(P, x):

    """
    Point on line function to use in orthogonal regression
    :param P: m, c values for line
    :param x: dataset values
    :return point on line values for x value in x
    """

    return P[0] * x + P[1]


def orthregress(x, y):

    """
    Perform an orthogonal distance regression,
    :param x: dataset 1 values
    :param y: dataset 2 values
    :return: gradient, intercept of orthogonal distance regression
    """

    # Run orthogonal regression

    model = Model(f)
    data = Data(x, y)
    odr = ODR(data, model, beta0=[1, 1])
    output = odr.run()
    m, c = output.beta

    return m, c


def probability(sample_magnitudes, sample_depths, sample_times,
                min_mag, max_mag, min_depth, max_depth,
                probability_time_period, number_of_events):

    """
    Calculate the Poisson probability that a number of earthquakes will occur with magnitude value
    within a given range within a given time period
    :param sample_magnitudes: list of magnitude values for earthquakes in the sample
    :param sample_depths: list of event depths for earthquakes in the sample
    :param sample_times: list of time values for earthquakes in the sample
    :param min_mag: minimum magnitude to include in the probability
    :param max_mag: maximum magnitude to include in the probability
    :param min_depth: minimum event depth to include in the probability (in km)
    :param max_depth: maximum event depth to include in the probability (in km)
    :param probability_time_period: length of time (in hours) to calculate the probability of at least 1
                                    earthquake occurring over
    :param number_of_events: number of events to calculate probability of occurrence for, use 'any' to calculate
                             the probability at least 1 event
    :return: probability of earthquake(s) happening in the length of time input with a magnitude in the range input
    """

    # Calculate set of magnitudes to calculate mean rate for

    set_magnitudes = []
    set_times = []
    for m in range(len(sample_magnitudes)):
        if min_mag <= sample_magnitudes[m] <= max_mag and min_depth <= sample_depths[m] <= max_depth:
            set_magnitudes.append(sample_magnitudes[m])
            set_times.append(sample_times[m])

    # Calculate rate of an earthquakes in the magnitude range per hour

    N = len(set_magnitudes)

    if len(set_times) < 2:
        print('Not enough data exists within the specified magnitude and depth limits for this time period!')
        return 'nan', 'nan'

    mean_rate = N / (max(set_times) - min(set_times)).total_seconds() * 3600.0

    # Calculate rate of an earthquake in the magnitude range for the given time period

    mean_rate_time_period = mean_rate * probability_time_period

    # Calculate Poisson probability of 1 or more earthquakes happening in this time period

    if number_of_events == 'any':
        p = 1 - (math.exp(-mean_rate_time_period) * (mean_rate_time_period ** 0) / math.factorial(0))
    else:
        p = (math.exp(-1 * mean_rate_time_period) * (mean_rate_time_period ** number_of_events) /
            math.factorial(number_of_events))

    return p, N


# Set data gathering parameters

minmagnitude = 5  # minimum event magnitude to get from catalog
minlatitude, maxlatitude = -67, -33  # minimum and maximum latitude for event search window
minlongitude, maxlongitude = 145, -175  # minimum and maximum longitude for event search window
starttime = '2001-05-18T00:00:00'  # event query starttime
endtime = 'now'  # event query endtime, 'now' will set it to the current time

if endtime == 'now':
    endtime = str(datetime.datetime.now().isoformat())[:19]

# Define catalogs and their associated FDSN webservice event URL

catalog_names = ['GeoNet_catalog', 'USGS_catalog']
services = ["https://service.geonet.org.nz/fdsnws/event/1/", "https://earthquake.usgs.gov/fdsnws/event/1/"]
catalogs = [[] for i in range(len(catalog_names))]

# Define comparison magnitudes: first nested list is from GeoNet catalog, second if from USGS
# Code will do all combinations across the two catalogs

comparison_magnitudes = [['MLv', 'ML', 'mB', 'Mw(mB)', 'M', 'Mw'], ['mww']]

# Set matching parameters

max_dt = 10  # maximum time (s) between events in separate catalogs for them to be considered records of the same earthquake
max_dist = 10000 # maximum distance (km) "

# Set probability parameters

probability_magnitude_types = comparison_magnitudes[0]  # magnitude types to find the largest magnitude for each event
min_mag = 6
max_mag = 10
min_depth = 0  # minimum depth of earthquake to include, in km
max_depth = 100  # maximum depth of earthquake to include, in km
probability_time_period = 24 * 7  # time period to calculate probability over, in hours
number_of_events = 'any'  # number of events to calculate probability for
time_bins = 19  # number of bins to split the query time period into
bin_overlap = 0.9  # percentage each time bin should overlap
# ^ need something to plot against to know if probabilities are reasonable, perhaps total # of events in time period,
# or # of stations? looks like pre-2010 data behaves very differently, suggesting something network-related influences
# this, and this data should not be included in a long-term dataset trying to determine likelihood of earthquakes!

# Set what level of processing you want the script to do

build_FDSN_timeseries = False
build_GeoNet_Mw_timeseries = False
probabilities = True
matching = False
show_matching = False

# Build event catalogs from FDSN

if build_FDSN_timeseries == True:

    print('\nSearching earthquake catalogs for events above magnitude ' + str(minmagnitude) +
          ' between ' + str(minlatitude) + ' and ' + str(maxlatitude) + ' degrees latitude and ' +
          str(minlongitude) + ' and ' + str(maxlongitude) + ' degrees longitude after ' + str(starttime))

    for n in range(len(catalogs)):

        catalogs[n] = FDSN_event_query(services[n], minmagnitude, minlongitude, maxlongitude,
                                       minlatitude, maxlatitude, starttime)
        print('\n' + str(len(catalogs[n])) + ' events were found in catalog ' + str(n + 1))

    # Create a timeseries of use of each magnitude type in comparison_magnitudes for each catalog

    build_magnitude_timeseries(catalogs, catalog_names, comparison_magnitudes)

# Build GeoNet Mw catalog

if build_GeoNet_Mw_timeseries == True:

    GeoNet_Mw(minmagnitude, starttime, endtime)

if probabilities:

    # Convert date strings to datetime objects

    starttime = datetime.datetime.strptime(starttime, '%Y-%m-%dT%H:%M:%S')
    endtime = datetime.datetime.strptime(endtime, '%Y-%m-%dT%H:%M:%S')

    # Load timeseries data from files

    magnitude_timeseries_files = glob.glob('./*timeseries.csv')
    bin_length = (endtime - starttime).total_seconds() / time_bins
    current_start_time = starttime
    current_end_time = starttime + datetime.timedelta(seconds=bin_length)

    start_times = []
    end_times = []
    probabilities = []
    sample_numbers = []
    while current_end_time <= endtime:

        magnitude_timeseries, timeseries_types = parse_data(magnitude_timeseries_files, '_timeseries',
                                                            current_start_time, current_end_time)

        print('\nCalculating probabilities...\n')

        # Calculate probability for each magnitude type

        for m in range(len(magnitude_timeseries[2])):

            if len(magnitude_timeseries[2][m]) == 0:
                continue

            sample_times = []
            sample_magnitudes = []
            sample_depths = []

            if magnitude_timeseries[2][m][0] in probability_magnitude_types:
                for n in range(len(magnitude_timeseries[2][m])):
                    sample_times.append(datetime.datetime.strptime(magnitude_timeseries[1][m][n],
                                                                   '%Y-%m-%dT%H:%M:%S.%fZ'))
                    sample_magnitudes.append(float(magnitude_timeseries[3][m][n]))
                    sample_depths.append(float(magnitude_timeseries[6][m][n]) / 1000.0)
            else:
                continue

            p, N = probability(sample_magnitudes, sample_depths, sample_times,
                               min_mag, max_mag, min_depth, max_depth,
                               probability_time_period, number_of_events)

            if p != 'nan':

                print('Probability of ' + str(number_of_events) + ' events of magnitude type ' + str(timeseries_types[m]) +
                      ' between magnitude values ' + str(min_mag) + ' - ' + str(max_mag) +
                      ' betweeen depths of ' + str(min_depth) + ' - ' + str(max_depth) +
                      ' km occurring in ' + str(probability_time_period) + ' hours is ' +
                      str(p) + ', from ' + str(N) + ' samples over the period ' +
                      str(current_start_time.isoformat()) + ' - ' + str(current_end_time.isoformat()) + '\n')

    # Calculate probability using maximum magnitude for each event

        all_events = []
        all_times = []
        all_magnitudes = []
        all_depths = []
        for m in range(len(magnitude_timeseries[2])):
            if len(magnitude_timeseries[2][m]) == 0:  # Do not process empty bins
                continue
            if magnitude_timeseries[2][m][0] in probability_magnitude_types:
                for n in range(len(magnitude_timeseries[2][m])):
                    all_events.append(magnitude_timeseries[0][m][n])
                    all_times.append(datetime.datetime.strptime(magnitude_timeseries[1][m][n], '%Y-%m-%dT%H:%M:%S.%fZ'))
                    all_magnitudes.append(float(magnitude_timeseries[3][m][n]))
                    all_depths.append(float(magnitude_timeseries[6][m][n]) / 1000.0)

        # Find the largest magnitude for each event

        sample_events = []
        sample_times = []
        sample_magnitudes = []
        sample_depths = []
        for m in range(len(all_events)):
            if all_events[m] not in sample_events:
                all_sample_magnitudes = []
                for n in range(len(all_events)):
                    if all_events[m] == all_events[n]:
                        all_sample_magnitudes.append(all_magnitudes[n])
                sample_events.append(all_events[m])
                sample_times.append(all_times[m])
                sample_depths.append(all_depths[m])
                sample_magnitudes.append(max(all_sample_magnitudes))

        # Calculate probability

        p, N = probability(sample_magnitudes, sample_depths, sample_times,
                           min_mag, max_mag, min_depth, max_depth,
                           probability_time_period, number_of_events)

        if p != 'nan':

            print('Probability of ' + str(number_of_events) + ' events of any of the above magnitude types' +
                  ' between magnitude values ' + str(min_mag) + ' - ' + str(max_mag) +
                  ' betweeen depths of ' + str(min_depth) + ' - ' + str(max_depth) +
                  ' km occurring in ' + str(probability_time_period) + ' hours is ' +
                  str(p) + ', from ' + str(N) + ' over the period ' +
                  str(current_start_time.isoformat()) + ' - ' + str(current_end_time.isoformat()) + '\n')

            probabilities.append(p)
            sample_numbers.append(N)
            start_times.append(current_start_time)
            end_times.append(current_end_time)

        current_start_time += datetime.timedelta(seconds=int(round(bin_length * (1 - bin_overlap))))
        current_end_time += datetime.timedelta(seconds=int(round(bin_length * (1 - bin_overlap))))

    plot_times = []
    for n in range(len(probabilities)):
        plot_time = start_times[n] # + datetime.timedelta(seconds=(end_times[n] - start_times[n]).total_seconds() / 2)
        plt.scatter(plot_time, probabilities[n], color='k', s=2)
        plt.text(plot_time, probabilities[n], str(sample_numbers[n]))
        plot_times.append(plot_time)
    plt.plot(plot_times, probabilities, linestyle='--', marker='o', color='k')

    plt.xlabel('window mid-point', labelpad=15)
    plt.ylabel('maximum probability of earthquake', labelpad=15)
    plt.tight_layout()
    plt.show()

    if matching:

        if probabilities:

            # Return start and end times to ISO8601 format

            starttime = str(starttime.isoformat())
            endtime = str(endtime.isoformat())

        # Parse all data

        magnitude_timeseries, timeseries_types = parse_data(magnitude_timeseries_files, '_timeseries',
                                                            starttime, endtime)

        # Do matching

        print("Matching events within temporal and spatial distance limits and with the desired magnitude types")

        match_magnitudes(magnitude_timeseries, timeseries_types, catalog_names, comparison_magnitudes, max_dt, max_dist,
                         show_matching)

        # Load magnitude match data

        with open('./magnitude_matches_all.csv', 'r') as openfile:
            rc = 0
            for row in openfile:
                if rc == 0:
                    data_types = row.split(',')
                    data_types[-1] = data_types[-1][:-1]
                    datalist = [[] for i in range(len(data_types))]
                    rc += 1
                else:
                    rowsplit = row.split(',')
                    for i in range(len(rowsplit)):
                        datalist[i].append(rowsplit[i])

        print('Saving plots...')

        # Magnitude value plotting

        complete_pairs = []
        for n in range(5, len(data_types)):
            for m in range(5, len(data_types)):

                # Don't repeat plotting

                if n == m or str(m) + ',' + str(n) in complete_pairs:
                    continue

                # Build plotting dataset

                plt.figure()
                x = []
                y = []
                for k in range(len(datalist[0])):
                    if datalist[n][k][:3] != 'nan' and datalist[m][k][:3] != 'nan':
                        x.append(float(datalist[n][k]))
                        y.append(float(datalist[m][k]))

                # Plot the data and a reference line for when magnitudes are equivalent at all values

                plt.scatter(x, y, s=2)
                plt.plot(range(11), range(11), color='k', linestyle='--', linewidth=0.5, alpha=0.5)

                # Do an orthogonal distance regression and plot this overtop of the data, showing the
                # linear relationship between the magnitudes

                slope, intercept = orthregress(x, y)
                plt.plot([0, 10], [intercept, slope * 10 + intercept], color='k', linewidth=0.5, alpha=0.5)

                # Add plot features for clarity

                plt.xlabel(data_types[n])
                plt.ylabel(data_types[m])
                plt.grid(which='major', axis='both', linestyle='-', alpha=0.5)
                plt.xlim(2, 10)
                plt.ylim(2, 10)
                plt.title('m=' + str(slope)[:5] + ', c=' + str(intercept)[:5], y=1.03)
                plt.gca().set_aspect('equal', adjustable='box')
                plt.savefig(data_types[n] + '_' + data_types[m] + '.png', format='png', dpi=300)
                plt.close()

                complete_pairs.append(str(n) + ',' + str(m))
