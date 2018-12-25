from obspy.io.quakeml.core import Unpickler
import pycurl
from io import BytesIO
import math
import matplotlib.pyplot as plt
import datetime
import pandas as pd

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

                print('\nAttempting catalog query for events between M ' + str(magnitude_limits[i - 1]) + ' and ' + str(magnitude_limits[i]))
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


def generate_lengths(event, reference_catalog, max_dt, max_dist):
    """
    :param event: obspy event object
    :param reference_catalog: list of obspy event objects to calculate lengths to
    :param max_dt: maximum seconds (absolute) between events for them to be matched
    :param max_dist: maximum distance (absolute) between events for them to be matched
    :return: list containing combined temporal and spatial length between the event and
                all reference events
    """

    # Initialise list to save lengths in

    lengths = [[], []]

    # Get event position in cartesian coordinates

    ETi, ELa, ELo, EDe = [event.origins[0].time, event.origins[0].latitude,
                          event.origins[0].longitude, event.origins[0].depth]
    Ex, Ey, Ez = to_cartesian(ELa, ELo, EDe)

    # Compare event temporal and spatial position to those of reference events

    for n in range(len(reference_catalog)):

        reference_event = reference_catalog[n]

        # Get reference event position in cartesian coordinates

        RETi, RELa, RELo, REDe = [reference_event.origins[0].time, reference_event.origins[0].latitude,
                                  reference_event.origins[0].longitude, reference_event.origins[0].depth]

        REx, REy, REz = to_cartesian(RELa, RELo, REDe)

        # Calculate temporal length and check temporal matching limit

        temporal_length = abs(ETi - RETi)
        if temporal_length > max_dt:
            continue

        # Calculate spatial length and check spatial matching limit

        spatial_length = math.sqrt((Ex - REx) ** 2 + (Ey - REy) ** 2 + (Ez - REz) ** 2) / 1000.0
        if spatial_length > max_dist:
            continue

        # Store values

        lengths[0].append(math.sqrt(temporal_length ** 2 + spatial_length ** 2))
        lengths[1].append(n)

    return lengths


def to_cartesian(latitude, longitude, depth):
    """
    Convert a point on the Earth in spherical coordinates to cartesian coordinates
    :param latitude: point latitude in decimal degrees, south is negative
    :param longitude: point longitude in decimal degrees, west is negative
    :param depth: point depth in metres, down is positive
    :return: cartesian coordinates of the point
    """

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


def match_magnitudes(event_catalog, reference_catalog, comparison_magnitudes, max_dt, max_dist):

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

    print('\nFinding event matches between the two catalogs')

    magnitude_types = [[] for i in range(len(comparison_magnitudes[0]) * len(comparison_magnitudes[1]))]
    datalist = [[[] for i in range(len(comparison_magnitudes[0]) * len(comparison_magnitudes[1]))] for j in range(4)]

    process_start = datetime.datetime.now()
    for n in range(len(event_catalog)):
        event = event_catalog[n]
        lengths = generate_lengths(event, reference_catalog, max_dt, max_dist)
        print(str(len(lengths[0])) + ' reference events found for the ' + str(n) + 'th event in the catalog')
        if len(lengths[0]) > 0:
            matched_event = reference_catalog[lengths[1][lengths[0].index(min(lengths[0]))]]
            cc = 0  # combination counter
            for dep_mag_type in comparison_magnitudes[0]:
                for indep_mag_type in comparison_magnitudes[1]:
                    magnitude_types[cc] = [dep_mag_type, indep_mag_type]
                    print('Finding event matches for magnitude types ' + dep_mag_type + ' and ' + indep_mag_type)
                    for event_magnitude in event.magnitudes:
                        if event_magnitude.magnitude_type == dep_mag_type:
                            for matched_event_magnitude in matched_event.magnitudes:
                                if matched_event_magnitude.magnitude_type == indep_mag_type:
                                    datalist[0][cc].append(event.resource_id)
                                    datalist[1][cc].append(event_magnitude.mag)
                                    datalist[2][cc].append(matched_event_magnitude.mag)
                                    datalist[3][cc].append(lengths[0][lengths[0].index(min(lengths[0]))])
                    print('    ' + str(
                        len(datalist[0][cc])) + ' matched events with desired magnitude types were found')
                    cc += 1

    process_end = datetime.datetime.now()
    print('That took ' + str((process_end - process_start).total_seconds()) + ' seconds')

    # Save data

    print('\nSaving data to file...')
    process_start = datetime.datetime.now()
    for i in range(cc):
        with open(magnitude_types[i][0] + '_' + magnitude_types[i][1] + '_' + 'magnitude_matches.csv', 'w') as outfile:
            outfile.write('eventID,' + magnitude_types[i][0] + ',' + magnitude_types[i][1] + ',' + 'length' + '\n')
        with open(magnitude_types[i][0] + '_' + magnitude_types[i][1] + '_' + 'magnitude_matches.csv', 'a') as outfile:
            for n in range(len(datalist[0][i])):
                outfile.write(
                    str(datalist[0][i][n]) + ',' + str(datalist[1][i][n]) + ',' + str(datalist[2][i][n])
                    + ',' + str(datalist[3][i][n]) + '\n')
    process_end = datetime.datetime.now()
    print('That took ' + str((process_end - process_start).total_seconds()) + ' seconds')


def magnitude_timeseries(catalogs, catalog_names, comparison_magnitudes):

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
        datalist = [[[] for j in range(len(comparison_magnitudes[i]))] for k in range(4)]
        for j in range(len(comparison_magnitudes[i])):
            for event in catalogs[i]:
                for magnitude in event.magnitudes:
                    if magnitude.magnitude_type == comparison_magnitudes[i][j]:
                        datalist[0][j].append(event.resource_id)
                        datalist[1][j].append(event.origins[0].time)
                        datalist[2][j].append(magnitude.magnitude_type)
                        datalist[3][j].append(magnitude.mag)

            print('Saving data to file for catalog ' + str(catalog_names[i]) + ' for magnitude type ' +
                  str(comparison_magnitudes[i][j]) + '...')

            with open(catalog_names[i] + '_' + comparison_magnitudes[i][j] + '_timeseries.csv', 'w') as outfile:
                outfile.write('eventID,' + 'origin_time' + ',' + 'magnitude_type' + ',' + 'magnitude' + '\n')
            with open(catalog_names[i] + '_' + comparison_magnitudes[i][j] + '_timeseries.csv', 'a') as outfile:
                for n in range(len(datalist[0][j])):
                    outfile.write(
                        str(datalist[0][j][n]) + ',' + str(datalist[1][j][n]) + ',' + datalist[2][j][n]
                        + ',' + str(datalist[3][j][n]) + '\n')


# Set script parameters

minmagnitude = 3     # minimum event magnitude to get from catalog
minlatitude, maxlatitude = -90, 90  # minimum and maximum latitude for event search window
minlongitude, maxlongitude = 0, -180  # minimum and maximum longitude for event search window

max_dt = 100  # maximum time (s) between events in separate catalogs for them to be considered records of the same earthquake
max_dist = 1000  # maximum distance (km) "

# Define catalogs and their associated FDSN webservice event URL

catalog_names = ['GeoNet_catalog', 'USGS_catalog']
services = ["https://service.geonet.org.nz/fdsnws/event/1/", "https://earthquake.usgs.gov/fdsnws/event/1/"]
catalogs = [[] for i in range(len(catalog_names))]

# Define comparison magnitudes: first nested list is from GeoNet catalog, second if from USGS
# Code will do all combinations across the two catalogs

comparison_magnitudes = [['MLv', 'mB', 'Mw(mB)', 'Mw', 'M', 'ML'], ['mw', 'mwc', 'mwb', 'mww', 'mB']]

# Build event catalogs

print('Searching earthquake catalogs for events above magnitude ' + str(minmagnitude) +
      ' between ' + str(minlatitude) + ' and ' + str(maxlatitude) + ' degrees latitude and ' +
      str(minlongitude) + ' and ' + str(maxlongitude) + ' degrees longitude')

for n in range(len(catalogs)):

    catalogs[n] = FDSN_event_query(services[n], minmagnitude, minlongitude, maxlongitude,
                                  minlatitude, maxlatitude)
    print('\n' + str(len(catalogs[n])) + ' events were found in catalog ' + str(n + 1))

# Create a timeseries of use of each magnitude type in comparison_magnitudes for each catalog

magnitude_timeseries(catalogs, catalog_names, comparison_magnitudes)

# Match events between catalogs and extract magnitudes then save the data

match_magnitudes(catalogs[0], catalogs[1], comparison_magnitudes, max_dt, max_dist)