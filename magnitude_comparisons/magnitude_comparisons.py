from obspy.io.quakeml.core import Unpickler
import pycurl
from io import BytesIO
import math
import matplotlib.pyplot as plt
import datetime

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
    :param event: event from an obspy Catalog object
    :param reference_catalog: obspy Catalog object to calculate lengths to
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


# Set script parameters

minmagnitude = 4 # minimum event magnitude to get from catalog
minlatitude, maxlatitude = -60, -13  # minimum and maximum latitude for event search window
minlongitude, maxlongitude = 145, -145  # minimum and maximum longitude for event search window

max_dt = 100  # maximum time (s) between events in separate catalogs for them to be considered records of the same earthquake
max_dist = 1000  # maximum distance (km) "

# Define comparison magnitudes: first nested list is from GeoNet catalog, second if from USGS
# Code will do all combinations across the two catalogs

comparison_magnitudes = [['MLv', 'mB', 'Mw(mB)', 'M'], ['mw', 'mwc', 'mwb', 'mww']]

# Build event catalogs

print('Searching earthquake catalogs for events above magnitude ' + str(minmagnitude) +
      ' between ' + str(minlatitude) + ' and ' + str(maxlatitude) + ' degrees latitude and ' +
      str(minlongitude) + ' and ' + str(maxlongitude) + ' degrees longitude')

# Build GeoNet event catalog

service = "https://service.geonet.org.nz/fdsnws/event/1/"
GeoNet_catalog = FDSN_event_query(service, minmagnitude, minlongitude, maxlongitude,
                                  minlatitude, maxlatitude)
print('\n' + str(len(GeoNet_catalog)) + ' events were found in the GeoNet catalog')

# Build USGS event catalog

service = "https://earthquake.usgs.gov/fdsnws/event/1/"
USGS_catalog = FDSN_event_query(service, minmagnitude, minlongitude, maxlongitude,
                                minlatitude, maxlatitude)
print('\n' + str(len(USGS_catalog)) + ' events were found in the USGS catalog')

# Match events between catalogs and extract magnitudes then save the data

print('\nFinding event matches between the two catalogs')
magnitude_types = [[] for i in range(len(comparison_magnitudes[0]) * len(comparison_magnitudes[1]))]
eventIDs = [[] for i in range(len(comparison_magnitudes[0]) * len(comparison_magnitudes[1]))]
dependent_magnitudes = [[] for i in range(len(comparison_magnitudes[0]) * len(comparison_magnitudes[1]))]
independent_magnitudes = [[] for i in range(len(comparison_magnitudes[0]) * len(comparison_magnitudes[1]))]
length = [[] for i in range(len(comparison_magnitudes[0]) * len(comparison_magnitudes[1]))]
process_start = datetime.datetime.now()
for n in range(len(GeoNet_catalog)):
    event = GeoNet_catalog[n]
    lengths = generate_lengths(event, USGS_catalog, max_dt, max_dist)
    print(str(len(lengths[0])) + ' reference events found for the ' + str(n) + 'th event in the catalog')
    if len(lengths[0]) > 0:
        matched_event = USGS_catalog[lengths[1][lengths[0].index(min(lengths[0]))]]
        cc = 0 # combination counter
        for dep_mag_type in comparison_magnitudes[0]:
            for indep_mag_type in comparison_magnitudes[1]:
                magnitude_types[cc] = [dep_mag_type, indep_mag_type]
                print('Finding event matches for magnitude types ' + dep_mag_type + ' and ' + indep_mag_type)
                for event_magnitude in event.magnitudes:
                    if event_magnitude.magnitude_type == dep_mag_type:
                        for matched_event_magnitude in matched_event.magnitudes:
                            if matched_event_magnitude.magnitude_type == indep_mag_type:
                                eventIDs[cc].append(str(event.resource_id).split('/')[-1])
                                dependent_magnitudes[cc].append(event_magnitude.mag)
                                independent_magnitudes[cc].append(matched_event_magnitude.mag)
                                length[cc].append(lengths[0][lengths[0].index(min(lengths[0]))])
                print('    ' + str(len(dependent_magnitudes[cc])) + ' matched events with desired magnitude types were found')
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
        for n in range(len(eventIDs[i])):
            outfile.write(
                eventIDs[i][n] + ',' + str(dependent_magnitudes[i][n]) + ',' + str(independent_magnitudes[i][n])
                + ',' + str(length[i][n]) + '\n')
process_end = datetime.datetime.now()
print('That took ' + str((process_end - process_start).total_seconds()) + ' seconds')