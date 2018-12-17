from obspy.io.quakeml.core import Unpickler
import pycurl
from io import BytesIO
import math
import matplotlib.pyplot as plt

quakeml_reader = Unpickler()

def FDSN_event_query(service, minmagnitude, minlongitude, maxlongitude,
                     minlatitude, maxlatitude, starttime='0000-01-01T00:00:00'):
    """
    Use obspy with pycurl to query event catalogs from FDSN members.

    :param service:
    :param minmagnitude:
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

    # Build query

    query = ""
    query = query.join((service, "query?", "minmagnitude=", str(minmagnitude),
                        "&minlatitude=", str(minlatitude), "&maxlatitude=", str(maxlatitude),
                        "&minlongitude=", str(minlongitude), "&maxlongitude=", str(maxlongitude),
                        "&starttime=", starttime))

    # Curl FDSN response

    buffer = BytesIO()
    c = pycurl.Curl()
    c.setopt(c.URL, query)
    c.setopt(c.WRITEDATA, buffer)
    c.perform()
    c.close()
    response = buffer.getvalue()

    # Parse quakeml bytes string through obspy

    catalog = quakeml_reader.loads(response)
    print(str(len(catalog)) + ' events were found in the catalog')

    return catalog


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

minmagnitude = 6  # minimum event magnitude to get from catalog
minlatitude, maxlatitude = -60, -13  # minimum and maximum latitude for event search window
minlongitude, maxlongitude = 145, -145  # minimum and maximum longitude for event search window

max_dt = 100  # maximum time (s) between events in separate catalogs for them to be considered records of the same earthquake
max_dist = 1000  # maximum distance (km) "

indep_mag_type = 'mww'  # first magnitude to compare (independent variable)
dep_mag_type = 'mB'  # second magnitude to compare (dependent variable)

# Build event catalogs

print('Searching earthquake catalogs for events above magnitude ' + str(minmagnitude) +
      ' between ' + str(minlatitude) + ' and ' + str(maxlatitude) + ' degrees latitude and ' +
      str(minlongitude) + ' and ' + str(maxlongitude) + ' degrees longitude')

# Build GeoNet event catalog

service = "https://service.geonet.org.nz/fdsnws/event/1/"
GeoNet_catalog = FDSN_event_query(service, minmagnitude, minlongitude, maxlongitude,
                                  minlatitude, maxlatitude)

# Build USGS event catalog

service = "https://earthquake.usgs.gov/fdsnws/event/1/"
USGS_catalog = FDSN_event_query(service, minmagnitude, minlongitude, maxlongitude,
                                minlatitude, maxlatitude)

# Match events between catalogs and extract magnitudes

dependent_magnitudes, independent_magnitudes = [], []
for n in range(len(GeoNet_catalog)):
    event = GeoNet_catalog[n]
    lengths = generate_lengths(event, USGS_catalog, max_dt, max_dist)
    print(str(len(lengths[0])) + ' reference events found for the ' + str(n) + 'th event in the catalog')
    if len(lengths[0]) > 0:
        matched_event = USGS_catalog[lengths[1][lengths[0].index(min(lengths[0]))]]
        for event_magnitude in event.magnitudes:
            if event_magnitude.magnitude_type == dep_mag_type:
                for matched_event_magnitude in matched_event.magnitudes:
                    if matched_event_magnitude.magnitude_type == indep_mag_type:
                        dependent_magnitudes.append(event_magnitude.mag)
                        independent_magnitudes.append(matched_event_magnitude.mag)

# Save data

print(str(len(dependent_magnitudes)) + ' matched magnitudes were found')

with open('magnitude_matches.csv', 'w') as outfile:
    outfile.write(indep_mag_type + ',' + dep_mag_type + '\n')
with open('magnitude_matches.csv', 'a') as outfile:
    for n in range(len(dependent_magnitudes)):
        outfile.write(str(independent_magnitudes[n]) + ',' + str(dependent_magnitudes[n]) + '\n')

# Plot data

plt.scatter(dependent_magnitudes, independent_magnitudes, s = 0.1)
plt.xlabel(dep_mag_type)
plt.ylabel(indep_mag_type)
plt.show()
