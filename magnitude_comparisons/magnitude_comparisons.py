from obspy.io.quakeml.core import Unpickler
quakeml_reader = Unpickler()
import pycurl
from io import BytesIO

def FDSN_event_query(service, minmagnitude, minlongitude, maxlongitude,
                     minlatitude, maxlatitude, starttime='0000-01-01T00:00:00'):

    """
    Use obspy with pycurl to query event catalogs from FDSN members.
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




# Set script parameters

minmagnitude = 8 # minimum event magnitude to get from catalog
minlatitude, maxlatitude = -60, -13 # minimum and maximum latitude for event search window
minlongitude, maxlongitude = 145, -145  # minimum and maximum longitude for event search window

indep_mag_type = 'Mw' # first magnitude to compare (independent variable)
dep_mag_type = 'MLv' # second magnitude to compare (dependent variable)

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

# Match events between catalogs

for event in GeoNet_catalog:
    print(event.origins[0])
    interevent_distances = []
    interevent_times = []
    for reference_event in USGS_catalog:
        pass
    break
print('More to come...')