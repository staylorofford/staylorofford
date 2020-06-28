from obspy.clients.fdsn import Client
from obspy.core.utcdatetime import UTCDateTime

# User-defined parameters of script
client = Client("https://service-nrt.geonet.org.nz")  # FDSN client to query for station information and responses
stations = ['WIZ', 'WSRZ']  # Give station codes as strings in a list
desired_locations = ['10']  # Set desired location codes to get response for, if all are wanted set to None
desired_channels = ['HHZ']  # As above, but for channels
starttime = '2019-09-01T00:00:00Z'  # Give start time as ISO8601 string
endtime = '2019-10-01T00:00:00Z'  # Give end time as ISO8601 string

# Prepare data lists for querying
location_codes = [[] for n in range(len(stations))]
channel_codes = [[] for n in range(len(stations))]
responses = [[] for n in range(len(stations))]
for n, station in enumerate(stations):
    # Get site details at station
    channels = client.get_stations(network='NZ',
                                   station=station,
                                   starttime=UTCDateTime(starttime),
                                   endtime=UTCDateTime(endtime),
                                   level='channel')[0][0].channels
    for channel in channels:

        if desired_locations and channel._location_code in desired_locations and \
                desired_channels and channel._code in desired_channels:

            # Query channel response for site
            response = client.get_stations(network='NZ',
                                           station=station,
                                           location=channel._location_code,
                                           channel=channel._code,
                                           starttime=UTCDateTime(starttime),
                                           endtime=UTCDateTime(endtime),
                                           level='response')[0][0].channels[0].response

            # Save response data
            if channel._location_code not in location_codes[n]:
                location_codes[n].append(channel._location_code)
                channel_codes[n].append([channel._code])
                responses[n].append([response])
            else:
                m = location_codes[n].index(channel._location_code)
                channel_codes[n][m].append(channel._code)
                responses[n][m].append(response)

# Return response storage information to operator
print(location_codes)
print(channel_codes)
print(responses)
# The responses in the responses lists are of the Response class in obspy
# (https://docs.obspy.org/master/packages/autogen/obspy.core.inventory.response.Response.html)
# and to use them to remove a stream parsed using obspy.read(), the remove_response() method should be used on each
# trace in the stream (https://docs.obspy.org/packages/autogen/obspy.core.trace.Trace.remove_response.html).
