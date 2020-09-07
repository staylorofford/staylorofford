"""
Working script for seismic site data quality. Initially, this script will have functionality to:
1. Query data over a given time window for site codes in a csv file, on the vertical channel,
2. Filter data for a given filter(s),
3. Calculate data completeness (% samples) over the time window, and produce average data values using the RMS.
Note: while FDSN functionality was initially supported, further development has not continued to test this functionality.
"""

# Import Python libraries
import datetime
import io
import math
import matplotlib.pyplot as plt
import multiprocessing
import numpy as np
from obspy.clients.fdsn import Client
from obspy.clients.fdsn.header import FDSNNoDataException
from obspy.signal import PPSD
from obspy.core.stream import Stream
from obspy.core.utcdatetime import UTCDateTime
from obspy import read
from os.path import isfile
from os import remove
import pandas as pd
import pycurl


# Import my Python functions
import sys
sys.path.append('../../S3_tools/')
from assume_role import assume_role
from create_key_list import create_miniseed_key_list


def calculate_noise(data, filter_type, minimum_frequency, maximum_frequency, starttime, endtime):

    """
    Calculate average noise for a seismic site in the recorded units.
    :param data: obspy Stream object containing single trace to calculate data quality for
    :param filter_type: obspy filter type
    :param minimum_frequency: minimum frequency to use in filter
    :param maximum_frequency: maximum frequency to use in filter
    :param starttime: start time of data query, for FDSN completeness calculation, as ISO8601 string
    :param endtime: end time of data query, for FDSN completeness calculation, as ISO8601 string
    :return: data RMS value
    """

    # Apply a filter to the data
    if filter_type:
        data = data.filter(type=filter_type,
                           freqmin=minimum_frequency,
                           freqmax=maximum_frequency)

    # Build probabilistic power spectral density objects for each trace

    client = Client("https://service.geonet.org.nz")
    try:
        metadata = client.get_stations(network='NZ',
                                       station=data.stats.station,
                                       location=data.stats.location,
                                       channel=data.stats.channel,
                                       starttime=UTCDateTime(starttime),
                                       endtime=UTCDateTime(endtime),
                                       level='response')
        ppsd = PPSD(data.stats, metadata)
        ppsd.add(data)
        ppsd.save_npz(data.stats.station + '_' + data.stats.location + '_' + data.stats.channel + '_' +
                      starttime.isoformat() + '_PPSD.npz')
    except FDSNNoDataException:
        # When no response data exists, do not generate the PPSD
        pass

    # Scale data into nm/s
    gain = query_gain(data.stats.station, data.stats.location, data.stats.channel, starttime, endtime)
    try:
        gain = float(gain)
        data.data = gain * data.data

        # Calculate RMS value
        RMS = np.sqrt(np.mean(data.data ** 2))

    except ValueError:
        print('Failed to find sensor gain from FDSN! No scaling will be applied to RMS value for station ' +
              data.stats.station + ', location ' + data.stats.location + ', channel ' + data.stats.channel +
              ', times ' + starttime.isoformat() + '-' + endtime.isoformat() + '. Will return NaN value for RMS.')
        RMS = np.nan

    return RMS


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


def query_fdsn(station, location, channel, starttime, endtime, client="https://service.geonet.org.nz"):

    """
    Query FDSN for timeseries data. Input parameter can be specific or include ? wildcard
    :param: station: site code to get data from (string)
    :param: location: location code to get data from (string)
    :param: channel: channel code to get data from (string)
    :param: starttime: UTC start time of desired data (ISO8601 string)
    :param: endtime: UTC end time of desired data (ISO8601 string)
    :param: client: FDSN webservice to use as the client:
            https://service.geonet.org.nz (archive) or https://service-nrt.geonet.org.nz (last week)
    :return: obspy stream object containing requested data
    """

    client = Client(client)
    stream = client.get_waveforms(network='NZ',
                                  station=station,
                                  channel=channel,
                                  location=location,
                                  starttime=starttime,
                                  endtime=endtime)
    return stream


def query_gain(station, location, channel, starttime, endtime):

    """
    Query FDSN for site metadata.
    :param: station: site code to get data from (string)
    :param: location: location code to get data from (string)
    :param: channel: channel code to get data from (string)
    :param: starttime: UTC start time of desired data (ISO8601 string)
    :param: endtime: UTC end time of desired data (ISO8601 string)
    :param: client: FDSN webservice to use as the client:
            https://service.geonet.org.nz (archive) or https://service-nrt.geonet.org.nz (last week)
    :return: gain value for site
    """

    query_string = 'https://service.geonet.org.nz/fdsnws/station/1/query?station=' + station + \
                   '&location=' + str(location) + \
                   '&channel=' + str(channel) + \
                   '&starttime=' + starttime.isoformat() + \
                   '&endtime=' + endtime.isoformat() + \
                   '&level=channel'
    curl_result = curl(query_string).decode('ascii')
    sti = curl_result.find('<InstrumentSensitivity><Value>')
    ste = curl_result.find('</Value>', sti)
    gain = curl_result[sti + len('<InstrumentSensitivity><Value>'):ste]
    return gain


def data_quality(site_code, location_code, channel_code, latitude, longitude, starttime, endtime, splittime=None):

    try:
        if int(location_code) >= 20:
            return None
    except ValueError:
        return None

    print('Downloading data for station ' + site_code + ' location code ' + location_code + ' channel code '
          + channel_code + " over time " + starttime.isoformat() + ' - ' + endtime.isoformat() + '...')

    if data_source == 'FDSN':
        # Query data from FDSN
        if client == 'both':
            data_arc, data_nrt = None, None
            try:
                data_arc = query_fdsn(station=site_code,
                                      location=location_code,
                                      channel=channel_code,
                                      starttime=UTCDateTime(starttime),
                                      endtime=UTCDateTime(splittime),
                                      client='https://service.geonet.org.nz')
            except FDSNNoDataException:
                print('No data exists for: ')
                print('Site code: ' + site_code)
                print('Location code: ' + location_code)
                print('Channel code: ' + channel_code)
                print('Starttime: ' + starttime.isoformat())
                print('Endtime: ' + splittime.isoformat())
                print('')
            try:
                data_nrt = query_fdsn(station=site_code,
                                      location=location_code,
                                      channel=channel_code,
                                      starttime=UTCDateTime(splittime),
                                      endtime=UTCDateTime(endtime),
                                      client='https://serive-nrt.geonet.org.nz')
            except FDSNNoDataException:
                print('No data exists for: ')
                print('Site code: ' + site_code)
                print('Location code: ' + location_code)
                print('Channel code: ' + channel_code)
                print('Starttime: ' + splittime.isoformat())
                print('Endtime: ' + endtime.isoformat())
                print('')

            # Combine data
            if data_arc is not None and data_nrt is not None:
                data = data_arc + data_nrt
            elif data_arc is not None and data_nrt is None:
                data = data_arc
            elif data_arc is None and data_nrt is not None:
                data = data_nrt
            else:
                data = Stream()
        else:
            try:
                data = query_fdsn(station=site_code,
                                  location=location_code,
                                  channel=channel_code,
                                  starttime=UTCDateTime(starttime),
                                  endtime=UTCDateTime(endtime),
                                  client=client)
            except FDSNNoDataException:
                print('No data exists for: ')
                print('Site code: ' + site_code)
                print('Location code: ' + location_code)
                print('Channel code: ' + channel_code)
                print('Starttime: ' + starttime.isoformat())
                print('Endtime: ' + endtime.isoformat())
                print('')
                data = Stream()

        if len(data) == 0:
            # Save non-results

            outstr = (site_code + ',' +
                      location_code + ',' +
                      channel_code + ',' +
                      latitude + ',' +
                      longitude + ',' +
                      starttime.isoformat() + ',' +
                      endtime.isoformat() + ',' +
                      str(filter_type) + ',' +
                      str(minimum_frequency) + ',' +
                      str(maximum_frequency) + ',' +
                      str(np.nan) + '\n')
            return outstr

    elif data_source == 'S3':

        # Download data

        possible_keys = create_miniseed_key_list([site_code], [location_code], [channel_code], starttime, endtime)
        files = []
        for file_key in possible_keys:
            if isfile(file_key.split('/')[-1]):
                files.append(file_key.split('/')[-1])
                continue
            try:
                client.download_file('geonet-archive',
                                     file_key,
                                     file_key.split('/')[-1])
                print(file_key + ' downloaded successfully.')
                files.append(file_key.split('/')[-1])
            except:  # Fails if the file doesn't exist
                print(file_key + ' failed to be downloaded.')
                pass

    print('Data download complete. Calculating data quality...')

    # Calculate data on/off times, RMS values

    if data_source == 'FDSN':

        mean_RMS, _ = calculate_noise(data, filter_type, minimum_frequency, maximum_frequency, starttime, endtime)

    elif data_source == 'S3':

        # Here RMS values are taken as mean spectral intensity values in each band of a PPSD

        RMSs, times = [], []
        for file in files:
            print('Calculating data quality for file ' + file)
            data = read(file)
            # Iterative over data if there are gaps
            for data_chunk in data:
                RMS = calculate_noise(data_chunk, filter_type, minimum_frequency, maximum_frequency,
                                      starttime, endtime)
                RMSs.append(RMS)
                times.append([data_chunk.stats.starttime, data_chunk.stats.endtime])
            # remove(file)

        if len(RMSs) == 0:
            print('There was no data, moving onto next station.')

            # Save non-results

            outstr = (site_code + ',' +
                      location_code + ',' +
                      channel_code + ',' +
                      latitude + ',' +
                      longitude + ',' +
                      starttime.isoformat() + ',' +
                      endtime.isoformat() + ',' +
                      str(filter_type) + ',' +
                      str(minimum_frequency) + ',' +
                      str(maximum_frequency) + ',' +
                      str(np.nan) + '\n')
            return outstr

        else:

            # Save results

            outstr = ''
            for idx, time in enumerate(times):

                outstr += (site_code + ',' +
                           location_code + ',' +
                           channel_code + ',' +
                           latitude + ',' +
                           longitude + ',' +
                           time[0].isoformat() + ',' +
                           time[1].isoformat() + ',' +
                           str(filter_type) + ',' +
                           str(minimum_frequency) + ',' +
                           str(maximum_frequency) + ',' +
                           str(RMSs[idx])[:6] + '\n')

            print('Data quality calculation complete for data from site ' + site_code + ' location code ' +
                  location_code + ' channel code ' + channel_code + ' over time ' + starttime.isoformat() + ' - ' +
                  endtime.isoformat())
            return outstr


# Give the path to the site file: a csv containing in each row the station code and location code of each site to use
# in the data quality calculations. The station codes should be in the first column of this file, and the location
# code can be in any other column so long as the file header specifies which by the location of the "Location" string
# in the header.
site_file = '/home/samto/PROCESSING/sites.csv'

# Give parameters for which data to query
channel_codes = ['HHZ']  # Channel code(s) to query, if only one is desired have a list with one entry
starttime = '2019-01-01T00:00:00Z'  # starttime can either be an ISO8601 string, or an integer number of seconds to
                                    # query data for before the endtime.
endtime = '2019-12-31T23:59:59Z'  # endtime can be an ISO8601 string or 'now', if 'now' then endtime
                                  # will be set to the current UTC time. If using S3 data_source, make endtime at least
                                  # 1 second into the final day of data you want to use.

# Set whether to get data from FDSN (short duration) or S3. If S3 data will be downloaded for all full days between
# starttime and endtime inclusive.
data_source = 'S3'

# If you set data_source = 'S3', provide role details for S3 use
role = 'arn:aws:iam::862640294325:role/tf-sod-team-s3-read-role'
user = 'arn:aws:iam::582058524534:mfa/samto'
duration = 43200

# Set how to filter the data, if at all. Use filter_type=None to negate filtering. Filter types are those in obspy.
filter_type = None
minimum_frequency = 2
maximum_frequency = 15

# Parse parameters
if endtime == 'now':
    endtime = datetime.datetime.utcnow()
    if 'T' not in starttime:
        starttime = endtime - datetime.timedelta(seconds=int(starttime))
elif endtime != 'now' and 'T' not in starttime:
    endtime = datetime.datetime.strptime(endtime, '%Y-%m-%dT%H:%M:%SZ')
    starttime = endtime - datetime.timedelta(seconds=int(starttime))
else:
    endtime = datetime.datetime.strptime(endtime, '%Y-%m-%dT%H:%M:%SZ')
    starttime = datetime.datetime.strptime(starttime, '%Y-%m-%dT%H:%M:%SZ')

# Parse site code metadata from file
site_metadata = pd.read_csv(site_file, header=0, index_col=0)

# Prepare data query variables depending on data_source
if data_source == 'FDSN':
    # Determine which FDSN client to use
    splittime = datetime.datetime.utcnow() - datetime.timedelta(days=7)
    if splittime <= starttime and splittime <= endtime:
        # All data is within the nrt-client window
        client = 'https://service-nrt.geonet.org.nz'
    elif starttime < splittime and splittime <= endtime:
        # Query starts before nrt client window, but endtime is in nrt client window
        client = 'both'
    elif starttime < splittime and endtime < splittime:
        # All data is within the archive-client window
        client = 'https://service.geonet.org.nz'

elif data_source == 'S3':
    # Prepare the boto3 client to work with S3
    client = assume_role(role, user, duration)

else:
    print('data_source needs to be set to either \'FDSN\' or \'S3\' for the script to perform work.')
    exit()

# Pump parameters into parallel processing pool and do data quality calculation

# Prepare the workers for parallel processing
pp_pool = multiprocessing.Pool(3)

jobs = []

for l in range(math.ceil((endtime - starttime).total_seconds() / 86400)):
    # Increment the date
    date = starttime + datetime.timedelta(days=l)

    for m in range(len(site_metadata)):

        # Get site information
        site_code = str(site_metadata.index[m])
        location_code = str(site_metadata.iloc[m]['Location'])
        latitude = str(site_metadata.iloc[m]['Longitude'])
        longitude = str(site_metadata.iloc[m]['Latitude'])

        for n in range(len(channel_codes)):

            # Send jobs to the pool of workers
            job = pp_pool.apply_async(data_quality, (site_code, location_code, channel_codes[n], latitude, longitude,
                                                     date, date + datetime.timedelta(days=1, microseconds=-1)))
            jobs.append(job)

# Collect results from the workers through the pool result queue
with open('seismic_site_data_quality_output.csv', 'w') as outfile:
    # ^ IMPORTANT: have to use 'w' when working with parallel processes even though 'a' behaviour is desired!
    outfile.write('Station,Location,Channel,Longitude,Latitude,Starttime,Endtime,Filter,Freqmin,Freqmax,RMS\n')
    for job in jobs:
        result = job.get()
        try:
            outfile.write(result)
        except TypeError:
            # Fails on result=None, as desired.
            pass

# Close pool and join workers

pp_pool.close()
pp_pool.join()

# Aggregate results per site

with open('seismic_site_data_quality_summary.csv', 'w') as outfile:
    outfile.write('Station,Location,Channel,Longitude,Latitude,Starttime,Endtime,Filter,Freqmin,Freqmax,Completeness,'
                  'RMS\n')
with open('seismic_site_data_quality_summary.csv', 'a') as outfile:
    # Load in input file and parse it into rows for each channel at each site
    entries = []
    with open('seismic_site_data_quality_output.csv', 'r') as infile:
        for row in infile:
            if ',nan' not in row and 'Station' not in row:
                cols = row.split(',')
                entries.append(cols[0] + ',' + cols[1] + ',' + cols[2] + ',' + cols[3] + ',' + cols[4])
        entries = list(set(entries))
        completenesses = [0] * len(entries)
        RMSs = [0] * len(entries)
        numentries = [0] * len(entries)
    with open('seismic_site_data_quality_output.csv', 'r') as infile:
        for row in infile:
            if ',nan' not in row and 'Station' not in row:
                cols = row.split(',')
                idx = entries.index(cols[0] + ',' + cols[1] + ',' + cols[2] + ',' + cols[3] + ',' + cols[4])
                completenesses[idx] += (datetime.datetime.strptime(cols[6][:19], '%Y-%m-%dT%H:%M:%S') -
                                        datetime.datetime.strptime(cols[5][:19], '%Y-%m-%dT%H:%M:%S')
                                        ).total_seconds()
                RMSs[idx] += float(cols[-1])
                numentries[idx] += 1
        for n in range(len(entries)):
            completenesses[n] /= (endtime - starttime).total_seconds()
            if numentries[n] > 0:
                RMSs[n] /= numentries[n]
            else:
                RMSs[n] = np.nan
            outfile.write(entries[n] + ',' +
                          starttime.isoformat() + ',' +
                          endtime.isoformat() + ',' +
                          str(filter_type) + ',' +
                          str(minimum_frequency) + ',' +
                          str(maximum_frequency) + ',' +
                          str(completenesses[n]) + ',' +
                          str(RMSs[n]) + '\n')
