"""
Working script for seismic site data quality. Initially, this script will have functionality to:
1. Query data over a given time window for site codes in a csv file, on the vertical channel,
2. Filter data for a given filter(s),
3. Calculate data completeness (% samples) over the time window, and produce average data values using the RMS.
"""

# Import Python libraries
import datetime
import math
import matplotlib.pyplot as plt
import numpy as np
from obspy.clients.fdsn import Client
from obspy.clients.fdsn.header import FDSNNoDataException
from obspy.core.stream import Stream
from obspy.core.utcdatetime import UTCDateTime
from obspy import read
from os import listdir
from os import remove
import pandas as pd
import queue
import threading

# Import my Python functions
import sys
sys.path.append('../../S3_tools/')
from assume_role import assume_role
from create_key_list import create_miniseed_key_list


def calculate_data_quality(data, filter_type, minimum_frequency, maximum_frequency, data_source, starttime, endtime):

    """
    Calculate data quality via data completeness and RMS value.
    :param data: obspy Stream object containing single trace to calculate data quality for
    :param filter_type: obspy filter type
    :param minimum_frequency: minimum frequency to use in filter
    :param maximum_frequency: maximum frequency to use in filter
    :param data_source: where the data is from, FDSN or S3: sets how to calculate completeness
    :param starttime: start time of data query, for FDSN completeness calculation, as ISO8601 string
    :param endtime: end time of data query, for FDSN completeness calculation, as ISO8601 string
    :return: data completeness as a percentage, data RMS value, number of data values
    """

    # Subtract a running mean of 10 second duration from the data to remove the influence of long periods
    if (data[0].stats.endtime - data[0].stats.starttime) < 10:
        running_mean = np.nanmean(data[0].data)
        for n in range(len(data[0].data)):
            data[0].data[n] -= running_mean
    else:
        running_mean = []
        for n in range(5 * int(data[0].stats.sampling_rate),
                       len(data[0].data) - 5 * int(data[0].stats.sampling_rate)):
            running_mean.append(np.nanmean(data[0].data[n - 5 * int(data[0].stats.sampling_rate):
                                                     n + 5 * int(data[0].stats.sampling_rate)]))
        for n in range(5 * int(data[0].stats.sampling_rate),
                       len(data[0].data) - 5 * int(data[0].stats.sampling_rate)):
            data[0].data[n] -= running_mean[n - 5 * int(data[0].stats.sampling_rate)]
        # Remove data that was not adjusted by the running mean subtraction
        data[0].data = data[0].data[5 * int(data[0].stats.sampling_rate):
                                    len(data[0].data) - 5 * int(data[0].stats.sampling_rate)]

    # Apply a filter to the data
    if filter_type:
        data = data.filter(type=filter_type,
                           freqmin=minimum_frequency,
                           freqmax=maximum_frequency)

    # Calculate the RMS value for the data
    RMS = 0
    data_N = 0
    mean = np.nanmean(np.abs(data[0].data))
    for n in range(len(data[0].data)):
        if not np.isnan(data[0].data[n]):
            RMS += math.pow((mean - data[0].data[n]), 2)
            data_N += 1
    if data_N != 0:
        RMS = math.sqrt(RMS / data_N)
    else:
        RMS = np.nan

    # Calculate completeness
    if data_source == 'FDSN':
        expected_data = (endtime - starttime).total_seconds() * data[0].stats.sampling_rate
        if expected_data > 0:
            completeness = 100 * data_N / expected_data
        else:
            completeness = 0
    elif data_source == 'S3':
        expected_data = data[0].stats.sampling_rate * 86400
        completeness = data_N / expected_data

    print('For data:')
    print(data)
    print('Completeness is: ' + str(completeness))
    print('RMS is: ' + str(RMS))
    print('N is: ' + str(data_N) + '\n')

    return completeness, RMS, data_N


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


def copy_files(q):

    """
    Copy files to the geonet-data bucket
    :param q: queue object containing file prefixes
    :return: none
    """

    while True:

        client, bucket, file_key = q.get()

        try:
            client.download_file('geonet-archive',
                                 file_key,
                                 file_key.split('/')[-1])
            print(file_key + ' downloaded successfully.')
        except:  # Fails if the file doesn't exist
            print(file_key + ' failed to be downloaded.')
            pass

        q.task_done()


# Give the path to the site file: a csv containing in each row the station code and location code of each site to use
# in the data quality calculations. The station codes should be in the first column of this file, and the location
# code can be in any other column so long as the file header specifies which by the location of the "Location" string
# in the header.
site_file = '/home/samto/git/delta/network/sites.csv'

# Give parameters for which data to query
channel_code = 'HHZ'  # Channel code to query
starttime = '2019-01-01T00:00:00Z'  # starttime can either be an ISO8601 string, or an integer number of seconds to
                                    # query data for before the endtime.
endtime = '2019-01-02T00:00:01Z'  # endtime can be an ISO8601 string or 'now', if 'now' then endtime
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
filter_type = 'bandpass'
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

# Prepare output file
with open('seismic_site_data_quality_output.csv', 'w') as outfile:
    outfile.write('Station,Longitude,Latitude,Starttime,Endtime,Filter,Freqmin,Freqmax,RMS,Mean_RMS,Completeness,'
                  'Mean_Completeness\n')

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

for m in range(len(site_metadata)):
    # Get site information
    site_code = site_metadata.index[m]
    location_code = site_metadata.iloc[m]['Location']

    # Hardcoded skip of strong motion sites
    try:
        if int(location_code) >= 20:
            continue
    except ValueError:
        continue

    print('Downloading data for station ' + site_code + ' location code ' + location_code + '...')

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
            continue

    elif data_source == 'S3':

        possible_keys = create_miniseed_key_list([site_code], [location_code], [channel_code], starttime, endtime)

        q = queue.Queue()
        num_threads = 100
        for i in range(num_threads):
            worker = threading.Thread(target=copy_files, args=(q,))
            worker.setDaemon(True)
            worker.start()

        for file_key in possible_keys:
            q.put([client, 'geonet-archive', file_key])

        q.join()

    print('Data download complete. Calculating data quality...')

    # Calculate data completeness, RMS value

    if data_source == 'FDSN':

        completeness, RMS, _ = calculate_data_quality(data, filter_type, minimum_frequency, maximum_frequency,
                                                      data_source, starttime, endtime)
        mean_RMS = np.nan  # There is no mean RMS using this source: the only RMS value is the one for the whole dataset

    elif data_source == 'S3':

        # Here completeness and RMS values are taken as averages of daily values. We assume that the number of
        # samples in each day is approximately the same, making this a statistically valid approach. Daily averaging
        # is done due to the large data volumes encountered using this data source.

        files = listdir('.')
        completenesses, RMSs, Ns = [], [], []
        for file in files:
            if file[-1:] == 'D' and site_code in file and location_code in file:
                print('Calculating data quality for file ' + file)
                data = read(file)
                daily_completeness, daily_RMS, daily_data_number = calculate_data_quality(data, filter_type,
                                                                                          minimum_frequency,
                                                                                          maximum_frequency,
                                                                                          data_source,
                                                                                          starttime, endtime)
                completenesses.append(daily_completeness)
                RMSs.append(daily_RMS)
                Ns.append(daily_data_number)
                remove(file)
        if len(completenesses) == 0:
            print('There was no data, moving onto next station.')
            continue

        # Calculate mean completeness, RMS
        mean_completeness = np.mean(completenesses)
        mean_RMS = np.nanmean(RMSs)  # The mean RMS is the mean of daily RMS values
        # Calculate overall completeness, RMS
        completeness = math.ceil((endtime - starttime).total_seconds() / 86400) * 86400 * \
                       data[0].stats.sampling_rate / sum(Ns)
        RMS = 0
        for n in range(len(RMSs)):
            if not np.isnan(RMSs[n]):
                RMS += math.pow(RMSs[n], 2) * Ns[n]
        # The overall RMS value is the sum of all scaled (by number of data in each respective day of data) squared
        # daily RMS values divided by the total amount of data, to the power of 1/2.
        # If the number of samples per day is roughly the same over the time window then the two RMS values should
        # be approximately equal. If that number varies, so too will the RMS values.
        RMS = math.sqrt(RMS / sum(Ns))
        print('Data quality calculation complete for data:')
        print(data)
        print('Mean RMS is: ' + str(mean_RMS))
        print('Reconstituted RMS is: ' + str(RMS))
        print('Mean completeness is: ' + str(mean_completeness))
        print('Overall completeness is: ' + str(completeness) + '\n')

    # Save results

    with open('seismic_site_data_quality_output.csv', 'a') as outfile:
        outfile.write(site_code + ',' +
                      str(site_metadata.iloc[m]['Longitude']) + ',' +
                      str(site_metadata.iloc[m]['Latitude']) + ',' +
                      starttime.isoformat() + ',' +
                      endtime.isoformat() + ',' +
                      filter_type + ',' +
                      str(minimum_frequency) + ',' +
                      str(maximum_frequency) + ',' +
                      str(RMS)[:6] + ',' +
                      str(mean_RMS)[:6] + ',' +
                      str(completeness)[:6] + ',' +
                      str(mean_completeness)[:6] + '\n')
