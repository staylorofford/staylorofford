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
from obspy.signal import PPSD
from obspy import read
from os import listdir
from os import remove
import pandas as pd
import queue
import threading
import time

# Import my Python functions
import sys
sys.path.append('../../S3_tools/')
from assume_role import assume_role
from create_key_list import create_miniseed_key_list


def calculate_PPSD_noise(data, filter_type, minimum_frequency, maximum_frequency, starttime, endtime):

    """
    Calculate data quality via data completeness and RMS value.
    :param data: obspy Stream object containing single trace to calculate data quality for
    :param filter_type: obspy filter type
    :param minimum_frequency: minimum frequency to use in filter
    :param maximum_frequency: maximum frequency to use in filter
    :param starttime: start time of data query, for FDSN completeness calculation, as ISO8601 string
    :param endtime: end time of data query, for FDSN completeness calculation, as ISO8601 string
    :return: data RMS value, number of data values
    """

    pst = time.perf_counter()

    # Apply a filter to the data
    if filter_type:
        data = data.filter(type=filter_type,
                           freqmin=minimum_frequency,
                           freqmax=maximum_frequency)

    # Build probabilistic power spectral density objects for each trace
    client = Client("https://service.geonet.org.nz")
    metadata = client.get_stations(network='NZ',
                                   station=data[0].stats.station,
                                   location=data[0].stats.location,
                                   channel=data[0].stats.channel,
                                   starttime=UTCDateTime(starttime),
                                   endtime=UTCDateTime(endtime),
                                   level='response')
    ppsd = PPSD(data[0].stats, metadata)
    ppsd.add(data[0])

    # Find RMS value from PPSD.
    # 1) Take the mean value of PPSD in given frequency window as the PSD value
    # 2) Calculate weighted mean of PSD values in all windows using frequency window width as weights and scaling the
    # acceleration squared values by the window centre frequency squared to convert the result into velocity squared.
    # Also convert the data values out of dB scale as precursor to this.
    # 3) Take sqrt of weighted mean, as data are squared when processed to produce PSD. This gives the RMS value.
    weighted_mean, weight_sum = 0, 0
    _, mean_psds = ppsd.get_mean()
    psd_widths = [1/ppsd.period_bin_left_edges[n] - 1/ppsd.period_bin_right_edges[n]
                  for n in range(len(ppsd.period_bin_left_edges))]
    psd_centres = [1/ppsd.period_bin_centers[n] for n in range(len(ppsd.period_bin_centers))]
    for n in range(len(mean_psds)):
        weighted_mean += math.sqrt(10**(mean_psds[n] / 10) / (psd_centres[n] ** 2)) * psd_widths[n]
        weight_sum += psd_widths[n]
    weighted_mean /= weight_sum
    pet = time.perf_counter()
    print('Processing that data took ' + str((pet - pst) / 60.0) + ' minutes.')

    return weighted_mean, ppsd


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
endtime = '2019-01-05T00:00:01Z'  # endtime can be an ISO8601 string or 'now', if 'now' then endtime
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
    outfile.write('Station,Longitude,Latitude,Starttime,Endtime,Filter,Freqmin,Freqmax,RMS\n')

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

    # Prepare the workers for parallel processing
    q = queue.Queue()
    num_threads = 100
    for i in range(num_threads):
        worker = threading.Thread(target=copy_files, args=(q,))
        worker.setDaemon(True)
        worker.start()

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
            # Save non-results

            with open('seismic_site_data_quality_output.csv', 'a') as outfile:
                outfile.write(site_code + ',' +
                              str(site_metadata.iloc[m]['Longitude']) + ',' +
                              str(site_metadata.iloc[m]['Latitude']) + ',' +
                              starttime.isoformat() + ',' +
                              endtime.isoformat() + ',' +
                              filter_type + ',' +
                              str(minimum_frequency) + ',' +
                              str(maximum_frequency) + ',' +
                              str(np.nan) + ',' +
                              str(np.nan) + ',' +
                              str(np.nan) + ',' +
                              str(np.nan) + '\n')
            continue

    elif data_source == 'S3':

        possible_keys = create_miniseed_key_list([site_code], [location_code], [channel_code], starttime, endtime)

        for file_key in possible_keys:
            q.put([client, 'geonet-archive', file_key])

        q.join()

    print('Data download complete. Calculating data quality...')

    # Calculate data completeness, RMS value

    if data_source == 'FDSN':

        mean_RMS, _ = calculate_PPSD_noise(data, filter_type, minimum_frequency, maximum_frequency,
                                                       starttime, endtime)

    elif data_source == 'S3':

        # Here completeness and RMS values are taken as averages of daily values. We assume that the number of
        # samples in each day is approximately the same, making this a statistically valid approach. Daily averaging
        # is done due to the large data volumes encountered using this data source.

        files = listdir('.')
        RMSs, ppsds = [], []
        for file in files:
            if file[-1:] == 'D' and site_code in file and location_code in file:
                print('Calculating data quality for file ' + file)
                data = read(file)
                daily_RMS, ppsd = calculate_PPSD_noise(data, filter_type, minimum_frequency, maximum_frequency,
                                                       starttime, endtime)
                RMSs.append(daily_RMS)
                # PPSD: need to find a way to compare the noise value returned from this (which is what, a mean?)
                # to the RMS value determined from timeseries data. May require some digging into literature, but
                # have timeseries mean value for comparison.
                ppsds.append(ppsd)
                remove(file)
        if len(RMSs) == 0:
            print('There was no data, moving onto next station.')
            # Save non-results

            with open('seismic_site_data_quality_output.csv', 'a') as outfile:
                outfile.write(site_code + ',' +
                              str(site_metadata.iloc[m]['Longitude']) + ',' +
                              str(site_metadata.iloc[m]['Latitude']) + ',' +
                              starttime.isoformat() + ',' +
                              endtime.isoformat() + ',' +
                              filter_type + ',' +
                              str(minimum_frequency) + ',' +
                              str(maximum_frequency) + ',' +
                              str(np.nan) + ',' +
                              str(np.nan) + ',' +
                              str(np.nan) + ',' +
                              str(np.nan) + '\n')
            continue

        # Calculate mean completeness, RMS
        mean_RMS = np.nanmean(RMSs)  # The mean RMS is the mean of daily RMS values
        print('Data quality calculation complete for data from site ' + site_code + ' location code ' +
              location_code + ' channel code ' + channel_code)
        print('Mean RMS is: ' + str(mean_RMS))

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
                      str(mean_RMS)[:6] + '\n')
