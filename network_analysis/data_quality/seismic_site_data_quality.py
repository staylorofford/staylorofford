"""
Working script for seismic site data quality. Initially, this script will have functionality to:
1. Query data over a given time window for site codes in a csv file, on the vertical channel,
2. Filter data for a given filter(s),
3. Calculate data completeness (% samples) over the time window, and produce average data values using the RMS.
Note: while FDSN functionality was initially supported, further development has not continued to test this functionality.
"""

# Import Python libraries
import datetime
import math
import matplotlib.pyplot as plt
import multiprocessing
import numpy as np
from obspy.clients.fdsn import Client
from obspy.clients.fdsn.header import FDSNNoDataException
from obspy.core.stream import Stream
from obspy.core.utcdatetime import UTCDateTime
from obspy.signal import PPSD
from obspy import read
from os import remove
import pandas as pd


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
    except FDSNNoDataException:
        # When no response data exists
        return np.nan, np.nan

    # Find RMS value from PPSD.
    # 1) Take the mean value of PPSD in given frequency window as the PSD value
    # 2) Calculate weighted mean of PSD values in all windows using frequency window width as weights and scaling the
    # acceleration squared values by the window centre frequency squared to convert the result into velocity squared.
    # Also convert the data values out of dB scale as precursor to this.
    # 3) Take sqrt of weighted mean, as data are squared when processed to produce PSD. This gives the RMS value.
    weighted_mean, weight_sum = 0, 0

    try:
        _, mean_psds = ppsd.get_mean()
    except Exception:
        # Fails when no data exists
        return np.nan, np.nan

    psd_widths = [1/ppsd.period_bin_left_edges[n] - 1/ppsd.period_bin_right_edges[n]
                  for n in range(len(ppsd.period_bin_left_edges))]
    psd_centres = [1/ppsd.period_bin_centers[n] for n in range(len(ppsd.period_bin_centers))]
    for n in range(len(mean_psds)):
        weighted_mean += math.sqrt(10**(mean_psds[n] / 10) / (psd_centres[n] ** 2)) * psd_widths[n]
        weight_sum += psd_widths[n]
    weighted_mean /= weight_sum

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
                      filter_type + ',' +
                      str(minimum_frequency) + ',' +
                      str(maximum_frequency) + ',' +
                      'NA' + ',' +
                      str(np.nan) + '\n')
            return outstr

    elif data_source == 'S3':

        # Download data

        possible_keys = create_miniseed_key_list([site_code], [location_code], [channel_code], starttime, endtime)
        files = []
        for file_key in possible_keys:
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

        mean_RMS, _ = calculate_PPSD_noise(data, filter_type, minimum_frequency, maximum_frequency, starttime, endtime)

    elif data_source == 'S3':

        # Here RMS values are taken as mean spectral intensity values in each band of a PPSD

        RMSs, times = [], []
        for file in files:
            print('Calculating data quality for file ' + file)
            data = read(file)
            # Iterative over data if there are gaps
            for data_chunk in data:
                RMS, ppsd = calculate_PPSD_noise(data_chunk, filter_type, minimum_frequency, maximum_frequency,
                                                 starttime, endtime)
                if np.isnan(RMS):
                    continue
                RMSs.append(RMS)
                chunk_st, chunk_et = ppsd.times_data[0]
                times.append([chunk_st, chunk_et])
            remove(file)

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
                      filter_type + ',' +
                      str(minimum_frequency) + ',' +
                      str(maximum_frequency) + ',' +
                      '0' + ',' +
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
                           filter_type + ',' +
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
# site_file = '/home/samto/git/delta/network/sites.csv'
site_file = './sites'
# site_file = './sites.csv'

# Give parameters for which data to query
channel_codes = ['EHZ', 'HHZ']  # Channel code(s) to query, if only one is desired have a list with one entry
starttime = '2019-03-06T00:00:00Z'  # starttime can either be an ISO8601 string, or an integer number of seconds to
                                    # query data for before the endtime.
endtime = '2020-12-31T00:00:01Z'  # endtime can be an ISO8601 string or 'now', if 'now' then endtime
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
    all_site_rows = [[[] for m in range(len(channel_codes))] for n in range(len(site_metadata))]
    for n, site in enumerate(site_metadata.index):
        for m, channel in enumerate(channel_codes):
            with open('seismic_site_data_quality_output.csv', 'r') as infile:
                for row in infile:
                    cols = row.split(',')
                    if cols[0] == site and cols[2] == channel:
                        all_site_rows[n][m].append(row)
        # Then, aggregate all the different rows into a summary row for each channel at each site. This collapses
        # the data by time.
        channel_rows = [[] for m in range(len(channel_codes))]
        for site_rows in all_site_rows:
            for channel_rows in site_rows:
                try:
                    first_row = channel_rows[0].split(',')
                except IndexError:
                    continue
                site = first_row[0]
                location_code = first_row[1]
                longitude = first_row[3]
                latitude = first_row[4]
                for channel in channel_codes:
                    mean_RMS = 0
                    mean_completeness = 0
                    numentries = 0
                    for row in channel_rows:
                        cols = row.split(',')
                        if cols[2] == channel:
                            if cols[-1][:3] != 'nan':
                                mean_RMS += float(cols[-1])
                                mean_completeness += (datetime.datetime.strptime(cols[6][:19], '%Y-%m-%dT%H:%M:%S') -
                                                      datetime.datetime.strptime(cols[5][:19], '%Y-%m-%dT%H:%M:%S')
                                                      ).total_seconds()
                                numentries += 1
                    if numentries > 0:
                        mean_RMS /= numentries
                    else:
                        mean_RMS = np.nan
                    mean_completeness /= (endtime - starttime).total_seconds()
                    outfile.write(site + ',' +
                                  location_code + ',' +
                                  channel + ',' +
                                  longitude + ',' +
                                  latitude + ',' +
                                  starttime.isoformat() + ',' +
                                  endtime.isoformat() + ',' +
                                  filter_type + ',' +
                                  str(minimum_frequency) + ',' +
                                  str(maximum_frequency) + ',' +
                                  str(mean_completeness) + ',' +
                                  str(mean_RMS) + '\n')
