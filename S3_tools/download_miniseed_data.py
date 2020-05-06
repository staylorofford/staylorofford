"""
This script is to be used in the fulfilment of large data requests as per the
procedure for such data requests. It has two modes:
    - Find the size of files in the archive for a large data request
    - Perform the required transfer of files for a large data request
All input is given via command line arguments when the script is run.

In version 1.0 functionality only exists for archived miniseed data.
"""

import argparse
import datetime
import getpass
import sys
import queue
import threading

# Import my Python functions
import sys
sys.path.append('../')
from assume_role import assume_role
from create_key_list import create_miniseed_key_list


def copy_files(q):

    """
    Copy files to the geonet-data bucket
    :param q: queue object containing file prefixes
    :return: none
    """

    while True:

        client, bucket, file_key = q.get()

        try:
            client.download_file(bucket,
                                 file_key,
                                 file_key.split('/')[-1])
        except:  # Fails if the file doesn't exist
            pass

        q.task_done()


# 1) Parse arguments from command line

parser = argparse.ArgumentParser()
parser.add_argument('--data-format', type=str, help='Format of desired data')
parser.add_argument('--sites', type=str, help='Comma separated list of site codes of desired data')
parser.add_argument('--start', type=str, help='YYYY/DOY format UTC date for start time of desired data')
parser.add_argument('--end', type=str, help='YYYY/DOY format UTC date for end time of desired data')
parser.add_argument('--user', type=str, help='AWS username, default is the username of the script operator')
args = parser.parse_args()

# Parse file details

if not args.data_format:
    print('A root prefix is required! Set the root prefix using the --root-prefix argument. '
          'Run the script again with the -h flag for help.')
    sys.exit()
else:
    data_format = args.data_format
    if data_format != 'miniseed':
        print('This script only has functionality for files of miniseed format, '
              'please expand the code functionality if another format is desired.')
        sys.exit()
if not args.sites:
    print('No sites given! Use --sites to give a comma-separated list of sites containing the desired data. '
          'Run the script again with the -h flag for help.')
    sys.exit()
else:
    sites = args.sites.split(',')
    # Ensure site codes are capitalised
    for s in range(len(sites)):
        sites[s] = str.upper(sites[s])
if not args.start:
    print('No start date given! Use --start to give one in a YYYY/DOY format. '
          'Run the script again with the -h flag for help.')
    sys.exit()
else:
    start_time = datetime.datetime.strptime(args.start, '%Y/%j')
if not args.end:
    print('No end date given! Use --end to give one in a YYYY/DOY format. '
          'Use the -h flag for help.')
    sys.exit()
else:
    end_time = datetime.datetime.strptime(args.end, '%Y/%j')
bucket = 'geonet-archive'  # The read/write role only has access to this one bucket

# Get the user for use in the AWS authorisation

if args.user:
    user = args.user
else:
    user = getpass.getuser()

# 2) Build list of desired file keys

desired_keys = []
if data_format == 'miniseed':

    # Query the operator for data format specific details

    print('The operator will now be queried for the desired location codes and channel codes. '
          'Please put all possible codes desired, the script will attempt to download data for '
          'all codes for each site and will present details on which were successful once complete.')
    location_codes = input('As a comma-separated list, what location code(s) are desired? ').split(',')
    channels = input('As a comma-separated list, what channel(s) codes are desired? ').split(',')

    # Build the list of files using the standard naming conventions for miniseed files
    desired_keys = create_miniseed_key_list(sites, location_codes, channels, start_time, end_time)

# 3) Query bucket for desired files and their sizes

role = 'arn:aws:iam::862640294325:role/tf-sod-team-s3-read-role'
user = 'arn:aws:iam::582058524534:mfa/' + user
duration = 43200
client = assume_role(role, user, 43200)

q = queue.Queue()
num_threads = 100
for i in range(num_threads):
    worker = threading.Thread(target=copy_files, args=(q,))
    worker.setDaemon(True)
    worker.start()

for file_key in desired_keys:
    q.put([client, bucket, file_key])

q.join()
