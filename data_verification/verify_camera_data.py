"""
Script to verify data from cameras over an equipment change.
"""

import datetime as dt
import io
import matplotlib.pyplot as plt
import pycurl


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


def propend_zero(string, desired_length):

    """
    Add a 0 to the front of a string if it is less than desired_length in length.
    :param string: any string; the string to add 0s to the front of
    :param desired_length: any integer; the length of the output string desired
    :return: The input string or the input string propended by 0s
    """

    while len(string) < desired_length:
        string = '0' + string
    return string


archive_root = 'ftp://ftp.geonet.org.nz/volcano_cams/'  # Root directory of data

# Give parameters for data verification
station = 'T'
session_endtime = '2020-02-18T00:40:01Z'  # End time of previous session, from delta
session_starttime = '2020-02-27T22:59:59Z'  # Start time of current session, from delta

# Parse parameters
set_dt = dt.datetime.strptime(session_endtime, '%Y-%m-%dT%H:%M:%SZ')
sst_dt = dt.datetime.strptime(session_starttime, '%Y-%m-%dT%H:%M:%SZ')

# Get start times for querying
st = set_dt - dt.timedelta(weeks=1)
et = sst_dt + dt.timedelta(weeks=1)

# Query file names and sizes from archive
ct = st
existing_files = []
data_times = []
file_sizes = []
while ct <= et:
    try:
        archive_files = curl(archive_root + str(ct.year) + '/' + str(ct.year) + propend_zero(str(ct.month), 2) +
                             propend_zero(str(ct.day), 2) + '/' + station + '/').decode('ascii').split('\n')
    except pycurl.error:
        # Occurs when directory does not exist, i.e. there is no data for the station on the current day
        ct += dt.timedelta(days=1)
        continue
    for file in archive_files:
        file_details = file.split()
        try:
            existing_files.append(file_details[8])
        except IndexError:
            # Empty row in curl result
            continue
        data_times.append(dt.datetime.strptime(file_details[8][:-5], '%Y%m%d%H%M%S'))
        file_sizes.append(file_details[4])
    ct += dt.timedelta(days=1)

# Return number of files pre- and post-change to operator
presum, postsum = 0, 0
for data_time in data_times:
    if data_time <= set_dt:
        presum += 1
    elif data_time >= sst_dt:
        postsum += 1
print('Number of files in the week pre-change is: ' + str(presum))
print('Number of files in the week post-change is: ' + str(postsum))

# Plot data existence timeseries
fig = plt.figure()
plt.scatter(data_times, [1] * len(data_times), s=1)
plt.xlim(st, et)
fig.autofmt_xdate()
plt.show()
