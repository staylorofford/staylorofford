"""
Use Dark Sky historical weather data to investigate weather relationships with data recorded at a given site.
"""

import datetime as dt
import math
import pandas as pd

def query_API(key, latitude, longitude, time):

    """
    Query the Dark Sky API for weather data and use pandas to parse the result into a dataframe.
    :param key: API secret key
    :param latitude: latitude of point to sample
    :param longitude: longitude of point to sample
    :param time: time of point to sample
    :return: dataframe containing hourly data at sampled point
    """

    if isinstance(time, str):

        # Get data for UTC day of time
        query = 'https://api.darksky.net/forecast/' + key + '/' + latitude + ',' + longitude + ',' + time
        queryresult = pd.read_json(query)
        all_data = queryresult['hourly']['data']
        for m, hour_data in enumerate(all_data):
            if m == 0:
                df = pd.DataFrame(hour_data, index=[m])
            else:
                df = pd.concat([df, pd.DataFrame(hour_data, index=[m])])
        return df

    elif isinstance(time, list):

        # Build data for first day
        query = 'https://api.darksky.net/forecast/' + key + '/' + latitude + ',' + longitude + ',' + time[0]
        queryresult = pd.read_json(query)
        all_data = queryresult['hourly']['data']
        for m, hour_data in enumerate(all_data):
            if m == 0:
                df = pd.DataFrame(hour_data, index=[m])
            else:
                df = pd.concat([df, pd.DataFrame(hour_data, index=[m])], sort=True)

        # Prepare time variables for all days
        st = dt.datetime.strptime(time[0], '%Y-%m-%dT%H:%M:%SZ')
        et = dt.datetime.strptime(time[1], '%Y-%m-%dT%H:%M:%SZ')
        doys = int(math.ceil((et - st).total_seconds() / 86400))

        # Build data for all days
        for doy in range(doys):
            time = (st + dt.timedelta(days=doy)).isoformat()
            query = 'https://api.darksky.net/forecast/' + key + '/' + latitude + ',' + longitude + ',' + time
            queryresult = pd.read_json(query)
            all_data = queryresult['hourly']['data']
            for n, hour_data in enumerate(all_data):
                m += n + 1
                df = pd.concat([df, pd.DataFrame(hour_data, index=[m])], sort=True)
        return df


# Give parameters for query: key, latitude, and longitude are single strings, but time can be a single string or
# a list of two strings framing the period to sample data over.
key = '5d40a297fa26ebd6b5e6777b7c7bc54c'  # Your secret API key
latitude = '-36.745228533'  # The latitude to sample data at
longitude = '175.720872942'  # The longitude to sample data at
# time = '2019-04-01T00:00:00Z'  # Time specifies which UTC day to return hourly data for
time = ['2019-04-03T00:00:00Z', '2019-05-05T00:00:00Z']

data = query_API(key, latitude, longitude, time)
print(data.to_string())