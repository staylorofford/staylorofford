"""
Functions to create a list file keys for a given start and end time over a given set of sites.
"""


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


def create_miniseed_key_list(site_codes, location_codes, channel_codes, start_time, end_time):

    """
    Create a list of all possible keys for day-long miniSEED files in the set of all site codes, location codes,
    and channel codes provided within the time set by the start and end times. No wildcard functionality exists.
    :param site_codes: list of site codes data is desired for
    :param location_codes: list of location codes data is desired for
    :param channel_codes: list of channel codes data is desired for
    :param start_time: datetime object describiing the desired data start date (UTC)
    :param end_time: as above, but for the end date
    :return: the list of all possible keys matching the input parameters
    """

    # Build the list of files using the standard naming conventions for miniseed files

    desired_keys = []
    for year in range(start_time.year, end_time.year + 1, 1):
        # Find the doy range to use
        if year == start_time.year and year != end_time.year:
            doy_min = start_time.timetuple().tm_yday
            if year % 4 == 0:
                doy_max = 366
            else:
                doy_max = 365
        elif year == start_time.year and year == end_time.year:
            doy_min = start_time.timetuple().tm_yday
            doy_max = end_time.timetuple().tm_yday
        elif year == end_time.year and year != start_time.year:
            doy_min = 1
            doy_max = end_time.timetuple().tm_yday
        else:
            doy_min = 1
            if year % 4 == 0:
                doy_max = 366
            else:
                doy_max = 365

        # Create desired keys for the current year
        for doy in range(doy_min, doy_max + 1, 1):
            doy = propend_zero(str(doy), 3)  # Ensure DOY is 0-padded appropriately
            for site in site_codes:
                for location_code in location_codes:
                    for channel in channel_codes:
                        channel = str.upper(channel)  # Ensure channel code is capitalised
                        desired_keys.append('miniseed/' + str(year) + '/' + str(year) + '.' + doy + '/' +
                                            site + '.NZ/' + str(year) + '.' + doy + '.' + site + '.' + location_code +
                                            '-' + channel + '.NZ.D')
    return desired_keys
