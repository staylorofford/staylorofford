"""
Transform GeoNet delta network tables to series of gis-ready files using
data from the GeoNet API to distinguish sensor types so that a network map
like that available on the Geonet website can be made.
Note: only  data for open sites is saved to file.
"""


import io
import os
import pycurl
import pandas

# Define sensor types in GeoNet API query by number and the name desired for that sensor type
site_sensor_types = [['1,10', 'Strong Motion Sensor'],
                     ['2,6', 'Air Pressure Sensor'],
                     ['3', 'Broadband Seismometer'],
                     ['7', 'Coastal Sea Level Gauge'],
                     ['8,9', 'Short Period Seismometer']]


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


def generate_sensor_types():

    """
    Create a nested list of sensor types for all sensors at sites
    :return: list of sensor types, nested list of site codes, loc codes, and associated sensor types
    """

    # Get site code, loc code for each currently installed network sensor at a station
    sensor_types = []
    site_sensor_type = []
    for network_sensor in site_sensor_types:
        sensor_types.append(network_sensor[1])
        query = 'https://api.geonet.org.nz/network/sensor?sensorType=' + network_sensor[0] + '&endDate=9999-01-01'
        network_sensor_data = pandas.read_json(curl(query))
        for sensor_data in network_sensor_data['features']:
            sensor_data = sensor_data['properties']
            site_sensor_type.append([sensor_data['Station'],
                                     sensor_data['Location'],
                                     network_sensor[1]])
    sensor_types = list(set(sensor_types))

    return sensor_types, site_sensor_type


def parse_csv(csv_path):

    """
    Parse the contents of CSV file at csv_path into nested lists
    :param csv_path: full path to CSV file
    :return: list of CSV header columns, nested lists of CSV row data
    """

    data = []
    with open(csv_path, 'r') as openfile:
        rc = -1
        for row in openfile:
            if rc == -1:
                header = row[:-1].split(',')
                rc = 0
            else:
                data.append(row[:-1].split(','))
    return header, data


def get_network_data(delta_root_path):

    """
    Parse all network data from GeoNet delta and return it in nested lists split by site type.
    :param delta_root_path: absolute path to local GeoNet delta git repository root
    :return: list of site types data was split into and nested list of network data and network data headers
    """

    # Build file paths
    if delta_root_path[-1] != '/':
        delta_root_path += '/'
    marks = delta_root_path + 'network/marks.csv'
    mounts = delta_root_path + 'network/mounts.csv'
    sites = delta_root_path + 'network/sites.csv'

    # Get data from delta network files
    marks_header, marks_data = parse_csv(marks)
    mounts_header, mounts_data = parse_csv(mounts)
    sites_header, sites_data = parse_csv(sites)

    # Parse sensor type data from GeoNet API
    sensor_types, site_sensor_type = generate_sensor_types()

    # Decompose station data into site types using station and site data
    # Plus sensor type data from the GeoNet API
    sites_header.append('Sensor Type')
    for m in range(len(sites_data)):
        for n in range(len(site_sensor_type)):
            # If the site code and location code match
            if sites_data[m][0] == site_sensor_type[n][0] and sites_data[m][1] == site_sensor_type[n][1]:
                # Add the sensor type data to the sites_data list
                sites_data[m].append(site_sensor_type[n][2])
                break  # Avoid duplicates

    # Add to other delta data a Sensor Type column
    marks_header.append('Sensor Type')
    for m in range(len(marks_data)):
        marks_data[m].append('GNSS')
    sensor_types.append('GNSS')
    mounts_header.append('Sensor Type')
    for m in range(len(mounts_data)):
        mounts_data[m].append('Camera')
    sensor_types.append('Camera')

    # Group data into lists
    network_data_headers = [marks_header, mounts_header, sites_header]
    network_data = [marks_data, mounts_data, sites_data]

    return sensor_types, network_data_headers, network_data


if __name__ == "__main__":

    # Parse data from config file
    parameters = ['delta_root_path']
    values = [0] * len(parameters)
    config_file = os.path.dirname(os.path.realpath(__file__)) + '/config-delta_to_gis.txt'
    with open(config_file, 'r') as openfile:
        for row in openfile:
            if row.split(' ')[0] in parameters:
                idx = parameters.index(row.split(' ')[0])
                values[idx] = row.split(' ')[2]
    print('Input parameter values are: ')
    for n in range(len(parameters)):
        print(parameters[n] + ' = ' + values[n])

    # Get network data from delta
    sensor_types, network_data_headers, network_data = get_network_data(values[0])

    # Prepare files for code output
    output_file_names = []
    outfile_file_headers = []
    for sensor_type in sensor_types:
        output_file_names.append(str.lower(sensor_type).replace(' ', '_') + '_sites.csv')
        if sensor_type not in ['GNSS', 'Camera']:
            outfile_file_headers.append(2)
        elif sensor_type == 'GNSS':
            outfile_file_headers.append(0)
        elif sensor_type == 'Camera':
            outfile_file_headers.append(1)

    # Write network data to files
    for m in range(len(output_file_names)):
        with open(output_file_names[m], 'w') as openfile:
            header = ''
            for col in network_data_headers[outfile_file_headers[m]][:-1]:
                header += col + ','
            header = header[:-1] + '\n'
            openfile.write(header)
        with open(output_file_names[m], 'a') as openfile:
            for n in range(len(network_data)):
                for k in range(len(network_data[n])):
                    if network_data[n][k][-1] == sensor_types[m]:
                        row = ''
                        if '9999-01-01T00:00:00Z' in network_data[n][k]:
                            for l in range(len(network_data[n][k]) - 1):
                                row += network_data[n][k][l] + ','
                            row = row[:-1] + '\n'
                            openfile.write(row)
