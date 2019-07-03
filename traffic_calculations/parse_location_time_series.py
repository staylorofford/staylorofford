"""
Before running this code:
- download .csv.xlsx files from google-drive
- using ssconvert, convert the file from xlsx to csv:
    for file in `ls *xlsx`; do ssconvert $file ${file:0:12}; rm $file; done
"""

import argparse
import datetime
import glob
import math

def parse_path(path_file):

    """
    :param path_file: path to path file
    :return: list of path point locations
    """

    path_points = [[], [], []]
    with open(path_file, 'r') as openfile:
        header = '0'
        for row in openfile:

            # Ignore csv header

            try:
                header += 1
            except:
                header = int(header)
                continue

            # Extract data into lists

            cols = row.split(',')
            for n in range(len(path_columns)):
                path_points[n].append(float(cols[n]))

    return path_points


def distance(p1, p2):

    """

    Calculate the distance between two points on a sphere, assuming no elevation difference exists between the points.

    :param p1: point 1 horizontal position [lon, lat]
    :param p2: point 2 horizontal position [lon, lat]
    :return: distance between p1 and p2
    """

    # Convert longitude and latitude to radians

    for n in range(len(p1)):
        p1[n] = math.radians(p1[n])
        p2[n] = math.radians(p2[n])

    # Calculate distance

    r = 6371000
    d = 2 * r * math.asin(math.sqrt((math.sin((p2[1] - p1[1]) / 2.0)) ** 2 + math.cos(p1[0]) * math.cos(p2[0]) *
                          (math.sin((p2[0] - p1[0]) / 2.0) ** 2)))
    return d


def cross_product(a, b):

    """
    Do cross product between a and b.
    """

    c = [a[1] * b[2] - a[2] * b[1],
         a[2] * b[0] - a[0] * b[2],
         a[0] * b[1] - a[1] * b[0]]

    return c


def dot_product(a, b):

    """
    Do dot product between a and b.
    """

    c = a[0] * b[0] + a[1] * b[1] + a[2] * b[2]

    return c


def normalise(a):

    """
    Return the normalised version of a and the magnitude of a.
    """

    a_mag = 0
    a_norm = []

    for value in a:
        a_mag += value ** 2
    a_mag = math.sqrt(a_mag)

    for value in a:
        a_norm.append(value / a_mag)

    return a_norm, a_mag



if __name__ == "__main__":

    # Parse arguments

    parser = argparse.ArgumentParser(description='What does this script do?')
    parser.add_argument('file_directory', type=str, help="Path to directory containing time series csv files")
    parser.add_argument('northgoing_path_csv', type=str, help="Path to csv file containing northgoing path points")
    parser.add_argument('southgoing_path_csv', type=str, help="Path to csv file containing southgoing path points")
    parser.add_argument('--threshold', type=str, help="Optional: location accuracy value above which to discard data " +
                                                    "points. Default = 10 m")
    args = parser.parse_args()

    file_dir = args.file_directory
    northgoing_path_file = args.northgoing_path_csv
    southgoing_path_file = args.southgoing_path_csv
    threshold = args.threshold

    # Parse path point locations from path file

    path_columns = ['lon', 'lat', 'elevation']
    northgoing_path_points = parse_path(northgoing_path_file)
    southgoing_path_points = parse_path(southgoing_path_file)
    point_lists = [northgoing_path_points, southgoing_path_points]

    # Calculate distance lists for each path

    distance_lists = []
    for point_list in point_lists:
        incremental_distance_sums = []
        incremental_distance_lists = []
        for n in range(len(point_list[0])):
            incremental_distances = [0]

            c = 0
            m = n
            while m >= 1:
                incremental_distance = -1 * distance([point_list[0][m], point_list[1][m]],
                                                 [point_list[0][m - 1], point_list[1][m - 1]])

                if int(abs(incremental_distance)) > 5:  # Ignore points out of path order
                    continue

                c += 1
                m = n - c

            c = 0
            m = n
            while m < len(point_list[0]) - 1:
                incremental_distance = distance([point_list[0][m], point_list[1][m]],
                                                     [point_list[0][m + 1], point_list[1][m + 1]])

                if int(abs(incremental_distance)) > 5:  # Ignore points out of path order
                    continue

                c += 1
                m = n + c

            incremental_distance_sums.append(sum(incremental_distances))
            incremental_distance_lists.append(incremental_distances)

        distance_lists.append(incremental_distance_lists[incremental_distance_sums.
                              index(max(incremental_distance_sums))])

    # Parse time series data in files

    data_columns = ['time', 'lat', 'lon', 'elevation', 'accuracy']
    northgoing_data_lists = [[], [], [], [], []]
    southgoing_data_lists = [[], [], [], [], []]
    data_lists = [northgoing_data_lists, southgoing_data_lists]

    if not threshold:  # If no threshold is given, set to default 10 m
        threshold = 10

    csv_files = glob.glob(file_dir + '*csv')
    for csv_file in csv_files:
        start_date = csv_file[:csv_file.index('.')]
        with open(csv_file, 'r') as openfile:
            header = '0'
            for row in openfile:

                # Ignore csv header

                try:
                    header += 1
                except:
                    header = int(header)
                    continue

                # Extract data into data_lists

                cols = row.split(',')

                if float(cols[4]) > threshold:  # Determine if the point's accuracy is sufficient
                    continue

                if cols[0][:10].replace('-', '') == start_date:  # Determine if data point is north- or south-going
                    m = 0
                else:
                    m = 1

                for n in range(len(data_columns)):
                    try:
                        data_lists[m][n].append(datetime.datetime.strptime(cols[n][:19], '%Y-%m-%dT%H:%M:%S'))
                    except:
                        data_lists[m][n].append(float(cols[n]))

    # Calculate distances of each data point along the path

    for m in range(len(data_lists)):
        for n in range(len(data_lists[m][0])):
            data_point = [data_lists[m][2][n], data_lists[m][1][n]]
            for k in range(len(point_lists[m][0])):  # Find the vertex which sits behind the data point
                vertex_point = [point_lists[m][0][k], point_lists[m][1][k]]
                if distance(vertex_point, data_point) < 0:
                    break
                # Transform vertex, next vertex, and data point into cartesian coordinates so dot and cross
                # products can be used to calculate distances along and off the line between vertices

                # something... something...

                # on_line_distance = dot_product(direction_vector, data_point) / direction_vector_magnitude
                # event_line_cross_product = cross_product(direction_vector, data_point)
                # off_line_distance = (event_line_cross_product[2] / abs(event_line_cross_product[2])) * \
                #                     normalise(event_line_cross_product)[1] / direction_vector_magnitude