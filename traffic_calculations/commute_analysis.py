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

def parse_sorted_path(sorted_path_file):

    path_points = [[], [], []]
    distances = []

    return path_points, distances


def distance(p1, p2):

    """

    Calculate the distance between two points in 2D space, assuming no elevation difference exists between the points.

    :param p1: point 1 horizontal position [X, Y]
    :param p2: point 2 horizontal position [X, Y]
    :return: distance between p1 and p2
    """

    # Calculate distance

    d = math.sqrt((p2[0] - p1[0]) ** 2 + (p2[1] - p1[1]) ** 2)

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
    parser.add_argument('sorted_northgoing_path_csv', type=str, help="Path to csv file containing sorted"
                                                                     " northgoing path points")
    parser.add_argument('sorted_southgoing_path_csv', type=str, help="Path to csv file containing sorted"
                                                                     "southgoing path points")
    parser.add_argument('file_directory', type=str, help="Path to directory containing time series csv files")
    parser.add_argument('--threshold', type=str, help="Optional: location accuracy value above which to discard data " +
                                                      "points. Default = 10 m")
    args = parser.parse_args()

    northgoing_path_file = args.northgoing_path_csv
    southgoing_path_file = args.northgoing_path_csv
    file_dir = args.file_directory
    threshold = args.threshold

    # Parse commute point lists

    northgoing_point_list, northgoing_distance = parse_sorted_path(northgoing_path_file)
    southgoing_point_list, southgoing_distance = parse_sorted_path(southgoing_path_file)

    point_lists = [northgoing_point_list, southgoing_point_list]
    distance_lists = [northgoing_distance, southgoing_distance]

    # Parse time series data in files

    data_columns = ['time', 'X', 'Y', 'accuracy']
    dc_idx = [3, 0, 1, 7]  # indices of data columns in data files
    northgoing_data_lists = [[] for i in range(len(data_columns))]
    southgoing_data_lists = [[] for i in range(len(data_columns))]
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

                # Determine if the point's accuracy is sufficient

                if float(cols[dc_idx[data_columns.index('accuracy')]].replace("\"", "")) > threshold:
                    continue

                # Determine if data point is north- or south-going

                if cols[dc_idx[data_columns.index('time')]][:10].replace('-', '') == start_date:
                    m = 0
                else:
                    m = 1

                for n in range(len(data_columns)):
                    try:
                        data_lists[m][n].append(datetime.datetime.strptime(cols[dc_idx[n]][:19], '%Y-%m-%dT%H:%M:%S'))
                    except:
                        data_lists[m][n].append(float(cols[dc_idx[n]].replace("\"", "")))

    # Calculate distances of each data point along the path

    for m in range(len(data_lists)):  # NOTE: all data is under data_lists[1] due to some error above
        # print(m)
        for n in range(len(data_lists[m][0])):
            # print(n)
            distance_sums = [999999999]
            data_point = [data_lists[m][1][n], data_lists[m][2][n]]
            for k in range(1, len(point_lists[m][0])):  # Find the vertices that surround the data point
                vertex_points = [[point_lists[m][0][k - 1], point_lists[m][1][k - 1]],
                                 [point_lists[m][0][k], point_lists[m][1][k]]]
                distance_sum = 0
                print(vertex_points)
                for vertex_point in vertex_points:
                    distance_sum += distance(data_point, vertex_point)

                # print(distance_sum)
                if distance_sum > distance_sums[-1]:
                    # At the set of vertices just beyond those containing the data point
                    # print('BREAK')
                    break
                else:
                    distance_sums.append(distance_sum)

            # For loop breaks at vertices surrounding data point, translate vertices
            # and data point into vectors relative to the "behind" vertex

            direction_vector = [vertex_points[1][0] - vertex_points[0][0],
                                vertex_points[1][1] - vertex_points[0][1],
                                0]
            relative_data_point = [data_point[0] - vertex_points[0][0],
                                   data_point[1] - vertex_points[0][1],
                                   0]
            print(direction_vector, relative_data_point)
            try:
                distance_on_line = dot_product(direction_vector, relative_data_point) / normalise(direction_vector)[1]

                data_vertex_cross_product = cross_product(direction_vector, relative_data_point)
                distance_from_line = ((data_vertex_cross_product[2] / abs(data_vertex_cross_product[2])) * \
                                      normalise(data_vertex_cross_product)[1] / normalise(direction_vector)[1])

                # print(distance_on_line, distance_from_line)
            except:  # Fails if direction vector has length 0
                pass