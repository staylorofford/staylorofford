"""
Before running this code:
- download .csv.xlsx files from google-drive
- using ssconvert, convert the file from xlsx to csv:
    for file in `ls *xlsx`; do ssconvert $file ${file:0:12}; rm $file; done
"""

import argparse
import datetime
import glob

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

    if not threshold:  # If no threshold is given, set to default 10 m
        threshold=10

    # Parse time series data in files

    data_columns = ['time', 'lat', 'lon', 'elevation', 'accuracy']
    northgoing_data_lists = [[], [], [], [], []]
    southgoing_data_lists = [[], [], [], [], []]
    data_lists = [northgoing_data_lists, southgoing_data_lists]

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

    # Parse path point locations from path file

    path_columns = ['lon', 'lat', 'elevation']
    northgoing_path_points = parse_path(northgoing_path_file)
    southgoing_path_points = parse_path(southgoing_path_file)
    point_lists = [northgoing_path_points, southgoing_path_points]

    # Generate distance lists for each path
