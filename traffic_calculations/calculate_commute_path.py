import argparse
import datetime
import glob
from itertools import permutations
import math
import matplotlib.pyplot as plt

def parse_path(path_file, path_columns):

    """
    :param path_file: path to path file
    :return: list of line segment vertices on path
    """

    vertex_list = []
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

            try:  # Check vertex is useful data
                float(cols[len(path_columns) - 1])
            except:  # Line identificator missing
                continue

            vertex_list.append([])

            for n in range(len(path_columns)):
                vertex_list[-1].append(float(cols[n]))

    return vertex_list


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


if __name__ == "__main__":

    # Parse arguments

    parser = argparse.ArgumentParser(description='What does this script do?')
    parser.add_argument('northgoing_path_csv', type=str, help="Path to csv file containing northgoing path points")
    parser.add_argument('southgoing_path_csv', type=str, help="Path to csv file containing southgoing path points")
    args = parser.parse_args()

    northgoing_path_file = args.northgoing_path_csv
    southgoing_path_file = args.southgoing_path_csv
    path_files = [northgoing_path_file, southgoing_path_file]

    # Parse path point locations from path file

    print('\nParsing path vertex data...\n')

    path_columns = ['X', 'Y', 'Z', 'fid']
    all_vertex_lists = [[], []]
    for m in range(len(path_files)):
        all_vertex_lists[m] = parse_path(path_files[m], path_columns)

    # Remove any duplicates from the data (these do exist, for whatever reason!)

    vertex_lists = [[], []]
    for m in range(len(all_vertex_lists)):
        for n in range(len(all_vertex_lists[m])):
            if all_vertex_lists[m][n] in vertex_lists[m]:
                continue
            else:
                vertex_lists[m].append(all_vertex_lists[m][n])

    # Order vertices into a commute path

    for m in range(len(vertex_lists)):

        end_points = []
        end_point_indices = []
        # Find line segment end points for the path
        for n in range(len(vertex_lists[m])):
            for k in range(len(vertex_lists[m])):
                if n != k and distance(vertex_lists[m][n][:-1], vertex_lists[m][k][:-1]) == 0:
                    end_points.append(vertex_lists[m][n][:-1])
                    end_point_indices.append(n)

        # Find along-line direction of indices for first line segment on path

        initial_idx = end_point_indices[0]
        path_points = [vertex_lists[m][initial_idx]]
        if vertex_lists[m][initial_idx + 1][3] == vertex_lists[m][initial_idx][3]:
            sign = 1
        elif vertex_lists[m][initial_idx - 1][3] == vertex_lists[m][initial_idx][3]:
            sign = -1

        # Build the path in the direction of the first line segment

        idx = initial_idx
        idx_visited = [idx]
        while True:
            idx += sign
            if vertex_lists[m][idx][:-1] in end_points:

                # Find the index of the associated end point

                idx_indices = []
                for i in range(len(end_points)):
                    if end_points[i] == vertex_lists[m][idx][:-1]:
                        idx_indices.append(i)
                idx = end_point_indices[idx_indices[idx_indices.index(end_point_indices.index(idx)) - 1]]

                # Find the counting direction on the associated line,
                # unless it is a endpoint of the data, or somewhere already visited,
                # in the latter case exist, in the former case simply do not add the
                # point to the path (this occurs due to overlapping line segments that
                # are from separate lines).

                if idx not in idx_visited:
                    duplicate = False
                    for path_point in path_points:
                        if vertex_lists[m][idx][:-1] == path_point[:-1]:
                            duplicate = True
                    if not duplicate:
                        path_points.append(vertex_lists[m][idx])
                    idx_visited.append(idx)
                else:
                    break

                if idx != 0 and idx != (len(vertex_lists[m]) - 1):
                    if vertex_lists[m][idx + 1][3] == vertex_lists[m][idx][3]:
                        sign = 1
                    elif vertex_lists[m][idx - 1][3] == vertex_lists[m][idx][3]:
                        sign = -1
                else:
                    break
            else:

                # Add the point to the path, unless it is somewhere already visited

                if idx not in idx_visited:
                    duplicate = False
                    for path_point in path_points:
                        if vertex_lists[m][idx][:-1] == path_point[:-1]:
                            duplicate = True
                    if not duplicate:
                        path_points.append(vertex_lists[m][idx])
                    idx_visited.append(idx)
                else:
                    break

        # Reverse the path so far

        path_points = list(reversed(path_points))

        # Build the path in the opposite direction

        idx = initial_idx
        idx_visited = [idx]
        while True:
            if vertex_lists[m][idx][:-1] in end_points:

                # Find the index of the associated end point

                idx_indices = []
                for i in range(len(end_points)):
                    if end_points[i] == vertex_lists[m][idx][:-1]:
                        idx_indices.append(i)
                idx = end_point_indices[idx_indices[idx_indices.index(end_point_indices.index(idx)) - 1]]

                # Find the counting direction on the associated line,
                # unless it is a endpoint of the data, or somewhere already visited,
                # in the latter case exist, in the former case simply do not add the
                # point to the path (this occurs due to overlapping line segments that
                # are from separate lines).

                if idx not in idx_visited:
                    duplicate = False
                    for path_point in path_points:
                        if vertex_lists[m][idx][:-1] == path_point[:-1]:
                            duplicate = True
                    if not duplicate:
                        path_points.append(vertex_lists[m][idx])
                    idx_visited.append(idx)
                else:
                    break

                if idx != 0 and idx != (len(vertex_lists[m]) - 1):
                    if vertex_lists[m][idx + 1][3] == vertex_lists[m][idx][3]:
                        sign = 1
                    elif vertex_lists[m][idx - 1][3] == vertex_lists[m][idx][3]:
                        sign = -1
                else:
                    break
            else:

                # Add the point to the path, unless it is somewhere already visited

                if idx not in idx_visited:
                    duplicate = False
                    for path_point in path_points:
                        if vertex_lists[m][idx][:-1] == path_point[:-1]:
                            duplicate = True
                    if not duplicate:
                        path_points.append(vertex_lists[m][idx])
                    idx_visited.append(idx)
                else:
                    break
            idx += sign

        if m == len(vertex_lists) - 1:
            # For the return direction, have the first point at work
            path_points = list(reversed(path_points))

        # Calculate the distance of each point along the path

        line_lengths = [0]
        distances = [0]
        for n in range(1, len(path_points)):
            line_lengths.append(distance(path_points[n - 1], path_points[n]))
            distances.append(int(sum(line_lengths)))

        # Write the path and distances to file

        with open('sorted_' + path_files[m], 'w') as outfile:
            outfile.write('X,Y,Z\n')
        with open('sorted_' + path_files[m], 'a') as outfile:
            for n in range(len(path_points)):
                outfile.write(str(path_points[n][0]) + ',' +
                              str(path_points[n][1]) + ',' +
                              str(path_points[n][2]) + ',' +
                              str(distances[n]) + '\n')
