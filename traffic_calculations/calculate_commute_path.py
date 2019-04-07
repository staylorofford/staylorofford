import argparse
import datetime
import glob
from itertools import permutations
import math
import matplotlib.pyplot as plt

def parse_path(path_file, path_columns):

    """
    :param path_file: path to path file
    :return: nested lists of path point coordinates, list of path point vertices
    """

    path_points = [[], [], []]
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
            vertex_list.append([])
            for n in range(len(path_columns)):
                path_points[n].append(float(cols[n]))
                vertex_list[-1].append(float(cols[n]))

    return path_points, vertex_list


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

    path_columns = ['X', 'Y', 'Z']
    vertex_lists = [[], []]
    for m in range(len(path_files)):
        vertex_lists[m] = parse_path(path_files[m], path_columns)

    # Calculate path through vertices such that each vertex is on the path and the path has the
    # minimum possible distance

    print('\nSorting path vertices into path order...\n')

    for m in range(len(vertex_lists)):
        vertex_list_permutations = permutations(vertex_lists[m])
        path_lengths = []
        path_length_sums = []
        n = 0
        for permutation in vertex_list_permutations:
            print('Calculating path length of permutation ' + str(n))
            n += 1
            # Calculate distance along path
            path_length = []
            for k in range(1, len(permutation)):
                path_length.append(distance(permutation[k - 1], permutation[k]))
            print('Path length is ' + str(sum(path_length)))
            path_lengths.append(path_length)
            path_length_sums.append(sum(path_length))
        print(path_lengths)
        print(min(path_length_sums))

        # Build sorted point list from the permutation with minimum path length

        print('\nBuilding vertex coordinate list for path...\n')

        idx = path_length_sums.index(min(path_length_sums))
        sorted_point_list = [[], [], []]
        for n in range(len(vertex_list_permutations[idx])):
            for k in range(len(vertex_list_permutations[idx][n])):
                sorted_point_list[k].append(vertex_list_permutations[idx][n][k])

        print('\nPlotting path vertices...\n')

        # Show the operator the sorted vertices

        plt.scatter(sorted_point_list[0], sorted_point_list[1])
        for k in range(len(sorted_point_list[0])):
            plt.text(sorted_point_list[0][k], sorted_point_list[1][k], str(k))
        plt.show()

        # Write the path and distances to file

        with open('sorted_' + path_files[m], 'w') as outfile:
            outfile.write('X,Y,Z,distance\n')
        with open('sorted_' + path_files[m], 'a') as outfile:
            for k in range(len(sorted_point_list[0])):
                outfile.write(str(sorted_point_list[0][k]) + ',' +
                              str(sorted_point_list[1][k]) + ',' +
                              str(sorted_point_list[2][k]) + ',' +
                              str(path_lengths[idx]) + '\n')

