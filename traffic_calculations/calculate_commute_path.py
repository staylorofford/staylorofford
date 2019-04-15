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

    point_lists = [[] for n in range(len(path_columns))]
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
                point_lists[n].append(float(cols[n]))
                vertex_list[-1].append(float(cols[n]))

    return point_lists, vertex_list


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
    points_lists = [[], []]
    vertex_lists = [[], []]
    for m in range(len(path_files)):
        points_lists[m], vertex_lists[m] = parse_path(path_files[m], path_columns)

    # Order vertices into a commute path

    end_points = []
    end_points_all = []  # used to ensure only 1 end point exists with a given location for each line
    end_point_indices = []
    for m in range(len(vertex_lists)):

        # Find line segment end points for the path
        for n in range(len(vertex_lists[m])):
            for k in range(len(vertex_lists[m])):
                if n != k and vertex_lists[m][n] not in end_points_all and distance(vertex_lists[m][n][:-1],
                                                                                    vertex_lists[m][k][:-1]) == 0:
                    end_points.append(vertex_lists[m][n][:-1])
                    end_points_all.append(vertex_lists[m][n])
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

        c = 0
        idx_visited = [idx]

        while True:
            idx += sign
            if vertex_lists[m][idx][:-1] in end_points:

                # Find the index of the associated end point

                print('')
                print(idx)
                # print(vertex_lists[m][idx])
                # idx_indices = [i for i, x in enumerate(end_points) if x == vertex_lists[m][idx][:-1]]
                idx_indices = []
                for i in range(len(end_points)):
                    # print(i)
                    if end_points[i] == vertex_lists[m][idx][:-1]:
                        print(i, idx)
                        print(end_points[i], vertex_lists[m][idx])
                        idx_indices.append(i)
                print(idx_indices)
                print(end_point_indices)
                print(end_point_indices.index(idx))
                print(idx_indices.index(end_point_indices.index(idx)))
                print(idx_indices.index(end_point_indices.index(idx)) - 1)
                print(idx_indices[idx_indices.index(end_point_indices.index(idx)) - 1])
                print(end_point_indices[idx_indices[idx_indices.index(end_point_indices.index(idx)) - 1]])
                idx = end_point_indices[idx_indices[idx_indices.index(end_point_indices.index(idx)) - 1]]
                # ^ this part works OK!

                # Find the counting direction on the associated line,
                # unless it is a path endpoint, in which case exit

                if idx in idx_visited:
                    print('BEEN HERE!')
                    print(c, idx)
                    for i in idx_indices:
                        print(vertex_lists[m][idx][:-1])
                        print(end_points[i])
                        print(vertex_lists[m][idx][:-1] == end_points[i])
                idx_visited.append(idx)

                path_points.append(vertex_lists[m][idx])

                for n in range(len(path_points)):
                    plt.scatter(path_points[n][0], path_points[n][1], color='k')
                    plt.text(path_points[n][0], path_points[n][1], str(n), color='red')
                plt.show()

                # Looks like line counter has a problem if the line consists of a bunch of stacked points!

                if idx != 0 and idx != (len(vertex_lists[m]) - 1):

                    if vertex_lists[m][idx + 1][3] == vertex_lists[m][idx][3]:
                        sign = 1
                    elif vertex_lists[m][idx - 1][3] == vertex_lists[m][idx][3]:
                        sign = -1
                else:
                    break
            else:

                if idx in idx_visited:
                    print('BEEN HERE!')
                    print(c, idx)
                idx_visited.append(idx)

                path_points.append(vertex_lists[m][idx])
            c += 1


        # Build the path in the opposite direction

        # print('')

        idx = initial_idx
        # print(idx)
        end_points = all_end_points
        end_point_indices = all_end_point_indices
        while idx != 0 and idx != (len(vertex_lists[m]) - 1):
            if vertex_lists[m][idx][:-1] in end_points:
                # Jump to the end point of the new line
                # ISSUE: poping seems to be removing required values to have this indexing method work!
                print(vertex_lists[m][idx][:-1])
                print('index is ' + str(idx))

                # Remove the end point from the reference lists

                end_points.pop(end_point_indices.index(idx))
                end_point_indices.pop(end_point_indices.index(idx))

                # Find the index of the associated end point in the data

                print('end points index is ' + str(end_points.index(vertex_lists[m][idx][:-1])))
                idx = end_point_indices[end_points.index(vertex_lists[m][idx][:-1])]
                print('associated index is ' + str(idx))
                print('data length is ' + str(len(vertex_lists[m])))

                # Remove the associated end point from the reference lists

                end_points.pop(end_point_indices.index(idx))
                end_point_indices.pop(end_point_indices.index(idx))

                # Find the counting direction on the associated line,
                # unless it is a path endpoint, in which case exit

                path_points.append(vertex_lists[m][idx])
                # if idx != 0 and idx != (len(vertex_lists[m]) - 1):

                if vertex_lists[m][idx + 1][3] == vertex_lists[m][idx][3]:
                    sign = 1
                elif vertex_lists[m][idx - 1][3] == vertex_lists[m][idx][3]:
                    sign = -1
                # else:
                #     break
            else:
                path_points.append(vertex_lists[m][idx])
            idx += sign

        print('DONE')
        print(len(path_points))
        print(len(vertex_lists[m]))
        for n in range(len(path_points)):
            plt.scatter(path_points[n][0], path_points[n][1], color='k')
            plt.text(path_points[n][0], path_points[n][1], str(n))
        plt.show()



            # Find which way the line approaches the end point




    # Calculate path through vertices such that each vertex is on the path and the path has the
    # minimum possible distance
    #
    # print('\nSorting path vertices into path order...\n')
    #
    # for m in range(len(vertex_lists)):
    #     vertex_list_permutations = permutations(vertex_lists[m])
    #     path_lengths = []
    #     path_length_sums = []
    #     n = 0
    #     for permutation in vertex_list_permutations:
    #         print('Calculating path length of permutation ' + str(n))
    #         n += 1
    #         # Calculate distance along path
    #         path_length = []
    #         for k in range(1, len(permutation)):
    #             path_length.append(distance(permutation[k - 1], permutation[k]))
    #         print('Path length is ' + str(sum(path_length)))
    #         path_lengths.append(path_length)
    #         path_length_sums.append(sum(path_length))
    #     print(path_lengths)
    #     print(min(path_length_sums))
    #     idx = path_length_sums.index(min(path_length_sums))
    #
    #     # Build sorted point list from the permutation with minimum path length
    #
    #     print('\nBuilding vertex coordinate list for path...\n')
    #
    #     sorted_point_list = [[], [], []]
    #     for n in range(len(vertex_list_permutations[idx])):
    #         for k in range(len(vertex_list_permutations[idx][n])):
    #             sorted_point_list[k].append(vertex_list_permutations[idx][n][k])
    #
    #     print('\nPlotting path vertices...\n')
    #
    #     # Show the operator the sorted vertices
    #
    #     plt.scatter(sorted_point_list[0], sorted_point_list[1])
    #     for k in range(len(sorted_point_list[0])):
    #         plt.text(sorted_point_list[0][k], sorted_point_list[1][k], str(k))
    #     plt.show()
    #
    #     # Write the path and distances to file
    #
    #     with open('sorted_' + path_files[m], 'w') as outfile:
    #         outfile.write('X,Y,Z,distance\n')
    #     with open('sorted_' + path_files[m], 'a') as outfile:
    #         for k in range(len(sorted_point_list[0])):
    #             outfile.write(str(sorted_point_list[0][k]) + ',' +
    #                           str(sorted_point_list[1][k]) + ',' +
    #                           str(sorted_point_list[2][k]) + ',' +
    #                           str(path_lengths[idx]) + '\n')