#!/usr/bin/env python3

"""
Produce earthquake locations
"""

import argparse
import datetime
import math
import numpy as np


def calculate_tt(grid_point, site_location, velocity_model):

    """
    Calculate the travel time of a seismic wave from a given grid point to a given site for a given 1D velocity model.
    Assumes ray travels in the plane containing the grid point, site, and the vertical axis.
    :param grid_point: x,y,z position of a point in the grid
    :param site_location: x,y,z position of a site in the grid
    :param velocity_model: velocity model as parsed from velocity file
    :return: travel time between grid point and site location in seconds
    """

    # Find which layer the grid point sits in, this being the bottom layer to use in ray tracing
    for n in range(len(velocity_model)):
        if velocity_model[n][0] <= grid_point[2] <= velocity_model[n][1]:
            nmax = n

    # Perform ray tracing to find the travel time between the grid point and site location for the first arriving ray
    horizontal_proximities = []
    travel_times = []
    for psi in range(1, 180):

        # Initialise ray tracing parameters
        n = nmax
        travel_time = 0
        current_angle = psi

        # Collapse the 3D positions of the two points into the plane containing both points and the vertical axis
        # with origin as the site location.
        hpos = math.sqrt((site_location[0] - grid_point[0]) ** 2 +
                         (site_location[1] - grid_point[1]) ** 2)
        vpos = grid_point[2] - site_location[2]

        # Trace the ray up through velocity layers until it reaches the depth of the site
        while True:

            # Find if the site is in the velocity layer containing the grid point
            # If it is, set the vertical distance traveled in the layer as that between the grid point and the site,
            # otherwise, set the vertical distance traveled in the layer as that between the grid point and the
            # layer roof.
            if velocity_model[n][0] < site_location[2]:
                dz = vpos - site_location[2]
            else:
                dz = vpos - velocity_model[n][0]

            # Increase the travel time by the time taken for the ray to traverse the hypotenuse of the
            # triangle with acute angle psi and opposite dz

            travel_time += (dz / math.sin(math.radians(current_angle))) / velocity_model[n][2]

            # Calculate horizontal distance traversed by ray over dz
            dh = dz / math.tan(math.radians(current_angle))

            # Adjust ray head position due to distance traversed
            hpos -= dh
            vpos -= dz

            # Once the ray reaches the same depth as the site, save the horizontal position and travel time of the ray
            if vpos == site_location[2]:

                # Once the ray hits the depth of the site, the ray tracing ends. The horizontal distance between
                # the site location and the piercing point of the ray and the plane perpendicular to the vertical
                # axis at the depth of the site indicates how close the ray comes to the site. The ray with the
                # minimum horizontal distance will be the best solution to the ray tracing problem and its travel
                # time is taken as the travel time of a ray from the grid point to the site.
                horizontal_proximities.append(abs(hpos))
                travel_times.append(travel_time)
                break

            # If not, adjust current angle to that of the refracted ray in the new velocity layer
            else:
                current_angle = math.degrees(
                                             math.asin(
                                                       (velocity_model[n - 1][2] / velocity_model[n][2]) *
                                                       math.sin(math.radians(90 - current_angle))
                                                      )
                                             )
                n -= 1

    # The travel time between the grid point and the site is that of the ray which is closest to the site when it
    # reaches the depth of the site.

    travel_time = travel_times[horizontal_proximities.index(min(horizontal_proximities))]

    return travel_time


def generate_tt_grid(xmin, xmax, ymin, ymax, zmin, zmax, xstep, ystep, zstep, network_data, velocity_model):

    """
    Generate travel times from each grid cell to each site in the network for the given velocity model.
    All distances are relative to the first site given in the network file.
    :param xmin: western distance (km) to extend grid to (-ve, or 0)
    :param xmax: eastern distance (km) to extend grid to (+ve, or 0)
    :param ymin: southern distance (km) to extend grid to (-ve, or 0)
    :param ymax: northern distance (km) to extend grid to (+ve, or 0)
    :param zmin: height (km) to extend grid to (-ve, or 0)
    :param zmax: depth (km) to extend grid to (+ve, or 0)
    :param xstep: distance between grid points in x direction (km)
    :param ystep: distance between grid points in y direction (km)
    :param zstep: distance between grid points in z direction (km)
    :param network_data: network model as parsed from network file
    :param velocity_model: velocity model as parsed from velocity file
    :return: grid points and travel times in a nested list
    """

    # Define grid points along axes

    gridx = np.linspace(xmin, xmax, int(round((xmax - xmin) / xstep + 1)))
    gridy = np.linspace(ymin, ymax, int(round((ymax - ymin) / ystep + 1)))
    gridz = np.linspace(zmin, zmax, int(round((zmax - zmin) / zstep + 1)))

    # Define site positions in grid

    site_grid_positions = []
    for site in network_data:
        site_grid_positions.append([])
        site_grid_positions[-1].append(site[1] - network_data[0][1])
        site_grid_positions[-1].append(site[2] - network_data[0][2])
        site_grid_positions[-1].append(site[3] - network_data[0][3])

    # Define every point in the grid and generate travel times for each grid point to each station
    # and append these in order behind the corresponding grid point in the nested lists

    grid_points = [[[[] for z in gridz] for y in gridy] for x in gridx]
    for i in range(len(gridx)):
        for j in range(len(gridy)):
            for k in range(len(gridz)):
                grid_points[i][j][k] = [gridx[i], gridy[j], gridz[k]]
                for m in range(len(site_grid_positions)):
                    grid_points[i][j][k].append(calculate_tt(grid_points[i][j][k], site_grid_positions[m],
                                                             velocity_model))

    return grid_points


# Parse arguments

parser = argparse.ArgumentParser()
parser.add_argument('--arrival-time-file', type=str, help='Comma-separated file containing event arrival times '
                                                          'details. First row is header with columns: '
                                                          'Event #,site1_P,site1_S,site2_P,...,siteN_S. '
                                                          'All other rows contain arrival time data as: '
                                                          'event number (e.g. 1),'
                                                          'P wave arrival time at given site in form '
                                                          '"YYYY-MM-DD HH:MM:SS" (UTC),'
                                                          'S wave arrival time at given site in form '
                                                          '"YYYY-MM-DD HH:MM:SS" (UTC). '
                                                          'If no data exist then a date prior '
                                                          'to 1900 should be entered.')
# Note, here the arrival time information is limited by the formatting in Excel. In the finished version of the code
# all time input should be ISO8601 or nan if no data exists.
parser.add_argument('--network-file', type=str, help='Comma-separated file containing network details. '
                                                     'First row is header with columns: '
                                                     'site,easting,northing,depth. '
                                                     'All other rows contain site details as: '
                                                     'site code,'
                                                     'easting (km),'
                                                     'northing (km),'
                                                     'depth (km, +ve direction is down). '
                                                     'Coordinate system is not important, so long as it is cartesian.')
parser.add_argument('--velocity-model', type=str, help='Comma-separated file containing velocity model. '
                                                       'First row is header with columns: '
                                                       'upper_depth,lower_depth,velocity. '
                                                       'All other rows contain site details as: '
                                                       'upper depth (km, +ve direction is down),'
                                                       'lower depth (km, +ve direction is down),'
                                                       'velocity (km/s).')
parser.add_argument('--grid-parameters', type=str, help='Comma-separated file containing parameters for search grid. '
                                                  'First row is header with columns: '
                                                  'xmin,xmax,ymin,ymax,zmin,zmax,xstep,ystep,zstep. '
                                                  'Second row contains grid parameters corresponding to header '
                                                  'columns with units in km.')
parser.add_argument('--grid-file', type=str, help='File saved by this script after travel time grid generation. '
                                                  'Giving this argument will disable travel time grid generation and'
                                                  'the travel time grid given in the file will be used.')
parser.add_argument('--method', type=str, help='Earthquake location method to use, options are: grid_search,'
                                               'RVT_semblance')
args = parser.parse_args()

# Parse files

arrival_time_data = []
with open(args.arrival_time_file, 'r') as openfile:
    header = -1
    for row in openfile:
        if header == -1:
            header = 0
            arrival_time_data_header = row[:-1].split(',')
        else:
            cols = row[:-1].split(',')
            arrival_time_data.append([])
            arrival_time_data[-1].append(cols[0])
            for col in cols[1:]:
                arrival_time_data[-1].append(datetime.datetime.strptime(col, '%Y-%m-%d %H:%M:%S'))

network_data = []
with open(args.network_file, 'r') as openfile:
    header = -1
    for row in openfile:
        if header == -1:
            header = 0
        else:
            cols = row[:-1].split(',')
            network_data.append([])
            network_data[-1] = [cols[0], float(cols[1]), float(cols[2]), float(cols[3])]

velocity_model = []
with open(args.velocity_model, 'r') as openfile:
    header = -1
    for row in openfile:
        if header == -1:
            header = 0
        else:
            cols = row[:-1].split(',')
            velocity_model.append([])
            for col in cols:
                velocity_model[-1].append(float(col))

# If no grid file is given, calculate the travel time grid from grid parameters
if not args.grid_file:

    grid_parameters = []
    with open(args.grid_parameters, 'r') as openfile:
        header = -1
        for row in openfile:
            if header == -1:
                header = 0
            else:
                cols = row[:-1].split(',')
                for col in cols:
                    grid_parameters.append(float(col))
    xmin, xmax, ymin, ymax, zmin, zmax, xstep, ystep, zstep = grid_parameters

    # Build travel time grid

    grid_points = generate_tt_grid(xmin, xmax, ymin, ymax, zmin, zmax, xstep, ystep, zstep,
                                   network_data, velocity_model)

    # Save travel time grid to file

    with open('tt_grid.csv', 'w') as outfile:
        outstr = 'x,y,z,'
        for n in range(len(network_data)):
            outstr += 'tt_' + network_data[n][0] + ','
        outstr = outstr[:-1] + '\n'
        outfile.write(outstr)
    with open('tt_grid.csv', 'a') as outfile:
        for i in range(len(grid_points)):
            for j in range(len(grid_points[i])):
                for k in range(len(grid_points[i][j])):
                    outstr = ''
                    for m in range(len(grid_points[i][j][k])):
                        outstr += str(grid_points[i][j][k][m]) + ','
                    outstr = outstr[:-1] + '\n'
                    outfile.write(outstr)

# Otherwise, build the travel time grid from a previously saved grid
# Note: haven't tested this thoroughly
else:
    print(args.grid_file)
    with open(args.grid_file, 'r') as infile:
        header = -1
        gridx = []
        gridy = []
        gridz = []
        travel_times = []
        for row in infile:
            if header == -1:
                header = 0
            else:
                cols = row.split(',')
                travel_times.append([])
                gridx.append(float(cols[0]))
                gridy.append(float(cols[1]))
                gridz.append(float(cols[2]))
                for col in cols[3:]:
                    travel_times[-1].append(float(col))

    gridx = list(set(gridx))
    gridy = list(set(gridy))
    gridz = list(set(gridz))
    grid_points = [[[[] for z in gridz] for y in gridy] for x in gridx]
    for i in range(len(gridx)):
        for j in range(len(gridy)):
            for k in range(len(gridz)):
                grid_points[i][j][k] = [gridx[i], gridy[j], gridz[k]]
                grid_points[i][j][k].extend(travel_times[i + j + k])
