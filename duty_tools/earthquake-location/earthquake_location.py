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
    :return: p-wave and s-wave travel times between grid point and site location in seconds
    """

    # Find which layer the grid point sits in, this being the bottom layer to use in ray tracing
    for n in range(len(velocity_model)):
        if velocity_model[n][0] <= grid_point[2] <= velocity_model[n][1]:
            nmax = n

    # Perform ray tracing to find the travel time between the grid point and site location for the first arriving ray
    horizontal_proximities = []
    p_travel_times = []
    s_travel_times = []
    for psi in range(1, 180):

        # Initialise ray tracing parameters
        n = nmax
        p_travel_time = 0
        s_travel_time = 0
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

            p_travel_time += (dz / math.sin(math.radians(current_angle))) / velocity_model[n][2]
            s_travel_time += (dz / math.sin(math.radians(current_angle))) / velocity_model[n][3]

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
                p_travel_times.append(p_travel_time)
                s_travel_times.append(s_travel_time)
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

    p_travel_time = p_travel_times[horizontal_proximities.index(min(horizontal_proximities))]
    s_travel_time = s_travel_times[horizontal_proximities.index(min(horizontal_proximities))]

    return p_travel_time, s_travel_time


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
                    grid_points[i][j][k].extend(calculate_tt(grid_points[i][j][k], site_grid_positions[m],
                                                             velocity_model))

    return grid_points, grid_header


# Parse arguments

parser = argparse.ArgumentParser()
parser.add_argument('--arrival-time-file', type=str, help='Comma-separated file containing event arrival times '
                                                          'details. First row is header with columns: '
                                                          'site1_P,site1_S,site2_P,...,siteN_S. '
                                                          'All other rows contain arrival time data as: '
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
                                                       'upper_depth,lower_depth,p_velocity,s_velocity. '
                                                       'All other rows contain site details as: '
                                                       'upper depth (km, +ve direction is down),'
                                                       'lower depth (km, +ve direction is down),'
                                                       'p-wave velocity (km/s),'
                                                       's-wave velocity (km/s).')
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
            for col in cols[1:]:
                if int(col[:4]) < 1900:
                    arrival_time_data[-1].append(float('nan'))
                else:
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
    grid_points, grid_header = generate_tt_grid(xmin, xmax, ymin, ymax, zmin, zmax, xstep, ystep, zstep,
                                                network_data, velocity_model)

    # Save travel time grid to file
    with open('tt_grid.csv', 'w') as outfile:
        outstr = 'x,y,z,'
        for n in range(len(network_data)):
            outstr += 'ptt_' + network_data[n][0] + ',stt_' + network_data[n][0] + ','
        outstr = outstr[:-1] + '\n'
        outfile.write(outstr)
        grid_header = outstr[:-1].split(',')  # Define the grid header in case the generated grid is used in this run
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
    with open(args.grid_file, 'r') as infile:
        header = -1
        gridx = []
        gridy = []
        gridz = []
        travel_times = []
        for row in infile:
            if header == -1:
                grid_header = row[:-1].split(',')
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

# Create a translation chart to find travel time data for a given column in the arrival time data.

translation_chart = []
for n in range(len(arrival_time_data_header)):
    for m in range(len(grid_header)):
        if m < 3:
            continue
        if ((arrival_time_data_header[n].split('_')[0] == grid_header[m].split('_')[1]) and
                arrival_time_data_header[n].split('_')[1] == str.upper(grid_header[m][:1])):  # site and wave match
            translation_chart.append(m)

# Locate all events (only location method currently available is grid search)

event_solutions = [[], [], [], [], [], [], [], []]
for n in range(len(arrival_time_data)):
    origin_grid = [[[[] for z in range(len(grid_points[0][0]))]
                    for y in range(len(grid_points[0]))]
                   for x in range(len(grid_points))]

    # For each event, consider each grid cell as a possible location
    for i in range(len(grid_points)):
        for j in range(len(grid_points[i])):
            for k in range(len(grid_points[i][j])):

                # Calculate all origin times from the arrival time minus the travel time from the grid cell to each site
                origin_times = []
                max_weight = 0
                weight_sum = 0
                for m in range(len(arrival_time_data[n])):
                    if arrival_time_data[m][n] != float('nan'):
                        origin_times.append(arrival_time_data[m][n] - grid_points[i][j][k][translation_chart[n]])

                # Calculate mean origin time and RMS error for the grid cell
                mean_ot = np.mean(origin_times)
                print(mean_ot)
                rms = 0
                for origin_time in origin_times:
                    rms += (origin_time - mean_ot) ** 2
                rms = math.sqrt(rms)

                # Calculate weight, maximum weight, sum of weights for confidence interval determination
                try:
                    weight = 1 / (rms ** 2)
                except:
                    weight = math.inf  # Catch zero division in the case of 0 RMS error
                if weight > max_weight:
                    max_weight = weight
                weight_sum += weight

                # Save weight, RMS error and mean origin time
                origin_grid[i][j][k] = [weight, rms, mean_ot]

    # Using the weights defined, find the 95% confidence region
    for c in np.linspace(0, 1, 100):
        c_weight_sum = 0
        ijk = []
        for i in range(len(grid_points)):
            for j in range(len(grid_points[i])):
                for k in range(len(grid_points[i][j])):
                    if origin_grid[i][j][k][0] > c * max_weight:
                        ijk.append([i, j, k])
                        c_weight_sum += origin_grid[i][j][k][0]
        c_weight_sum /= weight_sum
        P = 1 - c_weight_sum
        if P > 0.95:
            break

    # Define centre and uncertainty of 95% confidence region: x,y,z and origin time
    x, y, z, ots = [], [], [], []
    if len(ijk) > 0:
        for indices in ijk:
            x.append(grid_points[indices[0]][indices[1]][indices[2]][0])
            y.append(grid_points[indices[0]][indices[1]][indices[2]][1])
            z.append(grid_points[indices[0]][indices[1]][indices[2]][2])
            ots.append(mean_ot[indices[0]][indices[1]][indices[2]][2])
    else:
        print('No data exists in the 95% confidence region for this event.')
    x_err = (max(x) - min(x)) / 2
    y_err = (max(y) - min(y)) / 2
    z_err = (max(z) - min(z)) / 2
    ots_err = (max(ots) - min(ots)) / 2
    mid_x = min(x) + x_err
    mid_y = min(y) + y_err
    mid_z = min(z) + z_err
    mid_ot = min(ots) + ots_err

    # Store solution and error
    event_solutions[0] = mid_x
    event_solutions[1] = mid_y
    event_solutions[2] = mid_z
    event_solutions[3] = mid_ot
    event_solutions[4] = x_err
    event_solutions[5] = y_err
    event_solutions[6] = z_err
    event_solutions[7] = ots_err



