#!/usr/bin/env python3

"""
Produce earthquake locations
"""

import argparse
import datetime
import math
import numpy as np
from obspy.clients.fdsn import Client as FDSN_Client
import pyproj

# Set up FDSN client
client = FDSN_Client('https://service.geonet.org.nz')

# Set up pyproj coordinate system objects
proj_wgs84_geog = pyproj.Proj(init="epsg:4326")
proj_wgs84_geod = pyproj.Proj(init="epsg:4978")


def convert_wgs84_geo_geod(lat, lon, height):

    """
    Convert a WGS84 geographic coordinate to a 3D WGS84 cartesian coordinate.
    :param lat: coordinate latitude in decimal degrees
    :param lon: coordinate longitude in decimal degrees
    :param height: coordinate height in metres (+ve direction is up)
    :return: x,y,z: the respective directional representations of the input coordinate
    """

    x, y, z = pyproj.transform(proj_wgs84_geog, proj_wgs84_geod, lon, lat, height)

    return x, y, z


def get_origins(eventIDs):

    """
    Get origin parameters for all earthquakes with eventID in eventIDs.
    :param eventIDs: list of earthquake eventID
    :return: nested list of origin parameters (latitude, longitude, depth, origin time) of earthquakes, and
             nested lists of arrival time data headers and data for those picks used in the origin calculation, and
             nested lists containing network data with sites being all those used in the event origins.
    """

    all_event_origins = []
    all_event_ATDHs = []
    all_event_ATDs = []
    for eventID in eventIDs:

        # Get event details
        event = client.get_events(eventid=eventID)[0]

        # Get origin and picks
        origin = event.origins[0]
        arrivals = origin.arrivals  # Picks used in origin
        picks = event.picks  # All picks for event

        # Extract all arrival time data from picks, but only for those used in the origin
        arrival_data = [[], [], []]
        for n in range(len(picks)):
            for m in range(len(arrivals)):
                if arrivals[m].time_weight and arrivals[m].time_weight > 0:
                    if picks[n].resource_id == arrivals[m].pick_id:
                        arrival_data[0].append(picks[n].waveform_id['station_code'])
                        arrival_data[1].append(arrivals[m].phase)
                        arrival_data[2].append(picks[n].time)

        # Build arrival time data header for the event
        event_ATDH = []
        sites = list(set(arrival_data[0]))
        for site in sites:
            event_ATDH.append(site + '_P')
            event_ATDH.append(site + '_S')

        # Organise arrival time data into the required format
        event_ATD = [0] * len(event_ATDH)
        for n in range(len(arrival_data[0])):
            site = arrival_data[0][n]
            phase = arrival_data[1][n]
            event_ATD[event_ATDH.index(site + '_' + phase)] = arrival_data[2][n]

        # Store origin and arrival time data
        all_event_origins.append([origin.latitude, origin.longitude, origin.depth, origin.time])
        all_event_ATDHs.append(event_ATDH)
        all_event_ATDs.append(event_ATD)

    # Build arrival time header and data for all events
    all_event_ATDHs_unordered = []
    for n in range(len(all_event_ATDHs)):
        all_event_ATDHs_unordered.extend(all_event_ATDHs[n])
    arrival_time_data_header = list(set(all_event_ATDHs_unordered))

    arrival_time_data = [[float('nan') for n in range(len(arrival_time_data_header))] for m in range(len(all_event_ATDs))]
    for m in range(len(all_event_ATDs)):
        for n in range(len(all_event_ATDs[m])):
            header = all_event_ATDHs[m][n]
            idx = arrival_time_data_header.index(header)
            try:
                arrival_time_data[m][idx] = all_event_ATDs[m][n].datetime
            except:
                arrival_time_data[m][idx] = float('nan')

    # Build network data using site list and FDSN
    network_data = []
    for n in range(len(arrival_time_data_header)):
        network_data.append([])
        # Get site name and details
        site = arrival_time_data_header[n].split('_')[0]
        site_details = client.get_stations(network='NZ',
                                           station=site)[0][0]

        # Convert WGS84 site geographic coordinates to WGS84 geodetic coordinates
        x, y, _ = convert_wgs84_geo_geod(site_details.latitude, site_details.longitude, site_details.elevation)
        network_data[-1] = [site, x, y, site_details.elevation]  # Keep z data relative to spheroid

    return all_event_origins, arrival_time_data_header, arrival_time_data, network_data


def calculate_tt(grid_point, site_location, velocity_model, phase):

    """
    Calculate the travel time of a seismic wave from a given grid point to a given site for a given 1D velocity model.
    Assumes ray travels in the plane containing the grid point, site, and the vertical axis.
    :param grid_point: x,y,z position of a point in the grid
    :param site_location: x,y,z position of a site in the grid
    :param velocity_model: velocity model as parsed from velocity file
    :param phase: which phase to calculate travel time for: P or S
    :return: travel time between grid point and site location in seconds for the given phase
    """

    # Find which layer the grid point sits in, this being the bottom layer to use in ray tracing
    for n in range(len(velocity_model)):
        try:
            if velocity_model[n][0] <= grid_point[2] <= velocity_model[n + 1][0]:
                nmax = n
        except:  # If the current layer is the bottom layer in the velocity model
            if velocity_model[n][0] <= grid_point[2]:
                nmax = n

    # Perform ray tracing to find the travel time between the grid point and site location for the first arriving ray
    horizontal_proximities = []
    travel_times = []
    for psi in range(1, 180):

        # Initialise ray tracing parameters
        n = nmax
        travel_time = 0
        current_angle = psi
        v_idx = ['P', 'S'].index(phase)

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

            travel_time += (dz / math.sin(math.radians(current_angle))) / velocity_model[n][v_idx + 1]

            # Calculate horizontal distance traversed by ray over dz
            dh = dz / math.tan(math.radians(current_angle))

            # Adjust ray head position due to distance traversed
            hpos -= dh
            vpos -= dz

            # Once the ray reaches the same depth as the site, save the horizontal position and travel time of the ray
            if vpos - site_location[2] <= 0.01:

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
                try:
                    current_angle = 90 - math.degrees(
                                                      math.asin(
                                                                (velocity_model[n - 1][v_idx + 1] /
                                                                 velocity_model[n][v_idx + 1]) *
                                                                math.sin(math.radians(90 - current_angle))
                                                                )
                                                      )
                    n -= 1

                except ValueError:
                    print('ValueError! Is the start depth of your upper velocity inclusive of all the relative '
                          'locations of sites in your network? If not, the ray tracing will fail.')
                    print(current_angle, n, vpos, site_location[2], velocity_model[n - 1][v_idx + 1],
                          velocity_model[n][v_idx + 1])
                    exit()

    # The travel time between the grid point and the site is that of the ray which is closest to the site when it
    # reaches the depth of the site.
    travel_time = travel_times[horizontal_proximities.index(min(horizontal_proximities))]

    return travel_time


def generate_tt_grid(network_data, velocity_model, test_origins = None,
                     xmin=-1, xmax=1, ymin=-1, ymax=1, zmin=-1, zmax=1, xstep=1, ystep=1, zstep=1):

    """
    Generate travel times from each grid cell to each site in the network for the given velocity model.
    All site locations are relative to the first site given in the network file.
    Default values are given for grid parameters in case test_origins argument is given.
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
    :param test_origins: test origins locations as parsed from the test origins file. Note, if this argument is
                         given then the travel time grid will only contain travel times from the test origins.
    :return: grid points and travel times in a nested list
    """

    # Define grid points along axes
    if not test_origins:
        gridx = np.linspace(xmin, xmax, int(round((xmax - xmin) / xstep + 1)))
        gridy = np.linspace(ymin, ymax, int(round((ymax - ymin) / ystep + 1)))
        gridz = np.linspace(zmin, zmax, int(round((zmax - zmin) / zstep + 1)))
    else:
        gridx, gridy, gridz = [], [], []
        for n in range(len(test_origins)):
            gridx.append(float(test_origins[n][0]))
            gridy.append(float(test_origins[n][1]))
            gridz.append(float(test_origins[n][2]))

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
                                                             velocity_model, phase='P'))
                    grid_points[i][j][k].append(calculate_tt(grid_points[i][j][k], site_grid_positions[m],
                                                             velocity_model, phase='S'))

    return grid_points


def grid_search(arrival_time_data, arrival_time_data_header, grid_points, grid_header):

    """
    Use grid search location algorithm to produce earthquake location and origin time using 95% confidence interval.
    :param arrival_time_data:
    :param arrival_time_data_header:
    :param grid_points:
    :param grid_header:
    :return: nested list of event solutions: [x, y, z, origin time, xerr, yerr, zerr, oterr] for each solution
    """

    # Create a translation chart to find travel time data for a given column in the arrival time data.

    translation_chart = []
    for m in range(len(arrival_time_data_header)):
        for n in range(len(grid_header)):
            if n < 3:
                continue
            if ((arrival_time_data_header[m].split('_')[0] == grid_header[n].split('_')[1]) and
                    arrival_time_data_header[m].split('_')[1] == str.upper(grid_header[n][:1])):  # site and wave match
                translation_chart.append(n)

    # Locate all events (only location method currently available is grid search)

    event_solutions = [[], [], [], [], [], [], [], [], []]
    for m in range(len(arrival_time_data)):
        origin_grid = [[[[] for z in range(len(grid_points[0][0]))]
                        for y in range(len(grid_points[0]))]
                       for x in range(len(grid_points))]

        # For each event, consider each grid cell as a possible location
        for i in range(len(grid_points)):
            for j in range(len(grid_points[i])):
                for k in range(len(grid_points[i][j])):

                    # Calculate all origin times from the arrival time minus the travel time from
                    # the grid cell to each site
                    origin_times = []
                    max_weight = 0
                    weight_sum = 0
                    for n in range(len(arrival_time_data[m])):
                        if isinstance(arrival_time_data[m][n], datetime.datetime):
                            origin_times.append(arrival_time_data[m][n] -
                                                datetime.timedelta(seconds=grid_points[i][j][k][translation_chart[n]])
                                                )

                    # Calculate mean origin time
                    ot_diffs = []
                    for origin_time in origin_times:
                        ot_diffs.append(origin_time - origin_times[0])
                    mean_ot = origin_times[0] + np.mean(ot_diffs)

                    # Calculate RMS error
                    rms = 0
                    for origin_time in origin_times:
                        rms += (origin_time - mean_ot).total_seconds() ** 2
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
        if len(grid_points) * len(grid_points[0]) * len(grid_points[0][0]) > 1:  # Catch test origin case
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
                    ots.append(origin_grid[indices[0]][indices[1]][indices[2]][2])
            else:
                print('No data exists in the 95% confidence region for this event.')
            x_err = (max(x) - min(x)) / 2
            y_err = (max(y) - min(y)) / 2
            z_err = (max(z) - min(z)) / 2
            ots_err = ((max(ots) - min(ots)) / 2).total_seconds()
            mid_x = min(x) + x_err
            mid_y = min(y) + y_err
            mid_z = min(z) + z_err
            mid_ot = min(ots) + ots_err

            # Calculate RMS error
            rms = 0
            for origin_time in ots:
                rms += (origin_time - mid_ot).total_seconds() ** 2
            rms = math.sqrt(rms)

        else:
            x_err = float('nan')
            y_err = float('nan')
            z_err = float('nan')
            ots_err = float('nan')
            mid_x = grid_points[0][0][0][0]
            mid_y = grid_points[0][0][0][1]
            mid_z = grid_points[0][0][0][2]
            mid_ot = origin_grid[0][0][0][2]
            rms = float('nan')

        # Store solution and error
        event_solutions[0].append(mid_x)
        event_solutions[1].append(mid_y)
        event_solutions[2].append(mid_z)
        event_solutions[3].append(mid_ot)
        event_solutions[4].append(x_err)
        event_solutions[5].append(y_err)
        event_solutions[6].append(z_err)
        event_solutions[7].append(ots_err)
        event_solutions[8].append(rms)

    return event_solutions


if __name__ == "__main__":

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
    parser.add_argument('--eventid-file', type=str, help='File containing earthquake eventIDs for earthquakes in the '
                                                         'GeoNet catalogue. Arrival time data is pulled from the GeoNet'
                                                         'catalogue for each earthquake in the file. '
                                                         'First row is header with column:'
                                                         'eventID. '
                                                         'All other rows contain earthquake eventIDs. '
                                                         'If this argument is given, network data will be '
                                                         'generated automatically and the --network-file argument will '
                                                         'be ignored')
    parser.add_argument('--test-origins', type=str, help='Comma-separated file containing earthquake origins to test '
                                                         'against arrival times. Rows must be aligned with '
                                                         'that of arrival time data in arrival time file or eventIDs '
                                                         'in eventid-file file. '
                                                         'First row is header with columns: '
                                                         'latitude,longitude,depth,origin_time.'
                                                         'All other rows contain earthquake origins details as:'
                                                         'latitude (decimal degrees), '
                                                         'longitude (decimal degrees), '
                                                         'depth (km, +ve direction is down),'
                                                         'origin time (IS08601 format).')
    parser.add_argument('--network-file', type=str, help='Comma-separated file containing network details. '
                                                         'First row is header with columns: '
                                                         'site,latitude,longitude,elevation. '
                                                         'All other rows contain site details as: '
                                                         'site code,'
                                                         'latitude (decimal degrees),'
                                                         'longitude (decimal degrees),'
                                                         'elevation (km, +ve direction is up). ')
    parser.add_argument('--velocity-model', type=str, help='Comma-separated file containing velocity model. '
                                                           'First row is header with columns: '
                                                           'start_depth,p_velocity,s_velocity. '
                                                           'All other rows contain site details as: '
                                                           'start depth (km, +ve direction is down),'
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

    # If a file of eventIDs was given
    if args.eventid_file:

        # Build the list of eventIDs
        eventIDs = []
        with open(args.eventid_file, 'r') as openfile:
            header = -1
            for row in openfile:
                if header == -1:
                    header = 0
                else:
                    eventIDs.append(row[:-1])

        # Build the arrival time and network data lists using FDSN queries and the eventIDs
        all_event_origins, arrival_time_data_header, arrival_time_data, network_data = get_origins(eventIDs)

    # Otherwise, generate the arrival time and network data lists from the respective files
    else:
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

                    # Convert WGS84 geographic coordinates to WGS84 geodetic coordinates
                    x, y, _ = convert_wgs84_geo_geod(float(cols[1]), float(cols[2]), float(cols[3]))
                    network_data[-1] = [cols[0], x, y, float(cols[3])]  # Keep z data relative to spheroid

    # If desired, parse the origins to test
    if args.test_origins:
        test_origins = []
        with open(args.test_origins, 'r') as openfile:
            header = -1
            for row in openfile:
                if header == -1:
                    header = 0
                else:
                    cols = row[:-1].split(',')
                    test_origins.append(cols)

    # Build the velocity model lists from the velocity model file
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
    if not args.grid_file and not args.test_origins:

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
        grid_points = generate_tt_grid(network_data,
                                       velocity_model,
                                       xmin=xmin,
                                       xmax=xmax,
                                       ymin=ymin,
                                       ymax=ymax,
                                       zmin=zmin,
                                       zmax=zmax,
                                       xstep=xstep,
                                       ystep=ystep,
                                       zstep=zstep)

        # Save travel time grid to file
        with open('tt_grid.csv', 'w') as outfile:
            outstr = 'x,y,z,'
            for n in range(len(network_data)):
                outstr += 'ptt_' + network_data[n][0] + ',stt_' + network_data[n][0] + ','
            outstr = outstr[:-1] + '\n'
            outfile.write(outstr)

            # Define the grid header in case the generated grid is used in this run
            grid_header = outstr[:-1].split(',')

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
    elif not args.test_origins:
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

    # Otherwise, build the travel time grid and header for the test origins
    else:
        grid_points = generate_tt_grid(network_data,
                                       velocity_model,
                                       test_origins=test_origins)

        grid_header = 'x,y,z,'
        for n in range(len(network_data)):
            grid_header += 'ptt_' + network_data[n][0] + ',stt_' + network_data[n][0] + ','
        grid_header = grid_header[:-1].split(',')

    # Perform earthquake location
    method = args.method
    if method == 'grid_search':
        if args.test_origins:
            for n in range(len(test_origins)):
                earthquake_origins = grid_search([arrival_time_data[n]],
                                                 arrival_time_data_header,
                                                 [[[grid_points[n][n][n]]]],
                                                 grid_header)

                # Calculate RMS error as the time difference between the test origin time and that from the grid search
                rms = (datetime.datetime.strptime(test_origins[n][3],
                                                  '%Y-%m-%dT%H:%M:%SZ') - earthquake_origins[3][0]).total_seconds()
                print(earthquake_origins, rms)
