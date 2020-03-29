#!/usr/bin/env python3

"""
Produce earthquake locations. Only works for manually defined events, or those from the GeoNet catalogue.
"""

import argparse
import datetime
import io
import math
import numpy as np
from obspy.clients.fdsn import Client as FDSN_Client
from obspy.taup import TauPyModel
from obspy.io.quakeml.core import Unpickler
import pycurl
import pyproj
import xml.etree.ElementTree as ET

# Set up objects to use imported modules
client = FDSN_Client("https://service.geonet.org.nz")
spherical_velocity_model = TauPyModel(model="iasp91")
quakeml_reader = Unpickler()

# Set up pyproj coordinate system objects
proj_wgs84_geog = pyproj.Proj(init="epsg:4326")
proj_wgs84_geod = pyproj.Proj(init="epsg:4978")


def parse_files(arrival_time_file=None,
                eventid_file=None,
                test_origins=None,
                network_file=None,
                velocity_model=None,
                grid_parameters=None,
                grid_file=None,
                mode=None):

    """
    Parse parameters from files as described in the main execution of this code.
    """

    # If a file of eventIDs was given
    if eventid_file:

        # Build the list of eventIDs
        eventIDs = []
        with open(eventid_file, 'r') as openfile:
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
        with open(arrival_time_file, 'r') as openfile:
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
        with open(network_file, 'r') as openfile:
            header = -1
            for row in openfile:
                if header == -1:
                    header = 0
                else:
                    cols = row[:-1].split(',')
                    network_data.append([])

                    # Convert WGS84 geographic coordinates to WGS84 geodetic coordinates
                    x, y, _ = convert_wgs84_geo_geod(float(cols[1]), float(cols[2]), float(cols[3]))
                    network_data[-1] = [cols[0], x / 1000, y / 1000, -1 * float(cols[3]) / 1000,
                                        float(cols[1]), float(cols[2])]

    # If desired, parse the origins to test
    if test_origins:
        with open(test_origins, 'r') as openfile:
            test_origins = []
            header = -1
            for row in openfile:
                if header == -1:
                    header = 0
                else:
                    cols = row[:-1].split(',')
                    test_origins.append([])

                    # Convert WGS84 geographic coordinates to WGS84 geodetic coordinates
                    x, y, _ = convert_wgs84_geo_geod(float(cols[0]), float(cols[1]), -1 * float(cols[2]))
                    test_origins[-1] = [cols[3], x / 1000, y / 1000, float(cols[2]) / 1000,
                                        float(cols[0]), float(cols[1])]

    # If eventIDs are input, use the site positions from FDSN over any others.
    if not eventid_file:
        for n in range(len(all_event_origins)):
            x, y, _ = convert_wgs84_geo_geod(all_event_origins[n][0], all_event_origins[n][1],
                                                                 -1 * all_event_origins[n][2])
            all_event_origins[n] = [str(all_event_origins[n][3]), x / 1000, y / 1000, all_event_origins[n][2] / 1000,
                                    all_event_origins[n][0], all_event_origins[n][1]]
        test_origins = all_event_origins

    # Build the velocity model lists from the velocity model file
    if mode != 'spherical':
        with open(velocity_model, 'r') as openfile:
            velocity_model = []
            header = -1
            for row in openfile:
                if header == -1:
                    header = 0
                else:
                    cols = row[:-1].split(',')
                    velocity_model.append([])
                    for col in cols:
                        velocity_model[-1].append(float(col))
    else:
        velocity_model = None  # A velocity model is not required for spherical travel time calculation,
        # as this uses IASPEI 91.

    # If no grid file is given, calculate the travel time grid from grid parameters
    if not grid_file and not test_origins:
        with open(grid_parameters, 'r') as openfile:
            grid_parameters = []
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
                                       mode,
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
    elif not test_origins:
        with open(grid_file, 'r') as infile:
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
                                       mode,
                                       test_origins=test_origins)

        grid_header = 'x,y,z,'
        for n in range(len(network_data)):
            grid_header += 'ptt_' + network_data[n][0] + ',stt_' + network_data[n][0] + ','
        grid_header = grid_header[:-1].split(',')

    return arrival_time_data, arrival_time_data_header, grid_points, grid_header, test_origins


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


def FDSN_event_query(eventID):

    """
    Use obspy with pycurl to query event details via FDSN..

    :param eventID: FDSN earthquake eventID
    :return: obspy object containing earthquake details
    """

    query = 'https://service.geonet.org.nz/fdsnws/event/1/query?eventid=' + eventID
    queryresult = curl(query)
    event_details = quakeml_reader.loads(queryresult)

    return event_details


def FDSN_station_query(station):

    """
    Use pycurl to query station details via FDSN.

    :param station: station site code
    :return: station latitude (dec. deg.), longitude (dec. deg.), and depth (+ve direction is down, in m)
    """

    # Query station information
    query = 'https://service.geonet.org.nz/fdsnws/station/1/query?station=' + station
    queryresult = curl(query)

    # Parse station information
    root = ET.fromstring(queryresult.decode('utf-8'))
    latitude = float(root[4][3][2].text)
    longitude = float(root[4][3][3].text)
    depth = -1 * float(root[4][3][4].text)

    return latitude, longitude, depth


def GeoNet_delta_station_query(station):

    """
    Use pycurl to query station details from GeoNet delta
    
    :param station: station site code
    :return: station latitude (dec. deg.), longitude (dec. deg.), and depth (+ve direction is down, in m)
    """

    query = 'https://raw.githubusercontent.com/GeoNet/delta/master/network/stations.csv'
    queryresult = curl(query)
    for row in queryresult.decode('utf-8').split('\n')[1:]:
        cols = row.split(',')
        if cols[0] == station:
            latitude = float(cols[3])
            longitude = float(cols[4])
            depth = -1 * float(cols[5])
            break

    return latitude, longitude, depth


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
        event = FDSN_event_query(eventID)[0]

        # Get origin and picks
        origin = event.origins[0]
        arrivals = origin.arrivals  # Picks used in origin
        picks = event.picks  # All picks for event

        # Extract all arrival time data from picks, but only for those used in the origin
        arrival_data = [[], [], []]
        for n in range(len(picks)):
            for m in range(len(arrivals)):
                if arrivals[m].time_weight and arrivals[m].time_weight > 0 and arrivals[m].phase in ['P', 'S']:
                    if picks[n].resource_id == arrivals[m].pick_id:
                        arrival_data[0].append(picks[n].waveform_id['station_code'])
                        arrival_data[1].append(arrivals[m].phase)
                        arrival_data[2].append(picks[n].time)
                        break

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
    arrival_time_data_header.sort()

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
    sites = []
    found = []
    for n in range(len(arrival_time_data_header)):

        # Get site name and details
        site = arrival_time_data_header[n].split('_')[0]

        # Do not duplicate site details
        if site in sites:
            continue

        network_data.append([])
        # Search the GeoNet FDSN database for the station location
        try:
            latitude, longitude, depth = FDSN_station_query(site)
        except:
            # If this fails, search the GeoNet delta database for the station location
            try:
                latitude, longitude, depth = GeoNet_delta_station_query(site)
            except:
                # If this fails, save the location details as None for later removal
                latitude, longitude, depth = None, None, None

        if latitude:
            # Convert WGS84 site geographic coordinates to WGS84 geodetic coordinates
            x, y, _ = convert_wgs84_geo_geod(latitude, longitude, -1 * depth)
            network_data[-1] = [site, x / 1000, y / 1000, depth / 1000,
                                latitude, longitude]
            found.append(True)
        else:
            # As before, save data so it can be removed
            network_data[-1] = [site, None, None, None, None, None]
            found.append(False)
        sites.append(site)

    num_found = found.count(True)
    if num_found != len(sites):
        # Remake the arrival time header and network data lists
        temp_ATDH = []
        temp_ND = []
        for n in range(len(found)):
            if found[n]:
                temp_ATDH.append(sites[n] + '_P')
                temp_ATDH.append(sites[n] + '_S')
                temp_ND.append(network_data[n])

        # Remake the arrival time data list
        temp_ATD = [[] for n in range(len(temp_ATDH))]
        for n in range(len(temp_ATDH)):
            for m in range(len(arrival_time_data_header)):
                if temp_ATDH[n] == arrival_time_data_header[m]:
                    for k in range(len(arrival_time_data)):
                        temp_ATD[n].append(arrival_time_data[k][m])

        arrival_time_data_header = temp_ATDH
        arrival_time_data = temp_ATD
        network_data = temp_ND

    return all_event_origins, arrival_time_data_header, arrival_time_data, network_data


def calculate_tt(grid_point, site_location, velocity_model, phase, angle_step=0.1):

    """
    Calculate the travel time of a seismic wave from a given grid point to a given site for a given 1D velocity model.
    Assumes ray travels in the plane containing the grid point, site, and the vertical axis.
    Note: this calculation is only valid when grid points and sites are spatially close, i.e.
    when the Earth-radial (vertical) direction is parallel for all grid points and sites. When this is not
    true, the spherical mode should be used with the script.
    :param grid_point: x,y,z position of a point in the grid
    :param site_location: x,y,z position of a site in the grid
    :param velocity_model: velocity model as parsed from velocity file
    :param phase: which phase to calculate travel time for: P or S
    :param angle_step: step size on incident ray angle for ray tracing, default is 0.01 degrees.
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

    # Find horizontal angle between grid point and station. Force this into the positive x and positive y right
    # triangle. Appropriate signs for values are set elsewhere.
    horangle = math.degrees(math.atan((grid_point[0] - site_location[1]) / (grid_point[1] - site_location[2])))

    # Set up angles to iterate
    start_angle = angle_step
    end_angle = 180 - angle_step
    num_steps = int(round((end_angle - start_angle) / angle_step)) + 1

    # Perform ray tracing to find the travel time between the grid point and site location for the first arriving ray
    horizontal_proximities = []
    travel_times = []
    for psi in np.linspace(start_angle, end_angle, num_steps):

        # Initialise ray tracing parameters
        n = nmax
        travel_time = 0
        current_angle = psi
        v_idx = ['P', 'S'].index(phase)

        # Collapse the 3D positions of the two points into a) the plane containing both points and the vertical axis
        # and b) the horizontal plane. In both cases the origin is the site location in the given plane.
        hpos = [grid_point[0] - site_location[1], grid_point[1] - site_location[2]]
        vpos = grid_point[2] - site_location[3]

        # Trace the ray up through velocity layers until it reaches the depth of the site
        while True:

            # Find if the site is in the velocity layer containing the grid point
            # If it is, set the vertical distance traveled in the layer as that between the grid point and the site,
            # otherwise, set the vertical distance traveled in the layer as that between the grid point and the
            # layer roof.
            if velocity_model[n][0] < site_location[3]:
                dz = vpos
            else:
                dz = vpos - (velocity_model[n][0] - site_location[3])

            # Increase the travel time by the time taken for the ray to traverse the hypotenuse of the
            # triangle with acute angle psi and opposite dz
            travel_time += (dz / math.sin(math.radians(current_angle))) / velocity_model[n][v_idx + 1]

            # Calculate horizontal distance traversed by ray over dz.
            # Positive dh means ray travelled in bearing direction,
            # negative dh means ray travelled opposite to bearing direction.
            dh = dz / math.tan(math.radians(current_angle))

            # Adjust ray head position due to distance traversed
            # dz, dy will be of the opposite sign to x and y as they are "approaching" the origin
            # as the ray travels. If dh is negative then dx and dy will be travelling away from the origin,
            # i.e. in the direction of x,y from the origin, rather than in the direction of the origin from xy.
            hpos[0] += math.sin(math.radians(horangle)) * dh
            hpos[1] += math.cos(math.radians(horangle)) * dh
            vpos -= dz

            # Once the ray reaches the same depth as the site, save the horizontal position and travel time of the ray
            if vpos <= 0.01:

                # Once the ray hits the depth of the site, the ray tracing ends. The horizontal distance between
                # the site location and the piercing point of the ray and the plane perpendicular to the vertical
                # axis at the depth of the site indicates how close the ray comes to the site. The ray with the
                # minimum horizontal distance will be the best solution to the ray tracing problem and its travel
                # time is taken as the travel time of a ray from the grid point to the site.
                horizontal_proximities.append(math.sqrt(hpos[0] ** 2 + hpos[1] ** 2))
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
                    exit()

    # The travel time between the grid point and the site is that of the ray which is closest to the site when it
    # reaches the depth of the site.
    travel_time = travel_times[horizontal_proximities.index(min(horizontal_proximities))]

    return travel_time


def generate_tt_grid(network_data, velocity_model, mode, test_origins = None,
                     xmin=-1, xmax=1, ymin=-1, ymax=1, zmin=-1, zmax=1, xstep=1, ystep=1, zstep=1):

    """
    Generate travel times from each grid cell to each site in the network for the given velocity model.
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
    :param mode: whether to calculate travel times in a "flat" mode or in a "spherical" mode.
                 The "flat" mode calculates travel times using ray tracing in a cartesian coordinate system
                 where the Earth-radial (vertical) direction is parallel for all grid points and sites.
                 The "spherical" mode references IASPEI91 travel times from the obspy.taup module.
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
            if mode == 'flat':
                gridx.append(float(test_origins[n][1]))
                gridy.append(float(test_origins[n][2]))
            elif mode == 'spherical':
                gridx.append(float(test_origins[n][-2]))
                gridy.append(float(test_origins[n][-1]))
            gridz.append(float(test_origins[n][3]))

    # Define every point in the grid and generate travel times for each grid point to each station
    # and append these in order behind the corresponding grid point in the nested lists
    grid_points = [[[[] for z in gridz] for y in gridy] for x in gridx]
    for i in range(len(gridx)):
        for j in range(len(gridy)):
            if test_origins and i != j:
                continue
            for k in range(len(gridz)):
                if test_origins and j != k:
                        continue
                grid_points[i][j][k] = [gridx[i], gridy[j], gridz[k]]
                for m in range(len(network_data)):
                    if mode == 'flat':
                        grid_points[i][j][k].append(calculate_tt(grid_points[i][j][k], network_data[m],
                                                                 velocity_model, phase='P'))
                        grid_points[i][j][k].append(calculate_tt(grid_points[i][j][k], network_data[m],
                                                                 velocity_model, phase='S'))
                    elif mode == 'spherical':
                        # Find angle between grid point and site
                        delta = math.degrees(2 *
                                             (math.asin(((math.sin(1/2 *
                                                                   math.radians(abs(grid_points[i][j][k][0] -
                                                                                    network_data[m][-2])))) ** 2 +
                                                         math.cos(math.radians(grid_points[i][j][k][0])) *
                                                         math.cos(math.radians(network_data[m][-2])) *
                                                         (math.sin(1/2 * math.radians(abs(grid_points[i][j][k][1] -
                                                                                          network_data[m][-1])))) **
                                                         2) ** (1/2))))

                        # Set receiver depth to 0 if it is negative; function does not allow depth below 0.
                        # However, increase source depth by receiver height to try reduce error due to this.
                        if network_data[m][3] < 0:
                            receiver_depth = 0
                            source_depth = grid_points[i][j][k][2] - network_data[m][3]
                        else:
                            receiver_depth = network_data[m][3]
                            source_depth = grid_points[i][j][k][2]

                        # Calculate travel times for common phases
                        arrivals = spherical_velocity_model.get_travel_times(source_depth_in_km=source_depth,
                                                                             receiver_depth_in_km=receiver_depth,
                                                                             distance_in_degree=delta,
                                                                             phase_list=['p', 'P', 's', 'S'])
                        p_tt, s_tt = None, None  # Prepare travel time variables
                        for arrival in arrivals:
                            if arrival.name == 'p' or arrival.name == 'P':
                                p_tt = arrival.time
                                break
                        for arrival in arrivals:
                            if arrival.name == 's' or arrival.name == 'S':
                                s_tt = arrival.time
                                break
                        grid_points[i][j][k].extend([p_tt, s_tt])

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
                break

    # Locate all events (only location method currently available is grid search)
    event_solutions = [[], [], [], [], [], [], [], [], []]
    for m in range(len(arrival_time_data)):
        origin_grid = [[[[] for z in range(len(grid_points[0][0]))]
                        for y in range(len(grid_points[0]))]
                       for x in range(len(grid_points))]

        # For each event, consider each grid cell as a possible location
        max_weight = 0
        weight_sum = 0
        for i in range(len(grid_points)):
            for j in range(len(grid_points[i])):
                for k in range(len(grid_points[i][j])):

                    # Calculate all origin times from the arrival time minus the travel time from
                    # the grid cell to each site
                    origin_times = []
                    for n in range(len(arrival_time_data[m])):
                        if (isinstance(arrival_time_data[m][n], datetime.datetime) and
                                grid_points[i][j][k][translation_chart[n]]):  # Ensure data is not nan or None
                            origin_times.append(arrival_time_data[m][n] -
                                                datetime.timedelta(seconds=grid_points[i][j][k][translation_chart[n]])
                                                )
                    if len(origin_times) == 0:  # If no data exists that is not nan is None
                        # Save weight, RMS error and mean origin time as None
                        origin_grid[i][j][k] = [None, None, None]
                    else:
                        origin_times.sort()

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
                        except:  # If it fails, then there is only one data point, so it's weight should be very high
                            weight = 9999

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
                            if origin_grid[i][j][k][0]:
                                if origin_grid[i][j][k][0] > c * max_weight:
                                    ijk.append([i, j, k])
                                    c_weight_sum += origin_grid[i][j][k][0]
                try:
                    c_weight_sum /= weight_sum
                except:  # Catch when the weight sum is 0, i.e. there is no OT data at all for the event
                    break
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
                x_err = float('nan')
                y_err = float('nan')
                z_err = float('nan')
                ots_err = float('nan')
                mid_x = float('nan')
                mid_y = float('nan')
                mid_z = float('nan')
                mid_ot = float('nan')
                rms = float('nan')
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


def test_test_origins(method, arrival_time_data, arrival_time_data_header, grid_points, grid_header, test_origins):

    """
    Perform origin testing using the arrival time data, search grid, location methhod and test origins given.
    :return: earthquake origins as lists of latitude, longitude, elevation (i.e. the test origin) and origin
             times and the RMS error associated with each origin calculated as the time difference between its
             origin time and that of the test origin.
    """

    earthquake_origins, rms_errors = [], []
    for n in range(len(test_origins)):
        earthquake_origins.append(globals()[method]([arrival_time_data[n]],
                                                    arrival_time_data_header,
                                                    [[[grid_points[n][n][n]]]],
                                                    grid_header))

        # Calculate RMS error as the time difference between the test origin time and that from the grid search
        try:
            rms_errors.append(abs((datetime.datetime.strptime(test_origins[n][0],
                                                              '%Y-%m-%dT%H:%M:%S.%fZ') -
                                   earthquake_origins[-1][3][0]).total_seconds()))
        except:  # If this fails, assume it's due to a different format in the test origin
            try:
                rms_errors.append(abs((datetime.datetime.strptime(test_origins[n][0],
                                                                  '%Y-%m-%dT%H:%M:%SZ') -
                                       earthquake_origins[-1][3][0]).total_seconds()))
            except:  # If this fails, assume it is because there is no valid data for the event
                rms_errors.append(float('nan'))
    return earthquake_origins, rms_errors


if __name__ == "__main__":

    # If the code is executed directly, parse arguments from command line using argparse
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
                                                         'depth (m, +ve direction is down),'
                                                         'origin time (IS08601 format, including milliseconds).')
    parser.add_argument('--network-file', type=str, help='Comma-separated file containing network details. '
                                                         'First row is header with columns: '
                                                         'site,latitude,longitude,elevation. '
                                                         'All other rows contain site details as: '
                                                         'site code,'
                                                         'latitude (decimal degrees),'
                                                         'longitude (decimal degrees),'
                                                         'elevation (m, +ve direction is up). ')
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
    parser.add_argument('--mode', type=str, help='Which mode to use in travel time calculation: flat or spherical.')
    parser.add_argument('--method', type=str, help='Earthquake location method to use, options are: grid_search.')
    args = parser.parse_args()

    # Parse files and parameters
    arrival_time_data, arrival_time_data_header, grid_points, grid_header, test_origins = parse_files(
        eventid_file=args.eventid_file,
        test_origins=args.test_origins,
        network_file=args.network_file,
        velocity_model=args.velocity_model,
        grid_parameters=args.grid_parameters,
        grid_file=args.grid_file,
        search_mode=args.mode)
    method = args.method

    # Perform earthquake location
    if method == 'grid_search':  # Code only supports one method currently
        if args.test_origins:  # ...and has only been tested for one use case
            earthquake_origins, rms_errors = test_test_origins(method,
                                                               arrival_time_data,
                                                               arrival_time_data_header,
                                                               grid_points,
                                                               grid_header,
                                                               test_origins)
