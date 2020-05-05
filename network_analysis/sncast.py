#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ------------------------------------------------------------------
# Filename: sncast.py
#  Purpose: Seismic Network Capability Assessment Software Tool (SNCAST)
#   Author: Martin Möllhoff, DIAS
# Citation: Möllhoff, M., C.J. Bean and B. Baptie (2019), SN-CAST: Seismic Network Capability Assessment Software Tool
#           for Regional Networks – Examples from Ireland, Journal of Seismology, doi: 10.1007/s10950-019-09819-0.
#
#    Copyright (C) 2019 Martin Möllhoff
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#    email:         martin@dias.ie
#    snail mail:    Martin Möllhoff, DIAS, 5 Merrion Square, Dublin 2, Ireland
#    web:	          www.dias.ie/martin  www.insn.ie
# --------------------------------------------------------------------

import numpy as np
from obspy.signal.util import util_geo_km
from math import pow, log10, sqrt
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import geopandas as gpd
import pandas as pd


def minML(filename, dir_in='./', lon0=165, lon1=185, lat0=-50, lat1=-30, dlon=0.1,
          dlat=0.1, stat_num=10, snr=10, foc_depth=5, region='CAL', mag_min=2.0, mag_delta=0.1):

    """
    This routine calculates the geographic distribution of the minimum 
    detectable local magnitude ML for a given seismic network. Required 
    input is a file with .dat extension containg four comma separated
    columns containing for each seimic station:

         longitude, latitude, noise [nm], station name
    e.g.: -7.5100, 55.0700, 0.53, IDGL

    The output file *.grd lists in ASCII xyz format: longitud, latitude, ML
  
    Optional parameters are:

    :param  dir_in:	full path to input and output file
    :param  lon0:	minimum longitude of search grid
    :param  lon1:	maximum longitude of search grid
    :param  lat0:	minimum latitude of search grid
    :param  lat1:	maximum latitude of search grid
    :param  dlon:	longitude increment of search grid
    :param  dlat:	latitude increment of search grid
    :param  stat_num:	required number of station detections
    :param  snr:	required signal-to-noise ratio for detection
    :param  foc_depth:  assumed focal event depth
    :param  region:	locality for assumed ML scale parameters ('UK' or 'CAL')
    :param  mag_min:	minimum ML value for grid search
    :param  mag_delta:  ML increment used in grid search
    """

    # region specific ML = log(ampl) + a*log(hypo-dist) + b*hypo_dist + c
    if region == 'NZ':  # Ristau magnitude
        a = -1.49
        b = -1.27E-3
        c = 0.29
    elif region == 'UK':  # UK scale, Ottemöller and Sargeant (2013), BSSA, doi:10.1785/0120130085
        a = 0.95
        b = 0.00183
        c = -1.76
    elif region == 'CAL':  # South. California scale, IASPEI (2005),
                           # www.iaspei.org/commissions/CSOI/summary_of_WG_recommendations_2005.pdf
        a = 1.11
        b = 0.00189
        c = -2.09

    # read in data, file format: "LON, LAT, NOISE [nm], STATION"
    array_in = np.genfromtxt('%s/%s.dat' % (dir_in, filename), dtype=None, delimiter=",", filling_values=0)
    lon = ([t[0] for t in array_in])
    lat = [t[1] for t in array_in]
    noise = [t[2] for t in array_in]
    stat = [t[3] for t in array_in]
    corr = [t[4] for t in array_in]
    # grid size
    nx = int((lon1 - lon0) / dlon) + 1
    ny = int((lat1 - lat0) / dlat) + 1
    # open output file:
    f = open('%s/%s-stat%s-foc%s-snr%s-%s.grd' %(dir_in,
                                                 filename,
                                                 stat_num,
                                                 foc_depth,
                                                 snr,
                                                 region), 'wb')
    mag = []
    for ix in range(nx):  # Loop through longitude increments
        ilon = lon0 + ix * dlon
        for iy in range(ny):  # Loop through latitude increments
            ilat = lat0 + iy * dlat
            j = 0
            for jstat in stat:  # Loop through stations

                # Calculate hypocentral distance in km
                dx, dy = util_geo_km(ilon, ilat, lon[j], lat[j])
                hypo_dist = sqrt(dx ** 2 + dy ** 2 + foc_depth ** 2)

                # Find smallest detectable magnitude
                ampl = 0.0
                m = mag_min - mag_delta
                while ampl < snr * noise[j]:
                    m = m + mag_delta
                    #ampl = pow(10, (m - a * log10(hypo_dist) - b * hypo_dist - c ))
                    ampl = pow(10, (m - a * log10(hypo_dist) - b * hypo_dist - c - corr[j]))
                mag.append(m)
                j = j + 1   
            # Sort magnitudes in ascending order
            mag = sorted(mag)
            # Write out longitude, latitude and smallest detectable magnitude
            f.write("".join(str(ilon) + " " + str(ilat) + " " + str(mag[stat_num-1]) + "\n").encode())
            del mag[:]
    f.close()


def PlotminML(filename):

    """
    Plot the output on a  grid 
    """

    DAT = pd.read_csv(filename, sep=" ", header=None)
    x = DAT.iloc[:, 0].values
    y = DAT.iloc[:, 1].values
    z = DAT.iloc[:, 2].values

    def plot_contour(x, y, z, resolution=50, contour_method='linear'):

       resolution = str(resolution) + 'j'
       X, Y = np.mgrid[min(x):max(x):complex(resolution), min(y):max(y):complex(resolution)]
       points = [[a, b] for a, b in zip(x, y)]
       Z = griddata(points, z, (X, Y), method=contour_method)

       return X, Y, Z

    X, Y, Z = plot_contour(x, y, z, resolution=50, contour_method='linear')

    fig, ax = plt.subplots(figsize=(13, 8))
    cm = plt.cm.get_cmap('inferno_r')
    mesh = ax.pcolormesh(X, Y, Z, cmap=cm, zorder=1)
    CS = ax.contour(X, Y, Z, levels=[2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5], linewidths=1)
    ax.clabel(CS, inline=1, fontsize=10)
    ax.grid(color='k', alpha=0.2, zorder=2)
    cb = fig.colorbar(mesh)
    cb.set_label('magnitude', labelpad=15, fontsize=12, rotation=270)

    return fig, ax


# Calculate SN-CAST data
# minML('NZ_stations', foc_depth=5)
minML('NZ_CorrStation', foc_depth=5)

# Plot SN-CAST data
fig, ax = PlotminML('NZ_CorrStation-stat10-foc5-snr10-CAL.grd')
# Plot reference data
#station_metadata = pd.read_csv('NZ_stations.dat', header=None)
station_metadata = pd.read_csv('NZ_CorrStation.dat', header=None)
#station_metadata.columns = ['Longitude', 'Latitude', 'Noise', 'Station']
station_metadata.columns = ['Longitude', 'Latitude', 'Noise', 'Station', 'MLr_Corr']
outlines = gpd.read_file('nz-coastlines-and-islands-polygons-topo-150k.shp')
outlines.boundary.plot(color=None, edgecolor='k', ax=ax, linewidth=1, zorder=3)  # Plot NZ outline
ax.scatter(station_metadata['Longitude'],
           station_metadata['Latitude'],
           color='white',
           edgecolors='black',
           alpha=0.8,
           zorder=4)  # Plot station positions
plt.xlim([165, 185])
plt.ylim([-50, -30])
plt.xlabel('longitude', labelpad=10, fontsize=12)
plt.ylabel('latitude', rotation=90, labelpad=10, fontsize=12)
plt.title('SN-CAST: theoretical lowest magnitude of detection in New Zealand,\n' +
          'focal depth: 5 km, number of detections: 10, SNR for detection: 3',
          y=1.03)
plt.show()
