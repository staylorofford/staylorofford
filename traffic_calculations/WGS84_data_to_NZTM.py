"""
Convert csv file(s) containing WGS84 data coordinates into NZTM csv file(s) using ogr2ogr
"""

import argparse
import glob
import os
import shutil
import subprocess

if __name__ == "__main__":

    # Parse arguments

    parser = argparse.ArgumentParser(description='Convert csv file(s) in a folder from WGS84 to NZTM')
    parser.add_argument('file_directory', type=str, help="Absolute path to directory containing time series csv files "
                                                         "in WGS84")
    args = parser.parse_args()

    file_dir = args.file_directory
    data_files = glob.glob(file_dir + '*.csv')

    # Convert each data file to a NZTM csv file

    for data_file in data_files:
        os.chdir(file_dir)
        filename = data_file.split('/')[-1]
        p = subprocess.Popen('ogr2ogr -f \"ESRI Shapefile\" ' + filename[:-4] + '.dbf ' + filename, shell=True)
        p.wait()
        with open(filename[:-4] + '.vrt', 'w') as openfile:
            pass
        with open(filename[:-4] + '.vrt', 'a') as openfile:
            openfile.write('<OGRVRTDataSource>\n' +
                           '  <OGRVRTLayer name=\"' + filename[:-4] + '\">\n' +
                           '    <SrcDataSource>' + filename + '</SrcDataSource>\n' +
                           '    <SrcLayer>' + filename[:-4] + '</SrcLayer>\n' +
                           '    <LayerSRS>EPSG:4326</LayerSRS>\n' +
                           '    <GeometryField encoding=\"PointFromColumns\" x=\"lon\" y=\"lat\"/>' +
                           '  </OGRVRTLayer>\n' +
                           '</OGRVRTDataSource>\n')
        p = subprocess.Popen('ogr2ogr -f \"ESRI Shapefile\" ' + filename[:-4] + ' ' + filename[:-4] + '.vrt',
                             shell=True)
        p.wait()
        os.chdir(file_dir + filename[:-4])
        p = subprocess.Popen('ogr2ogr -t_srs epsg:2193 -s_srs epsg:4326 ./' + filename[:-4] + ' ./' + filename[:-4] +
                             '.shp', shell=True)
        p.wait()
        os.chdir('./' + filename[:-4])
        p = subprocess.Popen('ogr2ogr -f CSV ' + filename + ' ' + filename[:-4] + '.shp -lco GEOMETRY=AS_XYZ',
                             shell=True)
        p.wait()

        os.chdir(file_dir)
        os.remove(filename)
        os.remove(filename[:-4] + '.dbf')
        os.remove(filename[:-4] + '.vrt')
        os.rename('./' + filename[:-4] + '/' + filename[:-4] + '/' + filename, './' + filename)
        shutil.rmtree(filename[:-4])
