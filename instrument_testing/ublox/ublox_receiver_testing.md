# U-Blox Receiver Testing

## Project Overview

The ublox EVK-M8U GNSS receiver was identified as a potential low-cost GNSS receiver that could be used in GeoNet temporary equipment deployments. A project was established to test the receiver against the current CORS GNSS receiver following the testing process used in the 2018/19 GNSS receiver RFP to evaluate its suitability as an instrument in the GeoNet network.

### Project Purpose

During a geohazard response, GeoNet is occasionally involved in the temporary deployment of geodetic data collection equipment. Recent temporary deployments include those on the Fox Glacier landslide and those in the area surrounding the Kaikoura Earthquake fault ruptures. Currently such temporary deployments use the standard CORS GNSS equipment, this being suitable for the long-term, high-accuracy observations that GeoNet requires for its ongoing deformation monitoring. While capable of providing the highest quality data, this equipment often exceeds the requirements of such deployments, not to mention the high cost of such instruments and the risk associated with their deployment in active geohazard areas. This project aims to test a low-costⁱ GNSS receiver (the u-blox EVK-M8U) and Raspberry Pi as an alternative to the standard CORS GNSS receiver and datalogger (the datalogger in this case is integrated with the receiver) for use in temporary or high-risk deployments, as well as other applications that are not yet realised.

ⁱ Low cost is a field ready unit with a complete cost of under $1000.

### Project Scope

To test the low-cost GNSS receiver we propose the following tests:
1.	Configure receiver to closest equivalent configuration to standard GeoNet CORS GNSS receiver
2.	Log raw GNSS observables* on receiver when connected to a standard GeoNet CORS GNSS antenna.
3.	Log data from receiver on storage in the form of hourly files that are either in a RINEX file format or one that can be converted to it using existing software.
4.	Retrieve data from the receiver in a way that is compatible with current data reception methods.
5.	Derive a position using the logged data and differential-GNSS location techniques and check if the position is of the required accuracy^.
6.	Run the receiver in a field-style deployment for a week or more and check the logged data for quality and completeness.
7.	Record the power consumption of a fully operational receiver over a 24-hour period.
8.	Perform a security assessment of all devices used in the project to ensure compliance with GeoNet's network security standards.
9.	Assess and test only the Raspberry Pi and u-blox device only.

\* Raw GNSS observables are here defined as the L1, P1, S1 observables as a minimum for each tracked constellations.

^ The “required accuracy” is deployment-dependent but is here taken to be of 1 cm or greater as this would be appropriate for the “fast-moving” (>1 cm / day) environments the low-cost receiver is envisioned to be deployed in.

## Testing Documentation

### Equipment Setup

The equipment was setup as follows:
- The u-blox receiver and Raspberry Pi were connected to a mains-5V adapter via USB cable and powered over USB.
- The u-blox receiver was connected to the Raspberry Pi via the serial port (on the receiver) and a USB port (on the Raspberry Pi).
- The u-blox receiver was connected to a Zephyr 2 antenna via the RF IN connector.
- The Raspberry Pi was connected to the GeoNet laboratory LAN via an ethernet cable.
- The Zephyr 2 antenna was mounted on the roof of the Avalon campus of GNS Science and connected to the receiver via an established cable from the roof to the testing location.

This equipment setup was sufficient to perform most of the instrument testing, but for the accuracy testing and field testing alternative setups would be required - the testing location used in this setup is not representative of a field deployment and the antenna mount used is not stable enough to allow reliable estimates of receiver accuracy.

### Raspberry Pi Configuration

This section describes how the The Raspberry Pi was configured from the default state. 

#### Network Configuration

The first three points relate to the network settings of the Raspberry Pi and were completed using the GUI avaialble over HDMI from the Raspberry Pi due to the relative ease of changing such settings here and the inability to access the receiver by other methods until these changes were made. 
- To allow connection to the GeoNet laboratory LAN without disruption, the wireless capability of the Raspberry Pi was disabled.
- In the IPv4 networking, DHCP was disabled and a static IP was chosen. This IP address, subnet mask, and default gateway were chosen so that the Raspberry Pi would be accessible from outside the GeoNet laboratory LAN via jumphosts. This networking setup emulates that of a remote field site and allows testing team members to access the Raspberry Pi and receiver from anywhere with access to the jumphosts. 
- To allow connections to the Raspberry Pi over SSH, SSH was enabled in the Raspberry Pi Preferences.

All subsequent configuration changes were made over the command line by using SSH to connect to the Raspberry Pi from the jumphosts.

#### Python Setup and Data Logging

As the receiver has no logging capabilities, all logging must be performed externally. To do this Python was used.

The Raspberry Pi comes with Python3 installed, but to make sense of the data that the receiver outputs a non-standard Python module was used - the [pyUBX](https://github.com/mayeranalytics/pyUBX) module. This module was copied onto the Raspberry Pi from the jumphost used and all code that required functionality from the module was run in the module's root directory. The module itself was setup as a folder in the home directory of the Raspberry Pi with no alterations to its files or structure.  

Once `pyUBX` was setup on the Raspberry Pi, a Python script `read_serial.py` was executed in the `pyUBX` directory. This script continually read data from the Raspberry Pi's USB0 port which, in this case, corresponded to the port receiving data from the serial port on the receiver. To find which port on the Raspberry Pi was connected to the receiver, the command `dmesg | grep tty` was used. The interface at the end of the `usb 1-1.2: pl2303 converter now attached to ttyUSB0` row in the command output is that receiving the serial data from the receiver, i.e. `USB0`.

The `read_serial.py` script is:

```python
import serial
import time
from UBXManager import UBXManager
      
ser = serial.Serial(port='/dev/ttyUSB0',
                    baudrate=9600,
                    timeout=0
                   )

manager = UBXManager(ser, 
                     debug=True)
manager.start()
```
And is run with `python3 read_serial.py` in the `pyUBX` root directory. The script is ended by using `CTRL+C`. One can save the script output by executing it as `python 3 read_serial.py >> serial.log`. Note that the baud rate here is the receiver default. If the receiver is configured such that it sends many messages, the default baud rate may be insufficient. In this case the maximum baud rate (921600) can be used, but note that it must be configured on the receiver and in the `read_serial.py` script. If there is a mismatch in baud rate between the two ends of the serial communications then no data will be received. If not all messages are being received then it is possible that the baud rate is too low.

The receiver is capable of sending two data formats. The data format the `read_serial.py` reads is UBX messages, this being the proprietary format designed by the receiver manufacturer. UBX messages  can contain all data produced by the receiver, notably the raw GNSS observables. The other data format the receiver can send is NMEA messages. These contain information such as receiver position, tracked signals, signal quality, and other high-level GNSS data. NMEA messages cannot contain raw GNSS observables. By default, the receiver will send only NMEA messages. UBX Messages must be enabled in the receiver configuration before they will be sent and they are hexadecimally encoded within the UBX message structure. Parsing UBX messages from the receiver data stream is difficult, and rather than performing this directly in Python it proved easier to use the third-party module `pyUBX`. NMEA messages are comparatively easy to parse, but do not contain the data required for testing and high-accuracy applications.

A version of the `read_serial.py` script that can read NMEA messages is:

```python
import pynmea2
import serial
import time
      
ser = serial.Serial(port='/dev/ttyUSB0',
                    baudrate = 9600,
                    parity=serial.PARITY_NONE,
                    stopbits=serial.STOPBITS_ONE,
                    bytesize=serial.EIGHTBITS,
                    timeout=1
                   )

counter=0
while 1:
        line = ser.readline().decode('ascii', errors='replace')
        print('line is: ' + line)        
       
        nmeaobj = pynmea2.parse(line)
        try:
            for i in range(len(nmeaobj.fields)):
                print(str(nmeaobj.fields[i][0]) + ': ' + str(nmeaobj.data[i]))
        except:
            pass
```
which is presented without any further documentation. It is not intended that NMEA messages from the receiver will ever be used granted the intended use cases of the receiver.
 
Using `pyUBX` it is also possible to send UBX messages to the receiver over the serial connection, however this functionality is yet to be tested. It is imagined that `pyUBX` or equivalent Python functions will be used to configure the receiver in future.

### Receiver Configuration

The receiver was configured using [u-centre](https://www.u-blox.com/en/product/u-center), a Windows-only software that u-blox offers for configuring and logging data from the receiver. The receiver was configured from its default settings using a .ubx file that Robert Odolinski (Otago University) provided us by using u-centre on a Windows laptop connected to the receiver over the serial cable. Robert's .ubx file contained the UBX messages required to configure the receiver such that it would track L1 observables for GPS, GLONASS, and Galileo constellations at a 1 second sampling interval.

It is possible, and remarkably easy, to log data from UBX messages using u-centre. This highly desirable functionality was sadly not a realistic option for data logging given the Linux-based nature of GeoNet and our systems.     

While appropriate for configuring the receiver to a standard suitable for initial testing, any required use of u-centre to configure or interact with the receiver was considered unacceptable. As such, attempts were made to achieve the same functionality we required of u-centre - to configure the receiver and activate data logging on the receiver - in a Linux environment.

Data logging on the receiver did not function as expected from using geodetic-grade GNSS receivers. On such receivers data logging begins as soon as the receiver begins tracking satellites. On the u-blox receiver data logging (that is, the logging of raw GNSS observables) only began after the operator began a data logging session in u-centre. There is little transparency on what occurs after one begins such a session. It appeared that the receiver was not configured to produce UBX messages until a data logging session was started, and that starting such a session actually involved reconfiguring the receiver to log UBX messages containing the raw GNSS observables. As part of attempting to configure the receiver in a Linux environment we also tried to enable data logging from a Linux environment by emulating the configuration changes u-centre applied when a data logging session was enabled in that software. 

#### Linux-Based Receiver Configuration

By far the most challenging aspect of the testing performed was interacting with the receiver in a Linux-environment. The receiver only accepts UBX messages coming over the serial port as input. The principal issue with this is that UBX messages are not human-readable, or even human-intelligable: they are hexadecimally encoded and the packet structure varies with every message type. When one considers that 235 of the 400 pages of the u-blox receiver manual describe UBX messages the magnitude of the task becomes crushingly apparent. 

Attempts were made to translate the desired settings for each configuration element into UBX messages and send these to the receiver over the serial port but no success was had here. The principle issue was that the checksum included in each UBX message, this being the suffix to every message and being derived from the complete preceding message, was always at odds with that produced by the receiver for every message sent. Ultimately this problem was not able to be resolved in a timely manner and the `pyUBX` Python module was identified as an alternative to developing our own bespoke solution. Sadly the receiver testing did not progress to a point where `pyUBX` could be used in this capacity.

Granted that `pyUBX` could successfully send UBX messages to the receiver the issue of translating desired settings for each receiver configuration element remained. To do this one has to study the UBX message that can change the desired configuration element(s), figure out which numbers would produce the desired change, and then translate the numeric string corresponding to the desired configuration into the hexadecimal encoding the message requires. Rather than working through all such messages in this way we attempted to define the desired messages from those sent to the receiver from u-centre when we made changes to the receiver configuration using that software. This approach showed great promise. If we would be able to define all UBX messages required to configure the receiver to our desired standard then it would be possible to configure a receiver in seconds using a script that sends UBX messages over a serial connection.

In the process of trying to configure the receiver in this way it appears that some corruption of the receiver configuration occurred. The configuration approach was highly experimental, and many UBX messages were sent to the receiver in an attempt to reproduce the "enable data logging session" functionality that existed in u-centre. Prior to the conclusion of this testing, the `serial_read.py` script produces an output of the form:

```bash
pi@raspberrypi:~ $ python3 serial_read.py 
Writing log to UBX.log
UBX ERR 10:02 No parse, "'UBX'", payload=bd ba 22 01 18 18 00 00 e0 ff ff 10 b7 03 00 11 cf 27 00 12 bd ba 22 01
UBX ERR 10:02 No parse, "'UBX'", payload=b8 ba 22 01 18 20 00 00 e4 fd ff 0e cb fc ff 0d 80 ff ff 05 30 0a 00 0c b8 ba 22 01
UBX ERR 10:02 No parse, "'UBX'", payload=20 bb 22 01 18 18 00 00 de ff ff 10 b1 03 00 11 c2 27 00 12 20 bb 22 01
UBX ERR 10:02 No parse, "'UBX'", payload=1b bb 22 01 18 20 00 00 9c fe ff 0e 5b fd ff 0d a2 ff ff 05 30 0a 00 0c 1b bb 22 01
UBX ERR 10:02 No parse, "'UBX'", payload=82 bb 22 01 18 18 00 00 e1 ff ff 10 ba 03 00 11 d0 27 00 12 82 bb 22 01
UBX ERR 10:02 No parse, "'UBX'", payload=7d bb 22 01 18 20 00 00 48 fe ff 0e ea fc ff 0d 9f ff ff 05 31 0a 00 0c 7d bb 22 01
UBX ERR 10:02 No parse, "'UBX'", payload=e5 bb 22 01 18 18 00 00 dc ff ff 10 bb 03 00 11 c6 27 00 12 e5 bb 22 01
UBX ERR 10:02 No parse, "'UBX'", payload=e0 bb 22 01 18 20 00 00 06 fe ff 0e 0f fd ff 0d 7d ff ff 05 30 0a 00 0c e0 bb 22 01
UBX ERR 10:02 No parse, "'UBX'", payload=47 bc 22 01 18 18 00 00 e1 ff ff 10 b1 03 00 11 c7 27 00 12 47 bc 22 01
UBX ERR 10:02 No parse, "'UBX'", payload=42 bc 22 01 18 20 00 00 0c fe ff 0e e4 fc ff 0d a2 ff ff 05 31 0a 00 0c 42 bc 22 01
UBX ERR 10:02 No parse, "'UBX'", payload=a9 bc 22 01 18 18 00 00 df ff ff 10 b7 03 00 11 c6 27 00 12 a9 bc 22 01
UBX ERR 10:02 No parse, "'UBX'", payload=a4 bc 22 01 18 20 00 00 8d fd ff 0e 2f fd ff 0d d8 ff ff 05 30 0a 00 0c a4 bc 22 01
UBX ERR 10:02 No parse, "'UBX'", payload=0c bd 22 01 18 18 00 00 e4 ff ff 10 ae 03 00 11 c1 27 00 12 0c bd 22 01
UBX ERR 10:02 No parse, "'UBX'", payload=07 bd 22 01 18 20 00 00 54 fe ff 0e c2 fc ff 0d 80 ff ff 05 30 0a 00 0c 07 bd 22 01
``` 
for UBX messages of many types. It appeared that either `pyUBX` could not parse the UBX messages the receiver was configured to send, or that the receiver, or cable, or Raspberry Pi was somehow corrupting the messages sent. The receiver testing was terminated here.
  