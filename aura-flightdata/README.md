# AuraUAS aura-flightdata

Python libraries to load and interpolate a variety of flight data
formats.  These libs are primarily intended to support higher level
flight post processing scripts.

Included in this package is a weather/forecast front end to
forecast.io.  You will need to register for your own free apikey to
use this module.  This allows scripts to dump out a weather summary at
the time and location of the flight, even if the flight was weeks or
months or even years ago.

Data file formats supported include:
* AuraUAS (hdf5 and csv variants)
* PX4 sdlog2, ulog
* Sentera, Sentera2
* UMN Goldy 1 (matlab)
* UMN Goldy 3 (hdf5)

(With some effort) it is possible to extend this library to support
flight data from other autopilot systems.  Some of the tools built on
top of the flight data loader include:

* Running the UMN/AEM/UAV lab EKF and comparing it's output to the
  onboard native EKF.
* Generating HUD overlays on action cam flight footage.
* Various automatic magnetometer and IMU temp calibration tools.
* Mapping, stitching, and geotagging tools.
* Synthetic airspeed estimator.
* Wind estimator.

# Quick Install

    $ sudo python3 ./setup.py install

# Example Usage:

```python
    #!/usr/bin/python3
    from aurauas.flightdata import flight_loader, flight_interp
    data, flight_format = flight_loader.load("/flight/data/log/path")
    iter = flight_interp.IterateGroup(data)
    for i in range(iter.size()):
        record = iter.next()
        imu = record['imu']
        gps = record['gps']
        print(imu['time'], imu['p'], imu['q'], imu['r])
```