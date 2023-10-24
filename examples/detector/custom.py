import matplotlib.pyplot as plt
from pycbc.detector import add_detector_on_earth, Detector
import pycbc.psd
import numpy as np

# Set up potential Cosmic Explorer detector locations

# 40 km detector
lon = -125 / 180.0 * np.pi
lat = 46 / 180.0 * np.pi
yangle = 100.0 / 180.0 * np.pi 
# yangle is the rotation clockwise from pointing north at 0
# xangle can also be specified and allows for detectors that don't have
# 90 degree opening between arms. By default we assume xangle is yangle + pi/2
add_detector_on_earth("C4", lon, lat, yangle=yangle,
                      xlength=40000, ylength=40000)

# 20 km detector
# Arm length is optional, but if provided, you can accurately calcuale
# high-frequency corrects if you provide a frequency argument to the
# antenna pattern method
lon = -94 / 180.0 * np.pi
lat = 29 / 180.0 * np.pi
yangle = 160.0 / 180.0 * np.pi
add_detector_on_earth("C2", lon, lat, yangle=yangle,
                      xlength=20000, ylength=20000)
 
ra, dec = np.meshgrid(np.arange(0, np.pi*2.0, .1), 
                      np.arange(-np.pi / 2.0, np.pi / 2.0, .1))
ra = ra.flatten()
dec = dec.flatten()

pol = 0
time = 1e10 + 8000 # A time when ra ~ lines up with lat/lon coordinates

for d in [Detector("C4"), Detector("C2")]:
    fp, fc = d.antenna_pattern(ra, dec, pol, time)

    plt.figure()
    plt.subplot(111, projection="mollweide")
    ra[ra>np.pi] -= np.pi * 2.0
    plt.scatter(ra, dec, c=fp**2.0 + fc**2.0)
    plt.title("Mollweide")
    plt.grid(True)
    plt.show()
