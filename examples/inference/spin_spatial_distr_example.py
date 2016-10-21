import matplotlib as mpl
import matplotlib.pyplot as plt
from pycbc.inference import distributions
import pycbc.coordinates as co
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# We can choose any bounds between 0 and pi
# for this distribution but in units of pi
# so we use between 0 and 1
theta_low = 0.
theta_high = 1.

# Units of pi for the bounds of the azimuthal angle
# which goes from 0 to 2 pi
phi_low = 0.
phi_high = 2.

# Create a distribution object from distributions.py
# Here we are using the Uniform Solid Angle function
# which takes theta = polar_bounds(theta_lower_bound
# to a theta_upper_bound), and then phi = azimuthal_
# bound(phi_lower_bound to a phi_upper_bound).
uniform_solid_angle_distribution = distributions.UniformSolidAngle(
                                          polar_bounds=(theta_low,theta_high),
                                          azimuthal_bounds=(phi_low,phi_high))

# Now we can take a random variable sample from that
# distribution. In this case we want 50000 samples.
solid_angle_samples = uniform_solid_angle_distribution.rvs(size=50000)

# Make a spin 1 magnitude since solid angle is only
# 2 dimensions and we need a 3rd dimension for a 3D
# plot that we make later on
spin_mag = np.ndarray(shape=(50000), dtype=float)

for i in range(0,50000):
    spin_mag[i] = 1.

# Use the pycbc.coordinates as co
# spherical_to_cartesian function to convert from 
# spherical polar coordinates to cartesian coordinates
spinx, spiny, spinz = co.spherical_to_cartesian(spin_mag,
                                                solid_angle_samples['phi'],
                                                solid_angle_samples['theta'])

# Let's plot the spherical distribution of spins
# to make sure that we distributed properly across
# the surface of a sphere.

fig = plt.figure(figsize=(20,20))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(spinx, spiny, spinz, s=1)
plt.show()
