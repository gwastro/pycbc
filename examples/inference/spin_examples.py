import matplotlib as mpl
import matplotlib.pyplot as plt
from pycbc.inference import distributions
import pycbc.coordinates as co
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

# Let's plot what we've made so far with histograms.
# Choose 20 mass bins for the histograms.
n_bins = 20

plt.figure(figsize=(10,10))
plt.subplot(2, 2, 1)
plt.hist(spinx, bins = n_bins)
plt.title('Spin x samples')

plt.subplot(2, 2, 2)
plt.hist(spiny, bins = n_bins)
plt.title('Spin y samples')

plt.subplot(2, 2, 3)
plt.hist(spinz, bins = n_bins)
plt.title('Spin z samples')

plt.tight_layout()
plt.show()
