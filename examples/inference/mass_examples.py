import matplotlib.pyplot as plt
from pycbc.inference import distributions
import numpy as np

# Create a mass distribution object that is uniform between 0.5
# and 1.5 solar masses.
mass1_distribution = distributions.Uniform(mass1=(0.5, 1.5))
# Take 100000 random variable samples from this uniform mass distribution.
mass1_samples = mass1_distribution.rvs(size=1000000)

# Create another mass distribution object that is Gaussian between
# 0.5 and 1.5 solar masses with a mean of 1.2 solar masses and a
# standard deviation of 0.15 solar masses. Gaussian takes the variance
# as an output so square the standard deviation.
mass2_distribution = distributions.Gaussian(mass2=(0.5, 1.5, 1.2, 0.15*0.15))
# Take 100000 random variable samples from this gaussian mass distribution.
mass2_samples = mass2_distribution.rvs(size=1000000)

samples = np.ndarray(shape=(1000000), dtype=float)
for i in range(0,1000000):
    samples[i] = mass2_samples[i][0]
# We can make pairs of distributions together, instead of apart.
two_mass_distributions = distributions.Uniform(mass1=(1.6, 3.0),
                                               mass2=(1.6, 3.0))
two_mass_samples = two_mass_distributions.rvs(size=1000000)

# Let's plot what we've made so far with histograms.
# Choose 20 mass bins for the histograms.
n_bins = 20

# Make some subplot objects
fig, axes = plt.subplots(nrows=2, ncols=2)
ax0, ax1, ax2, ax3, = axes.flat

ax0.hist(mass1_samples['mass1'], bins = n_bins)
ax1.hist(samples, bins = n_bins)
ax2.hist(two_mass_samples['mass1'], bins = n_bins)
ax3.hist(two_mass_samples['mass2'], bins = n_bins)

ax0.set_title('Mass 1 samples')
ax1.set_title('Mass 2 samples')
ax2.set_title('Mass 3 samples')
ax3.set_title('Mass 4 samples')

plt.tight_layout()
plt.show()
