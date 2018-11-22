import matplotlib.pyplot as plt
from pycbc import distributions

# Create a mass distribution object that is uniform between 0.5 and 1.5
# solar masses.
mass1_distribution = distributions.Uniform(mass1=(0.5, 1.5))
# Take 100000 random variable samples from this uniform mass distribution.
mass1_samples = mass1_distribution.rvs(size=1000000)

# Draw another distribution that is Gaussian between 0.5 and 1.5 solar masses
# with a mean of 1.2 solar masses and a standard deviation of 0.15 solar
# masses. Gaussian takes the variance as an input so square the standard
# deviation.
variance = 0.15*0.15
mass2_gaussian = distributions.Gaussian(mass2=(0.5, 1.5), mass2_mean=1.2,
                                        mass2_var=variance)

# Take 100000 random variable samples from this gaussian mass distribution.
mass2_samples = mass2_gaussian.rvs(size=1000000)

# We can make pairs of distributions together, instead of apart.
two_mass_distributions = distributions.Uniform(mass3=(1.6, 3.0),
                                               mass4=(1.6, 3.0))
two_mass_samples = two_mass_distributions.rvs(size=1000000)

# Choose 50 bins for the histogram subplots.
n_bins = 50

# Plot histograms of samples in subplots
fig, axes = plt.subplots(nrows=2, ncols=2)
ax0, ax1, ax2, ax3, = axes.flat

ax0.hist(mass1_samples['mass1'], bins = n_bins)
ax1.hist(mass2_samples['mass2'], bins = n_bins)
ax2.hist(two_mass_samples['mass3'], bins = n_bins)
ax3.hist(two_mass_samples['mass4'], bins = n_bins)

ax0.set_title('Mass 1 samples')
ax1.set_title('Mass 2 samples')
ax2.set_title('Mass 3 samples')
ax3.set_title('Mass 4 samples')

plt.tight_layout()
plt.show()
