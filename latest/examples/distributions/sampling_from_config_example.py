import numpy as np
import matplotlib.pyplot as plt
from pycbc.distributions.utils import draw_samples_from_config


# A path to the .ini file.
CONFIG_PATH = "./pycbc_bbh_prior.ini"
random_seed = np.random.randint(low=0, high=2**32-1)

# Draw a single sample.
sample = draw_samples_from_config(
            path=CONFIG_PATH, num=1, seed=random_seed)

# Print all parameters.
print(sample.fieldnames)
print(sample)
# Print a certain parameter, for example 'mass1'.
print(sample[0]['mass1'])

# Draw 1000000 samples, and select all values of a certain parameter.
n_bins = 50
samples = draw_samples_from_config(
            path=CONFIG_PATH, num=1000000, seed=random_seed)

fig, axes = plt.subplots(nrows=3, ncols=2)
ax1, ax2, ax3, ax4, ax5, ax6 = axes.flat

ax1.hist(samples[:]['srcmass1'], bins=n_bins)
ax2.hist(samples[:]['srcmass2'], bins=n_bins)
ax3.hist(samples[:]['comoving_volume'], bins=n_bins)
ax4.hist(samples[:]['redshift'], bins=n_bins)
ax5.hist(samples[:]['distance'], bins=n_bins)
ax6.hist(samples[:]['mass1'], bins=n_bins)

ax1.set_title('srcmass1')
ax2.set_title('srcmass2')
ax3.set_title('comoving_volume')
ax4.set_title('redshift')
ax5.set_title('distance')
ax6.set_title('mass1 or mass2')

plt.tight_layout()
plt.show()
