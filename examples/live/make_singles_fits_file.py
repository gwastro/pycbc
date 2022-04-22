"""
This script will make a valid singles-fits file for use in the
pycbc_live CI tests. It doesn't have much physical meaning,
but will give broadly representative numbers for singles.
"""
import h5py
import numpy as np

f = h5py.File('single_trigger_fits.hdf','w')

# Some numbers to design the output
# These are loosely based on the O3a trigger fits file
n_days = 30
n_bins = 5
max_duration = 150.
min_duration = 6.
duty_cycle = 0.7
alpha = 4.
daily_counts_per_bin = 500.

# Set attributes to be something appropriate
f.attrs['start_date'] = 0
f.attrs['end_date'] = n_days
f.attrs['fit_threshold'] = 6.5

f['bins_edges'] = np.logspace(np.log10(min_duration),
                              np.log10(max_duration),
                              n_bins + 1)

# Create the fits groups for each ifo
for ifo in ['H1', 'L1', 'V1']:
    f.create_group(ifo)
    ifo_group = f[ifo]
    # Group attributes are set sensibly
    ifo_group.attrs['live_time'] = np.round(n_days * 0.7 * 86400.)
    ifo_group.attrs['mean_alpha'] = alpha
    ifo_group.attrs['total_counts'] = daily_counts_per_bin * n_days * n_bins

    # Make daily fits datasets, this doesn't atually matter for the test, but
    # means that the file format matches in case of future requirements
    ifo_group.create_group('daily_fits')
    daily_group = ifo_group['daily_fits']
    for bin_no in range(n_bins):
        key = "bin_%d" % bin_no
        daily_group.create_group(key)
        bin_group = daily_group[key]
        bin_group['date'] = np.arange(n_days)
        bin_group['counts'] = np.array([daily_counts_per_bin] * n_days)
        bin_group['fit_coeff'] = np.array([alpha] * n_days)

    # Again, not used for current tests, but keeps file format the same
    ifo_group.create_group('fixed')
    fixed_group = ifo_group['fixed']
    fixed_group['fit_coeff'] = np.array([0.] * n_bins)
    fixed_group['counts'] = np.array([1.] * n_bins)

    # As the daily fit are all the same, the mean and 95th %ile values
    # are the same
    for k in ['conservative', 'mean']:
        ifo_group.create_group(k)
        k_group = ifo_group[k]
        k_group['fit_coeff'] = np.array([alpha] * n_bins)
        k_group['counts'] = np.array([daily_counts_per_bin * n_days] * n_bins)
