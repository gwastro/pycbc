import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt

import pycbc
from pycbc.results import ifo_color


def generate_asd_plot(psddict, output_filename):
    """
    Generate an ASD plot as used for upload to GraceDB
    Parameters
    ----------

    psddict:
        A dictionary keyed on ifo containing the PSDs as
        FrequencySeries objects

    output_filename:
        The filename for the plot to be saved to
    """
    asd_fig, asd_ax = plt.subplots(1)
    for ifo in sorted(psddict.keys()):
        curr_psd = psddict[ifo]
        # Can't plot log(0) so start from point 1
        asd_ax.loglog(curr_psd.sample_frequencies[1:],
                      curr_psd[1:] ** 0.5 / pycbc.DYN_RANGE_FAC,
                      c=ifo_color(ifo), label=ifo)
    asd_ax.legend()
    asd_ax.set_xlim([10, 1300])
    asd_ax.set_ylim([3E-24, 1E-20])
    asd_ax.set_xlabel('Frequency (Hz)')
    asd_ax.set_ylabel('ASD')
    asd_fig.savefig(output_filename)
