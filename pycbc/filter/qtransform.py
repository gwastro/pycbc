# Copyright (C) 2017  Hunter A. Gabbard, Andrew Lundgren, Duncan Macleod
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#

"""
This module retrives a timeseries and then calculates the q-transform of that time series

Example
-------
    $ python q-transform.py -s 4096 -u test -o /Users/pycbc_qtransform

"""


from numpy import pi, ceil, log, exp
import numpy as np
from pycbc.strain  import next_power_of_2
from pycbc.types.timeseries import FrequencySeries, TimeSeries
from scipy.interpolate import (interp2d)
from pycbc.fft import ifft
from pycbc.types import zeros
from numpy import fft as npfft
import os


def plotter(interp, out_dir, now, frange, tres, fres):
    """Plotting mechanism for pycbc spectrograms

    Parameters
    ----------
    interp:
        2D interpolation function
    out_dir:
        path to output directory
    now:
        unique label for output directory
    frange:
        upper and lower bounds on frequency range plotted
    tres:
        desired time resolution
    fres:
        desired frequency resolution

    Returns
    -------
    plt:
        matplotlib spectrogram figure

    """

    from matplotlib import use
    use('Agg')
    from matplotlib import pyplot as plt

    # create directory where figure will be saved
    os.makedirs('%s/run_%s' % (out_dir,now))  # Fail early if the dir already exists

    # plot a spectrogram of the q-plane with the loudest normalized tile energy
    print 'plotting ...'
    # intitialize variables
    xarr=np.linspace(0, 1, 1. / tres)
    yarr=np.arange(int(frange[0]), int(frange[1]), fres)
    z = interp(xarr,yarr)

    # pick the desired colormap
    cmap = plt.get_cmap('viridis')

    plt.plot()

    p1 = plt.pcolormesh(xarr,
                        yarr,
                        z, cmap=cmap, norm=None)
    plt.colorbar(p1)
    plt.title('Pycbc q-transform')

    plt.savefig('%s/run_%s/spec.png' % (out_dir,now))

    return plt

def qplane(qplane_tile_dict, fseries, frange, normalized=True, tres=0.001, fres=1.):
    """Performs q-transform on each tile for each q-plane and selects
       tile with the maximum normalized energy. Q-transform is then
       interpolated to a desired frequency and time resolution.

    Parameters
    ----------
    qplane_tile_dict:
        Dictionary containing a list of q-tile tupples for each q-plane
    fseries: 'pycbc FrequencySeries'
        frequency-series data set
    normalized: 'bool'
        normalize energy time series?
    frange:
        upper and lower bounds on frequency range
    tres:
        desired time resolution
    fres:
        desired frequency resolution

    Returns
    -------
    out:
        2D interpolated q-transform
    interp:
        2D interpolation function

    """
    # store q-transforms of each tile in a dict
    qplane_qtrans_dict = {}
    dur = 1.0 / fseries.delta_f

    # check for sampling rate
    sampling = (len(fseries) - 1) * 2 * fseries.delta_f

    max_energy = []
    for i, key in enumerate(qplane_tile_dict):
        energies_lst=[]
        for tile in qplane_tile_dict[key]:
            energies = qtransform(fseries, tile[1], tile[0])
            if normalized:
                energies = energies[0]
            else:
                energies = energies[1]
            energies_lst.append(energies.numpy())
            maxen = float(energies.max())
            if i == 0:
                max_energy.append(maxen)
                max_energy.append(tile)
                max_energy.append(key)
            elif maxen > max_energy[0]:
                max_energy[0] = maxen
                max_energy[1] = tile
                max_energy[2] = key
                max_energy[3] = energies
        qplane_qtrans_dict[key] = np.array(energies_lst)

    # record peak q calculate above and q-transform output for peak q
    result = qplane_qtrans_dict[max_energy[2]]

    # then interpolate the spectrogram to increase the frequency resolution
    if fres is None:  # unless user tells us not to
        return result
    else:
        # initialize some variables
        frequencies = []

        for i in qplane_tile_dict[max_energy[2]]:
            frequencies.append(i[0])

        # 2-D interpolation
        time_array = np.linspace(0, dur, int(dur * sampling))
        interp = interp2d(time_array, frequencies, result)
 
    out = interp(np.arange(0, dur, tres), np.arange(int(frange[0]), int(frange[1]), fres))
    return out, interp

def qtiling(fseries, qrange, frange, mismatch):
    """Iterable constructor of QTile tuples

    Parameters
    ----------
    fseries: 'pycbc FrequencySeries'
        frequency-series data set
    qrange:
        upper and lower bounds of q range
    frange:
        upper and lower bounds of frequency range
    mismatch:
        percentage of desired fractional mismatch

    Returns
    -------
    qplane_tile_dict: 'dict'
        dictionary containing Q-tile tuples for a set of Q-planes
    frange: 'list'
        upper and lower bounds on frequency range
    """

    deltam = deltam_f(mismatch)
    qrange = (float(qrange[0]), float(qrange[1]))
    frange = [float(frange[0]), float(frange[1])]
    dur = 1.0 / fseries.delta_f
    qplane_tile_dict = {}

    # check for sampling rate
    sampling = (len(fseries) - 1) * 2 * fseries.delta_f

    qs = list(_iter_qs(qrange, deltam))
    if frange[0] == 0:  # set non-zero lower frequency
        frange[0] = 50 * max(qs) / (2 * pi * dur)
    if np.isinf(frange[1]):  # set non-infinite upper frequency
        frange[1] = sampling / 2 / (1 + 11**(1/2.) / min(qs))

    #lets now define the whole tiling (e.g. choosing all tiling in planes)
    for q in qs:
        qtilefreq = np.array(list(_iter_frequencies(q, frange, mismatch, dur)))
        qlst = np.empty(len(qtilefreq), dtype=float)
        qlst.fill(q)
        qtiles_array = np.vstack((qtilefreq, qlst)).T
        qplane_tiles_list = list(map(tuple, qtiles_array))
        qplane_tile_dict[q] = qplane_tiles_list

    return qplane_tile_dict, frange

def deltam_f(mismatch):
    """Fractional mismatch between neighbouring tiles

    Parameters
    ----------
    mismatch: 'float'
        percentage of desired fractional mismatch

    Returns
    -------
    :type: 'float'
    """
    return 2 * (mismatch / 3.) ** (1/2.)


def _iter_qs(qrange, deltam):
    """Iterate over the Q values
   
    Parameters
    ----------
    qrange:
        upper and lower bounds of q range
    deltam:
        Fractional mismatch between neighbouring tiles

    Returns
    -------
    Q-value:
        Q value for Q-tile
    """

    # work out how many Qs we need
    cumum = log(qrange[1] / qrange[0]) / 2**(1/2.)
    nplanes = int(max(ceil(cumum / deltam), 1))
    dq = cumum / nplanes
    for i in xrange(nplanes):
        yield qrange[0] * exp(2**(1/2.) * dq * (i + .5))
    raise StopIteration()

def _iter_frequencies(q, frange, mismatch, dur):
    """Iterate over the frequencies of this 'QPlane'

    Parameters
    ----------
    q:
        q value
    frange: 'list'
        upper and lower bounds of frequency range
    mismatch:
        percentage of desired fractional mismatch
    dur:
        duration of timeseries in seconds

    Returns
    -------
    frequencies:
        Q-Tile frequency
    """
    # work out how many frequencies we need
    minf, maxf = frange
    fcum_mismatch = log(maxf / minf) * (2 + q**2)**(1/2.) / 2.
    nfreq = int(max(1, ceil(fcum_mismatch / deltam_f(mismatch))))
    fstep = fcum_mismatch / nfreq
    fstepmin = 1. / dur
    # for each frequency, yield a QTile
    for i in xrange(nfreq):
        yield (minf *
               exp(2 / (2 + q**2)**(1/2.) * (i + .5) * fstep) //
               fstepmin * fstepmin)
    raise StopIteration()

def qtransform(fseries, Q, f0):
    """Calculate the energy 'TimeSeries' for the given fseries

    Parameters
    ----------
    fseries: 'pycbc FrequencySeries'
        frequency-series data set
    Q:
        q value
    f0:
        central frequency

    Returns
    -------
    norm_energy: '~pycbc.types.aligned.ArrayWithAligned'
        A 'TimeSeries' of the normalized energy from the Q-transform of
        this tile against the data.
    cenergy: '~pycbc.types.aligned.ArrayWithAligned'
        A 'TimeSeries' of the complex energy from the Q-transform of 
        this tile against the data.
    """
    # q-transform data for each (Q, frequency) tile
    # initialize parameters
    qprime = Q / 11**(1/2.) # ... self.qprime
    dur = 1.0 / fseries.delta_f

    # check for sampling rate
    sampling = (len(fseries) - 1) * 2 * fseries.delta_f

    # choice of output sampling rate
    output_sampling = sampling # Can lower this to highest bandwidth
    output_samples = int(dur * output_sampling)

    # window fft
    window_size = 2 * int(f0 / qprime * dur) + 1

    # get start and end indices
    start = int((f0 - (f0 / qprime)) * dur)
    end = int(start + window_size)

    # apply window to fft
    # normalize and generate bi-square window
    norm = np.sqrt(315. * qprime / (128. * f0))
    windowed = fseries[start:end].numpy() * (bisquare(window_size) * norm)

    # pad data, move negative frequencies to the end, and IFFT
    padded = np.pad(windowed, padding(window_size, output_samples), mode='constant')
    wenergy = npfft.ifftshift(padded)

    # return a 'TimeSeries'
    wenergy = FrequencySeries(wenergy, delta_f=1./dur)
    cenergy = TimeSeries(zeros(output_samples, dtype=np.complex128),
                            delta_t=1./sampling)
    ifft(wenergy, cenergy)
    energy = cenergy.squared_norm()
    medianenergy = np.median(energy.numpy())
    norm_energy = energy / float(medianenergy)
    return norm_energy, cenergy

def padding(window_size, desired_size):
    """The (left, right) padding required for the IFFT

    Parameters
    ----------
    window_size: int
        Size of window
    desired_size: int
        Desired size of window

    Returns
    -------
    tuple
       Number of values padded to the edges of each axis. 
 
    """
    pad = desired_size - window_size
    return (int(pad/2.), int((pad + 1)/2.))


def _get_indices(window_size):
    """ Windows indices for fft
    
    Parameters
    ---------
    window_size: int
        size of window

    Returns
    -------
    numpy.ndarray
        Window indices for fft using total duration of segment

    """
    half = int((window_size - 1) / 2.)
    return np.arange(-half, half + 1)

def bisquare(size):
    """Generate the bi-square window for this row
 
    Parameters
    ---------
    size: int
        size of window

    Returns
    -------
    window : numpy.ndarray
    """
    # dimensionless frequencies
    xfrequencies = np.linspace(-1., 1., size)

    # ported from https://github.com/gwpy/gwpy/blob/master/gwpy/signal/qtransform.py
    return (1 - xfrequencies ** 2) ** 2

def n_tiles(dur, f0, Q):
    """The number of tiles in this row 
    
    Parameters
    ----------
    dur: int
        Duration of timeseries in seconds
    f0: int
        Central frequency
    Q: int
        q value

    Returns
    -------
    :type: 'int'
    """

    # ported from https://github.com/gwpy/gwpy/blob/master/gwpy/signal/qtransform.py
    tcum_mismatch = dur * 2 * pi * f0 / Q  
    return next_power_of_2(tcum_mismatch / deltam())

def deltam():
    """Fractional mismatch between neighbouring tiles

    Parameters
    ----------
    None

    Returns
    -------
    :type: 'float'
    """
    mismatch = 0.2
    return 2 * (mismatch / 3.) ** (1/2.)
