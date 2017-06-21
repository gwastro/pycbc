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
import h5py
from pycbc.filter import highpass_fir, matched_filter
from pycbc.psd import welch, interpolate

def pycbc_insp_main(segments, filename='q_info.hdf'):
    """Main function for pycbc_inspiral implementation of qtransform.py
    Parameters
    ----------
    segments:
        list of pycbc frequency series segments
    filename:
        name of output file for qtransform info to be stored

    """

    comb_q_dict = {}
    for s_num, stilde in enumerate(segments):
        # getting q-tiles for segments
        Qbase, q_frange, q_data = pycbc_insp_tiling(stilde, stilde.psd)

        # getting q-plane for segment
        Qplane, interp, qs_time, qe_time = qplane(Qbase, q_data, q_frange, fres=None, seg=stilde)
        
        del q_frange, q_data

        # write q info to an hdf file
        Qbase_tmp = {}
        Qplane_tmp = {}
        print 'start time: %s, end time: %s' % (qs_time, qe_time)
        if not 'qtiles' in comb_q_dict:
            Qbase_tmp['seg_%s-%s' % (str(qs_time),str(qe_time))] = Qbase
            comb_q_dict['qtiles'] = Qbase_tmp
        if not 'qplanes' in comb_q_dict:
            Qplane_tmp['seg_%s-%s' % (str(qs_time),str(qe_time))] = Qplane
            comb_q_dict['qplanes'] = Qplane_tmp
        else:
            comb_q_dict['qtiles']['seg_%s-%s' % (str(qs_time),str(qe_time))] = Qbase
            comb_q_dict['qplanes']['seg_%s-%s' % (str(qs_time),str(qe_time))] = Qplane

    # save qtransform info to hdf5 file
    #path = '/'
    #recursively_save_dict_contents_to_group(out_vals, path, comb_q_dict)

    return comb_q_dict

def pycbc_insp_tiling(seg, frange=(0,1024), qrange=(4,64)):
    """Iterable constructor of QTile tuples
    Parameters
    ----------
    seg:
        fft'd analysis chunck
    frange:
        frequency range
    qrange:
        q range

    Returns
    -------
    Qbase: 'dict'
        dictionary containing Q-tile tuples for a set of Q-planes
    frange: 'list'
        upper and lower bounds on frequency range
    data:
        whitened strain of analysis chunk

    """

    # assume segment is fft'd
    # check for sampling rate
    sampling = (len(seg) - 1) * 2 * seg.delta_f
    sampling = 64.
    mismatch = 0.2
    qrange=(4,64)
    frange=(0,np.inf)

    valid_time = seg.to_timeseries()

    # highpass and suppress frequencies below 20Hz
    data = highpass_fir(valid_time, 20, 8)
    dur = data.duration

    # calculate the noise spectrum
    seg_psd = interpolate(welch(data), 1.0 / dur)

    # retrieve whitened strain
    white_strain = (data.to_frequencyseries() / seg_psd ** 0.5 * seg_psd.delta_f)            
    data = white_strain
 
    # perform Q-tiling
    Qbase, frange = qtiling(data, qrange, frange, mismatch)

    return Qbase, frange, data

def save_dict_to_hdf5(dic, filename):
    """
    Parameters
    ----------
    dic:
        python dictionary to be converted to hdf5 format
    filename:
        desired name of hdf5 file
    """
    with h5py.File(filename, 'w') as h5file:
        recursively_save_dict_contents_to_group(h5file, '/', dic)

def recursively_save_dict_contents_to_group(h5file, path, dic):
    """
    Parameters
    ----------
    h5file:
        h5py file to be written to
    path:
        path within h5py file to saved dictionary
    dic:
        python dictionary to be converted to hdf5 format
    """
    for key, item in dic.items():
        if isinstance(item, (np.ndarray, np.int64, np.float64, str, bytes, tuple, list)):
            h5file[path + str(key)] = item
        elif isinstance(item, dict):
            recursively_save_dict_contents_to_group(h5file, path + key + '/', item)
        else:
            raise ValueError('Cannot save %s type'%type(item))

def load_dict_from_hdf5(filename):
    """
    Parameters
    ----------
    filename:
        name of h5py file to be loaded

    Returns
    -------
    recursively_load_dict_contents_from_group():
        python dictionary object
    """
    with h5py.File(filename, 'r') as h5file:
        return recursively_load_dict_contents_from_group(h5file, '/')

def recursively_load_dict_contents_from_group(h5file, path):
    """
    Parameters
    ----------
    h5file:
        h5py file to be loaded from
    path:
        path within h5py file to items
    """
    ans = {}
    for key, item in h5file[path].items():
        if isinstance(item, h5py._hl.dataset.Dataset):
            ans[key] = item.value

    

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

def qplane(qplane_tile_dict, fseries, frange, normalized=True, tres=1., fres=1., seg=None):
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
    dur = fseries.to_timeseries().duration

    # check for sampling rate
    sampling = (len(fseries) - 1) * 2 * fseries.delta_f
    #sampling = 64.

    max_energy = []
    for i, key in enumerate(qplane_tile_dict):
        energies_lst=[]
        for tile in qplane_tile_dict[key]:
            energies = qtransform(fseries, tile[1], tile[0])
            if normalized:
                energies = energies[0]
            else:
                energies = energies[1]
            energies_lst.append(energies)
            if i == 0:
                max_energy.append(max(energies))
                max_energy.append(tile)
                max_energy.append(key)
            elif max(energies) > max_energy[0]:
                max_energy[0] = max(energies)
                max_energy[1] = tile
                max_energy[2] = key
                max_energy[3] = energies 
        qplane_qtrans_dict[key] = energies_lst

    # record q-transform output for peak q
    result = qplane_qtrans_dict[max_energy[2]]
    qtile_max = qplane_tile_dict[key]

    # then interpolate the spectrogram to increase the frequency resolution
    if fres is None:  # unless user tells us not to
        if not seg:
            return result, qtile_max
        elif seg:
            #result[0].start_time
            s_time = seg.epoch + (seg.analyze.start*(1. / sampling))
            e_time = seg.epoch + (seg.analyze.stop*(1. / sampling))
            result = np.array(result)
            s_idx = int(((len(result[:,0]) / dur) / sampling) * seg.analyze.start)
            e_idx = int(((len(result[:,0]) / dur) / sampling) * seg.analyze.stop)
            result = result[:,s_idx:e_idx]
            return result, qtile_max, s_time, e_time
    else:
        # initialize some variables
        frequencies = []

        for i in qplane_tile_dict[max_energy[2]]:
            frequencies.append(i[0])

        # 2-D interpolation
        time_array = np.linspace(-(dur / 2.),(dur / 2.),int(dur * sampling))
        interp = interp2d(time_array, frequencies, np.array(result))

    out = interp(np.linspace(0,1, 1. / tres), np.arange(int(frange[0]), int(frange[1]), fres))

    # if segment, downsample to valid time
    if not seg:
        seg = None
    elif seg:
        s_idx = int((((1. / tres) / dur) / sampling) * seg.analyze.start)
        e_idx = int((((1. / tres) / dur) / sampling) * seg.analyze.stop)
        out = out[s_idx:e_idx,:]

    return out, interp

def qtiling(fseries, qrange, frange, mismatch=0.2):
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
    dur = fseries.to_timeseries().duration
    qplane_tile_dict = {}

    # check for sampling rate
    sampling = (len(fseries) - 1) * 2 * fseries.delta_f
    #sampling = 64.

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
        qtiles_array = np.vstack((qtilefreq,qlst)).T
        qplane_tiles_list = list(map(tuple,qtiles_array))
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
    dur = fseries.to_timeseries().duration

    # check for sampling rate
    sampling = (len(fseries) - 1) * 2 * fseries.delta_f
    #sampling = 64.

    # window fft
    window_size = 2 * int(f0 / qprime * dur) + 1

    # get start and end indices
    start = int((f0 - (f0 / qprime)) * dur)
    end = int(start + window_size)

    # apply window to fft
    # normalize and generate bi-square window
    norm = np.sqrt(315. * qprime / (128. * f0))
    windowed = fseries[start:end] * (bisquare(window_size) * norm)

    # choice of output sampling rate
    output_sampling = sampling # Can lower this to highest bandwidth
    output_samples = dur * output_sampling

    # pad data, move negative frequencies to the end, and IFFT
    padded = np.pad(windowed, padding(window_size, output_samples), mode='constant')
    wenergy = npfft.ifftshift(padded)

    # return a 'TimeSeries'
    wenergy = FrequencySeries(wenergy, delta_f=1./dur)
    tdenergy = TimeSeries(zeros(int(output_samples), dtype=np.complex128),
                            delta_t=1./sampling)
    ifft(wenergy, tdenergy)
    cenergy = TimeSeries(tdenergy,
                         delta_t=tdenergy.delta_t, copy=False)
    energy = type(cenergy)(
        cenergy.real() ** 2. + cenergy.imag() ** 2.,
        delta_t=1, copy=False)
    medianenergy = np.median(energy)
    norm_energy = energy / medianenergy
 
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
