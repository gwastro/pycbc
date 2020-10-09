# Copyright (C) 2020  Daniel Finstad
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
"""This module provides model classes and functions for implementing
a relative binning likelihood for parameter estimation.
"""


import logging
import numpy
from scipy.interpolate import interp1d
from scipy import special

from pycbc.waveform import get_fd_waveform_sequence
from pycbc.detector import Detector
from pycbc.types import Array

from .gaussian_noise import BaseGaussianNoise


def setup_bins(f_full, f_lo, f_hi, chi=1.0, eps=0.5, gammas=None):
    """Construct frequency bins for use in a relative likelihood
    model. For details, see [Barak, Dai & Venumadhav 2018].

    Parameters
    ----------
    f_full : array
        The full resolution array of frequencies being used in the analysis.
    f_lo : float
        The starting frequency used in matched filtering. This will be the
        left edge of the first frequency bin.
    f_hi : float
        The ending frequency used in matched filtering. This will be the right
        edge of the last frequency bin.
    chi : float, optional
        Tunable parameter, see [Barak, Dai & Venumadhav 2018]
    eps : float, optional
        Tunable parameter, see [Barak, Dai & Venumadhav 2018]. Lower values
        result in larger number of bins.
    gammas : array, optional
        Frequency powerlaw indices to be used in computing bins.

    Returns
    -------
    nbin : int
        Number of bins.
    fbin : numpy.array of floats
        Bin edge frequencies.
    fbin_ind : numpy.array of ints
        Indices of bin edges in full frequency array.
    """
    f = numpy.linspace(f_lo, f_hi, 10000)
    # f^ga power law index
    ga = gammas if gammas is not None else \
        numpy.array([-5./3, -2./3, 1., 5./3, 7./3])
    logging.info("Using powerlaw indices: %s", ga)
    dalp = chi * 2.0 * numpy.pi / numpy.absolute((f_lo ** ga) - (f_hi ** ga))
    dphi = numpy.sum(numpy.array([numpy.sign(g) * d * (f ** g) for
                                  g, d in zip(ga, dalp)]), axis=0)
    dphi_diff = dphi - dphi[0]
    # now construct frequency bins
    nbin = int(dphi_diff[-1] / eps)
    dphi2f = interp1d(dphi_diff, f, kind='slinear', bounds_error=False,
                      fill_value=0.0)
    dphi_grid = numpy.linspace(dphi_diff[0], dphi_diff[-1], nbin+1)
    # frequency grid points
    fbin = dphi2f(dphi_grid)
    # indices of frequency grid points in the FFT array
    fbin_ind = numpy.unique([numpy.argmin(numpy.absolute(f_full - ff)) for
                             ff in fbin])
    # make sure grid points are precise
    fbin = numpy.array([f_full[i] for i in fbin_ind])
    nbin = len(fbin)
    return nbin, fbin, fbin_ind


class Relative(BaseGaussianNoise):
    r"""Model that assumes the likelihood in a region around the peak
    is slowly varying such that a linear approximation can be made, and
    likelihoods can be calculated at a coarser frequency resolution. For
    more details on the implementation, see https://arxiv.org/abs/1806.08792.

    This model requires the use of a fiducial waveform whose parameters are
    near the peak of the likelihood. The fiducial waveform and all template
    waveforms used in likelihood calculation are currently generated using
    the SPAtmplt approximant.

    For more details on initialization parameters and definition of terms, see
    :py:class:`BaseGaussianNoise`.

    Parameters
    ----------
    variable_params : (tuple of) string(s)
        A tuple of parameter names that will be varied.
    data : dict
        A dictionary of data, in which the keys are the detector names and the
        values are the data (assumed to be unwhitened). All data must have the
        same frequency resolution.
    low_frequency_cutoff : dict
        A dictionary of starting frequencies, in which the keys are the
        detector names and the values are the starting frequencies for the
        respective detectors to be used for computing inner products.
    figucial_params : dict
        A dictionary of waveform parameters to be used for generating the
        fiducial waveform. Keys must be parameter names in the form
        'PARAM_ref' where PARAM is a recognized extrinsic parameter or
        an intrinsic parameter compatible with the chosen approximant.
    gammas : array of floats, optional
        Frequency powerlaw indices to be used in computing frequency bins.
    epsilon : float, optional
        Tuning parameter used in calculating the frequency bins. Lower values
        will result in higher resolution and more bins.
    \**kwargs :
        All other keyword arguments are passed to
        :py:class:`BaseGaussianNoise`.
    """
    name = "relative"

    def __init__(self, variable_params, data, low_frequency_cutoff,
                 fiducial_params=None, gammas=None, epsilon=0.5, **kwargs):
        super(Relative, self).__init__(
            variable_params, data, low_frequency_cutoff, **kwargs)
        # check that all of the frequency cutoffs are the same
        # FIXME: this can probably be loosened at some point
        kmins = list(self.kmin.values())
        kmaxs = list(self.kmax.values())
        if any(kk != kmins[0] for kk in kmins):
            raise ValueError("All lower frequency cutoffs must be the same")
        if any(kk != kmaxs[0] for kk in kmaxs):
            raise ValueError("All high frequency cutoffs must be the same")
        # store data and frequencies
        d0 = list(self.data.values())[0]
        self.f = numpy.array(d0.sample_frequencies)
        self.df = d0.delta_f
        self.end_time = float(d0.end_time)
        self.det = {ifo: Detector(ifo) for ifo in self.data}
        self.epsilon = float(epsilon)
        # store data and psds as arrays for faster computation
        self.comp_data = {ifo: d.numpy() for ifo, d in self.data.items()}
        self.comp_psds = {ifo: p.numpy() for ifo, p in self.psds.items()}
        # store fiducial waveform params
        self.fid_params = fiducial_params

        # get detector-specific arrival times relative to end of data
        dt = {ifo:
              self.det[ifo].time_delay_from_earth_center(
                  self.fid_params['ra'], self.fid_params['dec'],
                  self.fid_params['tc'])
              for ifo in self.data}
        self.ta = {ifo: self.fid_params['tc'] + dt[ifo] - self.end_time
                   for ifo in self.data}

        # generate fiducial waveform
        f_lo = kmins[0] * self.df
        f_hi = kmaxs[0] * self.df
        logging.info("Generating fiducial waveform from %s to %s Hz",
                     f_lo, f_hi)
        # prune low frequency samples to avoid waveform errors
        nbelow = sum(self.f < 10)
        fpoints = Array(self.f.astype(numpy.float64))[nbelow:]
        approx = self.static_params['approximant']
        fid_hp, fid_hc = get_fd_waveform_sequence(approximant=approx,
                                                  sample_points=fpoints,
                                                  **self.fid_params)
        # check for zeros at high frequencies
        numzeros = list(fid_hp[::-1] != 0j).index(True)
        n_above_fhi = (len(self.f) - 1) - kmaxs[0]
        # make sure only nonzero samples are included in bins
        if numzeros > n_above_fhi:
            nremove = numzeros - n_above_fhi
            new_kmax = kmaxs[0] - nremove
            f_hi = new_kmax * self.df
            logging.info("WARNING! Fiducial waveform terminates below "
                         "high-frequency-cutoff, final bin frequency "
                         "will be %s Hz", f_hi)
        self.h00 = {}
        for ifo in self.data:
            # make copy of fiducial wfs, adding back in low frequencies
            hp0 = numpy.concatenate([[0j] * nbelow, fid_hp.copy()])
            hc0 = numpy.concatenate([[0j] * nbelow, fid_hc.copy()])
            fp, fc = self.det[ifo].antenna_pattern(
                self.fid_params['ra'], self.fid_params['dec'],
                self.fid_params['polarization'], self.fid_params['tc'])
            tshift = numpy.exp(-2.0j * numpy.pi * self.f * self.ta[ifo])
            self.h00[ifo] = numpy.array(hp0 * fp + hc0 * fc) * tshift

        # compute frequency bins
        logging.info("Computing frequency bins")
        nbin, fbin, fbin_ind = setup_bins(f_full=self.f, f_lo=f_lo, f_hi=f_hi,
                                          gammas=gammas, eps=self.epsilon)
        logging.info("Using %s bins for this model", nbin)
        # store bins and edges in sample and frequency space
        self.edges = fbin_ind
        self.fedges = numpy.array(fbin).astype(numpy.float64)
        self.bins = numpy.array([(self.edges[i], self.edges[i+1]) for
                                 i in range(len(self.edges) - 1)])
        self.fbins = numpy.array([(fbin[i], fbin[i+1]) for
                                  i in range(len(fbin) - 1)])
        # store low res copy of fiducial waveform
        self.h00_sparse = {ifo: self.h00[ifo].copy().take(self.edges) for ifo
                           in self.h00}

        # compute summary data
        logging.info("Calculating summary data at frequency resolution %s Hz",
                     self.df)
        self.sdat = self.summary_data()

    def summary_data(self):
        """Compute summary data bin coefficients encoding linear
        approximation to full resolution likelihood.

        Returns
        -------
        dict of dicts
            Dictionary keyed by detector name, whose values are dictionaries
            containing bin coefficients a0, b0, a1, b1, for each frequency
            bin.
        """
        # calculate coefficients
        sdat = {}
        for ifo in self.data:
            hd = numpy.conjugate(self.comp_data[ifo]) * self.h00[ifo]
            hd /= self.comp_psds[ifo]
            hh = (numpy.absolute(self.h00[ifo]) ** 2.0) / self.comp_psds[ifo]
            # constant terms
            a0 = numpy.array([4. * self.df * numpy.sum(hd[l:h]) for
                              l, h in self.bins])
            b0 = numpy.array([4. * self.df * numpy.sum(hh[l:h]) for
                              l, h in self.bins])
            # linear terms
            bin_lefts = [fl for fl, fh in self.fbins]
            a1 = numpy.array([4. * self.df
                              * numpy.sum(hd[l:h] * (self.f[l:h] - bl)) for
                              (l, h), bl in zip(self.bins, bin_lefts)])
            b1 = numpy.array([4. * self.df
                              * numpy.sum(hh[l:h] * (self.f[l:h] - bl)) for
                              (l, h), bl in zip(self.bins, bin_lefts)])

            sdat[ifo] = {'a0': a0, 'a1': a1,
                         'b0': b0, 'b1': b1}
        return sdat

    def _loglr(self):
        r"""Computes the log likelihood ratio,

        .. math::

            \log \mathcal{L}(\Theta) = \sum_i
                \left<h_i(\Theta)|d_i\right> -
                \frac{1}{2}\left<h_i(\Theta)|h_i(\Theta)\right>,

        at the current parameter values :math:`\Theta`.

        Returns
        -------
        float
            The value of the log likelihood ratio.
        """
        # get model params
        p = self.current_params.copy()
        p.update(self.static_params)

        hh = 0.
        hd = 0j
        for ifo in self.data:
            # get detector antenna pattern
            fp, fc = self.det[ifo].antenna_pattern(p['ra'], p['dec'],
                                                   p['polarization'],
                                                   p['tc'])
            # get timeshift relative to end of data
            dt = self.det[ifo].time_delay_from_earth_center(p['ra'], p['dec'],
                                                            p['tc'])
            dtc = p['tc'] + dt - self.end_time
            tshift = numpy.exp(-2.0j * numpy.pi * self.fedges * dtc)
            # generate template and calculate waveform ratio
            hp, hc = get_fd_waveform_sequence(sample_points=Array(self.fedges),
                                              **p)
            htilde = numpy.array(fp * hp + fc * hc) * tshift
            r = (htilde / self.h00_sparse[ifo]).astype(numpy.complex128)
            r0 = r[:-1]
            r1 = (r[1:] - r[:-1]) / (self.fedges[1:] - self.fedges[:-1])

            # <h, d> is sum over bins of A0r0 + A1r1
            hd += numpy.sum(self.sdat[ifo]['a0'] * r0
                            + self.sdat[ifo]['a1'] * r1)
            # <h, h> is sum over bins of B0|r0|^2 + 2B1Re(r1r0*)
            hh += numpy.sum(self.sdat[ifo]['b0'] * numpy.absolute(r0) ** 2.
                            + 2. * self.sdat[ifo]['b1']
                            * (r1 * numpy.conjugate(r0)).real)
        hd = abs(hd)
        llr = numpy.log(special.i0e(hd)) + hd - 0.5 * hh
        return float(llr)

    def write_metadata(self, fp):
        """Adds writing the fiducial parameters and epsilon to file's attrs.

        Parameters
        ----------
        fp : pycbc.inference.io.BaseInferenceFile instance
            The inference file to write to.
        """
        super(Relative, self).write_metadata(fp)
        fp.attrs['epsilon'] = self.epsilon
        for p, v in self.fid_params.items():
            fp.attrs['{}_ref'.format(p)] = v

    @staticmethod
    def extra_args_from_config(cp, section, skip_args=None, dtypes=None):
        """Adds reading fiducial waveform parameters from config file.
        """
        # add fiducial params to skip list
        skip_args += [option for option in cp.options(section) if
                      option.endswith('_ref')]
        # get frequency power-law indices if specified
        # NOTE these should be supplied in units of 1/3
        gammas = None
        if cp.has_option(section, 'gammas'):
            skip_args.append('gammas')
            gammas = numpy.array([float(g) / 3. for g in
                                  cp.get(section, 'gammas').split()])
        args = super(Relative, Relative).extra_args_from_config(
            cp, section, skip_args=skip_args, dtypes=dtypes)
        # get fiducial params from config
        fid_params = {p.replace('_ref', ''): float(cp.get('model', p)) for p
                      in cp.options('model') if p.endswith('_ref')}
        # add optional params with default values if not specified
        opt_params = {'ra': numpy.pi, 'dec': 0.0, 'inclination': 0.0,
                      'polarization': numpy.pi}
        fid_params.update({p: opt_params[p] for p in opt_params if p
                           not in fid_params})
        args.update({'fiducial_params': fid_params, 'gammas': gammas})
        return args
