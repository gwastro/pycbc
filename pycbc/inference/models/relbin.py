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


import numpy
import logging
from scipy.interpolate import interp1d

from pycbc.waveform.spa_tmplt import spa_tmplt
from pycbc.detector import Detector

from .base_data import BaseDataModel


def setup_bins(f_full, f_lo, f_hi, chi=1.0, eps=0.5):
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
    ga = numpy.array([-5./3, -2./3, 1., 5./3, 7./3])
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
    fbin_ind = numpy.array([numpy.argmin(numpy.absolute(f_full - ff)) for
                            ff in fbin])
    # make sure grid points are precise
    fbin = numpy.array([f_full[i] for i in fbin_ind])
    return nbin, fbin, fbin_ind

class Relative(BaseDataModel):
    r"""Model that assumes the likelihood in a region around the peak
    is slowly varying such that a linear approximation can be made, and
    likelihoods can be calculated at a coarser frequency resolution. For
    more details on the implementation, see https://arxiv.org/abs/1806.08792.

    This model requires the use of a fiducial waveform whose parameters are
    near the peak of the likelihood. The fiducial waveform and all template
    waveforms used in likelihood calculation are currently generated using
    the SPAtmplt approximant.

    For more details on initialization parameters and definition of terms, see
    :py:class:`models.BaseDataModel`.
    Parameters
    ----------
    data : dict
        A dictionary of data, in which the keys are the detector names and the
        values are the data (assumed to be unwhitened).
    psds : dict
        A dictionary of FrequencySeries keyed by the detector names. The
        dictionary must have a psd for each detector specified in the data
        dictionary. The inner products in each detector will be
        weighted by 1/psd of that detector.
    mass1 : float
        The primary mass in solar masses used for generating the fiducial
        waveform.
    mass2 : float
        The secondary mass in solar masses used for generating the fiducial
        waveform.
    spin1z : float
        The component of primary dimensionless spin along the orbital angular
        momentum used for generating the fiducial waveform.
    spin2z : float
        The component of secondary dimensionless spin along the orbital angular
        momentum used for generating the fiducial waveform.
    ra : float
        The right ascension in radians used for generating the fiducial
        waveform.
    dec : float
        The declination in radians used for generating the fiducial waveform.
    tc : float
        The GPS time of coalescence used for generating the fiducial waveform.
    low_frequency_cutoff : float
        The starting frequency used in computing inner products. This will be
        the same for all detectors.
    high_frequency_cutoff : float
        The ending frequency used in computing inner products. This will be
        the same for all detectors. Note that no checks are made to ensure
        that a template waveform does not end below this frequency, so
        choose this value to be lower than any expected final frequency
        for template waveforms.
    epsilon : float, optional
        Tuning parameter used in calculating the frequency bins. Lower values
        will result in higher resolution and more bins.
    \**kwargs :
        All other keyword arguments are passed to ``BaseDataModel``.
    """
    name = "relative"

    def __init__(self, data, psds, mass1, mass2, spin1z, spin2z,
                 ra, dec, tc, low_frequency_cutoff,
                 high_frequency_cutoff, epsilon=0.5, **kwargs):
        super(Relative, self).__init__(data=data, **kwargs)
        # store data and frequencies
        f_lo = float(low_frequency_cutoff)
        f_hi = float(high_frequency_cutoff)
        self.f = numpy.array(data[data.keys()[0]].sample_frequencies)
        self.df = data[data.keys()[0]].delta_f
        self.data = data
        self.end_time = float(data[data.keys()[0]].end_time)
        self.det = {ifo: Detector(ifo) for ifo in data}
        epsilon = float(epsilon)
        # store data and psds as arrays for faster computation
        self.comp_data = {ifo: numpy.array(data[ifo].data) for ifo in data}
        self.comp_psds = {ifo: numpy.array(psds[ifo].data) for ifo in data}
        self._psds = psds
        # store fiducial waveform params
        mass1 = float(mass1)
        mass2 = float(mass2)
        spin1z = float(spin1z)
        spin2z = float(spin2z)
        ra = float(ra)
        dec = float(dec)
        tc = float(tc)

        # get detector-specific arrival times relative to end of data
        dt = {ifo: self.det[ifo].time_delay_from_earth_center(ra, dec, tc) for
              ifo in data}
        self.ta = {ifo: tc + dt[ifo] - self.end_time for ifo in data}

        # generate fiducial waveform
        logging.info("Generating fiducial waveform")
        hp = spa_tmplt(f_lower=f_lo, f_upper=f_hi+self.df,
                       delta_f=self.df, mass1=mass1,
                       mass2=mass2, spin1z=spin1z,
                       spin2z=spin2z, distance=1.,
                       spin_order=-1, phase_order=-1)
        hp.resize(len(self.f))
        self.h00 = numpy.array(hp)

        # compute frequency bins
        logging.info("Computing frequency bins")
        nbin, fbin, fbin_ind = setup_bins(f_full=self.f, f_lo=f_lo, f_hi=f_hi,
                                          eps=epsilon)
        logging.info("Using %s bins for this model", nbin)
        # store bins and edges in sample and frequency space
        self.edges = fbin_ind
        self.fedges = numpy.array(fbin).astype(numpy.float32)
        self.bins = numpy.array([(self.edges[i], self.edges[i+1]) for
                                 i in range(len(self.edges) - 1)])
        self.fbins = numpy.array([(fbin[i], fbin[i+1]) for
                                  i in range(len(fbin) - 1)])
        # store low res copy of fiducial waveform
        self.h00_sparse = self.h00.copy().take(self.edges)

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
        # timeshift the fiducial waveform for each detector
        shift = {ifo: numpy.exp(-2.0j * numpy.pi * self.f * self.ta[ifo]) for
                 ifo in self.data}
        h0 = {ifo: self.h00.copy() * shift[ifo] for ifo in self.data}
        # calculate coefficients
        sdat = {}
        for ifo in self.data:
            hd = numpy.conjugate(self.comp_data[ifo]) * h0[ifo]
            hd /= self.comp_psds[ifo]
            hh = (numpy.absolute(h0[ifo]) ** 2.0) / self.comp_psds[ifo]
            # constant terms
            a0 = numpy.array([4. * self.df * numpy.sum(hd[l:h]) for
                              l, h in self.bins])
            b0 = numpy.array([4. * self.df * numpy.sum(hh[l:h]) for
                              l, h in self.bins])
            # linear terms
            bin_centers = [0.5 * (fl + fh) for fl, fh in self.fbins]
            a1 = numpy.array([4. * self.df * \
                              numpy.sum(hd[l:h] * (self.f[l:h] - bc)) for
                              (l, h), bc in zip(self.bins, bin_centers)])
            b1 = numpy.array([4. * self.df * \
                              numpy.sum(hh[l:h] * (self.f[l:h] - bc)) for
                              (l, h), bc in zip(self.bins, bin_centers)])

            sdat[ifo] = {'a0': a0, 'a1': a1,
                         'b0': b0, 'b1': b1}
        return sdat

    def waveform_ratio(self, p, htf, dtc=0.0):
        """Calculate waveform ratio between template and fiducial
        waveforms.
        """
        # generate template
        hp = spa_tmplt(sample_points=self.fedges,
                       mass1=p['mass1'], mass2=p['mass2'],
                       spin1z=p['spin1z'], spin2z=p['spin2z'],
                       distance=1., spin_order=-1, phase_order=-1)
        htarget = numpy.array(hp)
        # apply antenna pattern, inclination, and distance factor
        htarget *= htf
        # compute waveform ratio and timeshift
        shift = numpy.exp(-2.0j * numpy.pi * self.fedges * dtc)
        r = htarget / self.h00_sparse * shift
        r0 = 0.5 * (r[:-1] + r[1:])
        r1 = (r[1:] - r[:-1]) / (self.fedges[1:] - self.fedges[:-1])
        return numpy.array([r0, r1], dtype=numpy.complex128)

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

        llr = 0.
        for ifo in self.data:
            # get detector antenna pattern
            fp, fc = self.det[ifo].antenna_pattern(p['ra'], p['dec'],
                                                   p['polarization'],
                                                   p['tc'])
            ip = numpy.cos(p['inclination'])
            ic = 0.5 * (1.0 + ip * ip)
            htf = (fp * ip + 1.0j * fc * ic) / p['distance']
            # get timeshift relative to fiducial waveform
            dt = self.det[ifo].time_delay_from_earth_center(p['ra'], p['dec'],
                                                            p['tc'])
            dtc = p['tc'] + dt - self.end_time - self.ta[ifo]
            # generate template and calculate waveform ratio
            r0, r1 = self.waveform_ratio(p, htf, dtc=dtc)
            # <h, d> is real part of sum over bins of A0r0 + A1r1
            hd = numpy.sum(self.sdat[ifo]['a0'] * r0
                           + self.sdat[ifo]['a1'] * r1).real
            # <h, h> is real part of sum over bins of B0|r0|^2 + 2B1Re(r1r0*)
            hh = numpy.sum(self.sdat[ifo]['b0'] * numpy.absolute(r0) ** 2.
                           + 2. * self.sdat[ifo]['b1']
                           * (r1 * numpy.conjugate(r0)).real).real
            # increment loglr
            llr += (hd - 0.5 * hh)
        return float(llr)

    def _loglikelihood(self):
        r"""Computes the log likelihood of the parameters,

        .. math::

            \log p(d|\Theta, h) = \log \alpha -\frac{1}{2}\sum_i
                \left<d_i - h_i(\Theta) | d_i - h_i(\Theta)\right>,

        at the current parameter values :math:`\Theta`.

        Returns
        -------
        float
            The value of the log likelihood evaluated at the given point.
        """
        return self.loglr + self.lognl

    def _lognl(self):
        """Calculate the log of the noise likelihood. Currently not
        implemented, so just returns zero.
        """
        return 0.

    def write_metadata(self, fp):
        """Adds writing the psds and lognl, since it's a constant.

        The lognl is written to the sample group's ``attrs``.

        Parameters
        ----------
        fp : pycbc.inference.io.BaseInferenceFile instance
            The inference file to write to.
        """
        super(Relative, self).write_metadata(fp)
        if self._psds is not None:
            fp.write_psd(self._psds)
        try:
            attrs = fp[fp.samples_group].attrs
        except KeyError:
            # group doesn't exist, create it
            fp.create_group(fp.samples_group)
            attrs = fp[fp.samples_group].attrs
        attrs['lognl'] = self.lognl
